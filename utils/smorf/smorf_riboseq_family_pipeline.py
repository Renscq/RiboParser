#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.008
# Function: Evaluate family-aware smORF evidence with compatible task progress reporting.
# Input: smorf_cluster family tables, density tracks, and genePred annotation.
# Output: Unit-level, family-level, member-level, calibration, and summary tables.

"""Family-aware adaptive Ribo-seq evidence pipeline."""

from __future__ import annotations

import csv
import gzip
import math
import multiprocessing as mp
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import Final, Iterator, Mapping, Sequence, TextIO

import numpy as np
import pandas as pd

from utils.ribo.ArgsParser import message_print, progress_print, warning_print

from .smorf_riboseq_adaptive import (
    AdaptiveEvidenceConfig,
    AdaptiveEvidenceResult,
    EVIDENCE_LEVELS,
    evaluate_adaptive_evidence,
    evidence_max,
)
from .smorf_riboseq_constants import EvidenceThresholds
from .smorf_riboseq_density import (
    ChromDensity,
    PreparedDensity,
    iter_density_chromosomes,
    prepare_density,
)
from .smorf_riboseq_density_cache import (
    DensityCacheSet,
    available_memory_bytes,
    build_density_cache,
    chromosome_cache_bytes,
    load_chromosome_map,
)
from .smorf_riboseq_family_io import (
    DensitySample,
    FamilyInput,
    read_density_design,
    read_family_input,
    resolve_analysis_mode,
)
from .smorf_riboseq_io import GenePredRecord, read_genepred
from .smorf_riboseq_metrics import calculate_release
from .smorf_riboseq_profile import (
    extract_transcript_downstream_profile,
    extract_transcript_profile,
    get_coding_nt_length,
    make_codon_profile,
)

UNIT_COLUMNS: Final[tuple[str, ...]] = (
    "evidence_scope",
    "evidence_unit",
    "group",
    "family_id",
    "structural_primary",
    "evidence_primary",
    "start_resolution",
    "family_translation_evidence",
    "family_evidence_reason",
    "family_common_body_orf",
    "family_common_body_codons",
    "representative_count",
    "supported_representative_count",
    "rpf_sum",
    "rpf_per_codon",
    "covered_codon",
    "covered_codon_ratio",
    "frame0_density",
    "frame1_density",
    "frame2_density",
    "frame0_ratio",
    "length_class",
    "window_count",
    "supported_window_count",
    "separated_supported_window_count",
    "supported_window_fraction",
    "best_window_rpf",
    "best_window_frame0_ratio",
    "top3_window_score_median",
    "coverage_span_ratio",
    "window_coverage_span_ratio",
    "signal_span_ratio",
    "signal_bin_count",
    "top_window_rpf_fraction",
    "start_support",
    "body_support",
    "release_support",
    "distributed_support",
    "distributed_sparse_support",
    "boundary_support",
    "high_depth_support",
    "localized_only",
    "evidence_evaluability",
)

UNIT_MEMBER_COLUMNS: Final[tuple[str, ...]] = (
    "evidence_scope",
    "evidence_unit",
    "group",
    "family_id",
    "representative_orf_id",
    "structural_primary",
    "evidence_primary",
    "start_resolution",
    "member_category",
    "member_start_codon",
    "member_aa_length",
    "member_translation_evidence",
    "member_status",
    "evidence_reason",
    "extension_codons",
    "extension_rpf",
    "extension_covered_codon",
    "extension_frame0_ratio",
    "extension_support",
    "start_support",
    "localized_only",
    "coverage_span_ratio",
)

FAMILY_COLUMNS: Final[tuple[str, ...]] = (
    "family_id",
    "gene_id",
    "category",
    "structural_primary",
    "evidence_primary",
    "family_translation_evidence",
    "start_resolution",
    "consensus_scope",
    "consensus_unit",
    "supported_sample_count",
    "supported_group_count",
    "high_sample_count",
    "medium_sample_count",
    "max_sample_evidence",
    "max_group_evidence",
    "pooled_evidence",
    "family_common_body_orf",
    "family_common_body_codons",
    "coverage_span_ratio",
    "supported_window_count",
    "localized_only",
    "selection_reason",
)

MEMBER_COLUMNS: Final[tuple[str, ...]] = (
    "family_id",
    "member_orf_id",
    "representative_orf_id",
    "structural_primary",
    "evidence_primary",
    "member_category",
    "member_start_codon",
    "member_aa_length",
    "family_role",
    "collapse_reason",
    "evidence_inherited_from",
    "member_translation_evidence",
    "member_status",
    "extension_codons",
    "extension_rpf",
    "extension_covered_codon",
    "extension_frame0_ratio",
    "extension_support",
    "start_support",
    "member_selection_reason",
)

CALIBRATION_COLUMNS: Final[tuple[str, ...]] = (
    "evidence_scope",
    "evidence_unit",
    "group",
    "positive_control_count",
    "positive_control_eligible",
    "calibration_status",
    "base_moderate_periodicity",
    "base_strong_periodicity",
    "positive_frame0_q20",
    "positive_frame0_median",
    "positive_rpf_per_codon_q20",
    "calibrated_moderate_periodicity",
    "calibrated_strong_periodicity",
    "calibrated_min_window_rpf",
)

_BUFFER_BYTES: Final[int] = 8 * 1024 * 1024


@dataclass(frozen=True, slots=True)
class AnalysisUnit:
    """Describe one sample, group-pooled, or globally pooled evidence unit."""

    scope: str
    name: str
    group: str
    samples: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PreparedTrack:
    """Store one prepared density track reusable across multiple passes."""

    sample: str
    strand: str
    file_format: str
    prepared: PreparedDensity


@dataclass(slots=True)
class ProfileScore:
    """Store one representative profile and its adaptive evidence."""

    row: pd.Series
    nt_profile: np.ndarray
    codon_profile: np.ndarray
    downstream_profile: np.ndarray
    release: dict[str, object]
    adaptive: AdaptiveEvidenceResult


@dataclass(frozen=True, slots=True)
class EvidencePipelineResult:
    """Store output paths and analysis counts."""

    mode: str
    sample_count: int
    family_count: int
    representative_count: int
    outputs: tuple[str, ...]


class AtomicTableWriter:
    """Write one TSV or TSV.GZ file atomically."""

    def __init__(self, path: str | Path, columns: Sequence[str]) -> None:
        self.final_path = Path(path)
        self.columns = tuple(columns)
        self.temp_path = self.final_path.with_name(self.final_path.name + ".tmp")
        self.handle: TextIO | None = None
        self.writer: csv.DictWriter | None = None

    def __enter__(self) -> "AtomicTableWriter":
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        compressed = self.final_path.name.lower().endswith(".gz")
        if compressed:
            self.handle = gzip.open(
                self.temp_path,
                "wt",
                encoding="utf-8",
                newline="",
                compresslevel=1,
            )
        else:
            self.handle = self.temp_path.open(
                "w",
                encoding="utf-8",
                newline="",
                buffering=_BUFFER_BYTES,
            )
        self.writer = csv.DictWriter(
            self.handle,
            fieldnames=self.columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        self.writer.writeheader()
        return self

    def write(self, row: Mapping[str, object]) -> None:
        """Write one row in fixed column order."""
        if self.writer is None:
            raise RuntimeError("Writer is not active.")
        self.writer.writerow({column: row.get(column, "") for column in self.columns})

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        if self.handle is not None:
            self.handle.close()
        self.handle = None
        self.writer = None
        if exception_type is None:
            self.temp_path.replace(self.final_path)
        else:
            self.temp_path.unlink(missing_ok=True)
        return False


class PreparedDensitySet:
    """Prepare density tracks once and reopen iterators for each pass."""

    def __init__(
        self,
        samples: Sequence[DensitySample],
        work_directory: str | Path,
    ) -> None:
        self.samples = tuple(samples)
        self.work_directory = Path(work_directory)
        self.tracks: tuple[PreparedTrack, ...] = ()

    def __enter__(self) -> "PreparedDensitySet":
        prepared_tracks: list[PreparedTrack] = []
        for sample in self.samples:
            for track in sample.tracks:
                prepared = prepare_density(
                    path=track.path,
                    file_format=track.file_format,
                    work_directory=self.work_directory,
                )
                prepared_tracks.append(
                    PreparedTrack(
                        sample=sample.sample,
                        strand=track.strand,
                        file_format=track.file_format,
                        prepared=prepared,
                    )
                )
        self.tracks = tuple(prepared_tracks)
        return self

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        for track in self.tracks:
            track.prepared.cleanup()
        return False


class DensityCursor:
    """Read one prepared density file chromosome by chromosome."""

    def __init__(self, track: PreparedTrack) -> None:
        self.track = track
        self.iterator = iter_density_chromosomes(
            track.prepared.path,
            track.file_format,
        )
        self.current: ChromDensity | None = None
        self.exhausted = False
        self.previous_chrom: str | None = None
        self._advance()

    def _advance(self) -> None:
        if self.exhausted:
            return
        try:
            next_density = next(self.iterator)
        except StopIteration:
            self.current = None
            self.exhausted = True
            return
        if self.previous_chrom is not None and next_density.chrom < self.previous_chrom:
            raise ValueError(
                "Density chromosomes are not lexicographically ordered: "
                f"{self.track.prepared.path}. Convert WIG to ordered bedGraph "
                "or sort chromosome blocks before hybrid analysis."
            )
        self.previous_chrom = next_density.chrom
        self.current = next_density

    def get(self, chrom: str) -> ChromDensity | None:
        """Return one chromosome density or ``None`` when absent."""
        while self.current is not None and self.current.chrom < chrom:
            self._advance()
        if self.current is not None and self.current.chrom == chrom:
            return self.current
        return None


class DensityPass:
    """Create synchronized chromosome cursors for all prepared tracks."""

    def __init__(self, prepared_set: PreparedDensitySet) -> None:
        self.cursors = {
            (track.sample, track.strand): DensityCursor(track)
            for track in prepared_set.tracks
        }

    def chromosome_map(
        self,
        chrom: str,
    ) -> dict[tuple[str, str], ChromDensity | None]:
        """Return all sample/strand densities for one chromosome."""
        return {
            key: cursor.get(chrom)
            for key, cursor in self.cursors.items()
        }


def _build_base_config(args: object) -> AdaptiveEvidenceConfig:
    """Build adaptive parameters from old evidence presets and new options."""
    old = EvidenceThresholds.from_mode(
        str(getattr(args, "evidence_mode", "balanced"))
    ).with_overrides(
        min_rpf_sum=getattr(args, "min_rpf_sum", None),
        min_covered_codon=getattr(args, "min_covered_codon", None),
        min_codon_coverage=getattr(args, "min_codon_coverage", None),
        moderate_periodicity=getattr(args, "moderate_periodicity", None),
        strong_periodicity=getattr(args, "strong_periodicity", None),
    )
    config = AdaptiveEvidenceConfig(
        min_rpf_sum=old.min_rpf_sum,
        min_covered_codon=old.min_covered_codon,
        min_codon_coverage=old.min_codon_coverage,
        moderate_periodicity=old.moderate_periodicity,
        strong_periodicity=old.strong_periodicity,
        short_max_codons=int(getattr(args, "short_max_codons", 30)),
        long_min_codons=int(getattr(args, "long_min_codons", 100)),
        window_codons=int(getattr(args, "window_codons", 20)),
        window_step_codons=int(getattr(args, "window_step_codons", 5)),
        min_window_rpf=float(getattr(args, "min_window_rpf", 3.0)),
        min_window_covered_codons=int(
            getattr(args, "min_window_covered_codons", 3)
        ),
        min_supported_windows=int(getattr(args, "min_supported_windows", 2)),
        min_window_gap_codons=int(
            getattr(args, "min_window_gap_codons", 15)
        ),
        min_coverage_span=float(getattr(args, "min_coverage_span", 0.35)),
        localized_span_max=float(getattr(args, "localized_span_max", 0.20)),
        localized_top_window_fraction=float(
            getattr(args, "localized_top_window_fraction", 0.70)
        ),
        positive_quantile=float(getattr(args, "positive_quantile", 0.20)),
        positive_min_controls=int(
            getattr(args, "positive_min_controls", 30)
        ),
        positive_min_rpf_sum=float(
            getattr(args, "positive_min_rpf_sum", 10.0)
        ),
        pseudocount=old.pseudocount,
    )
    config.validate()
    return config


def _analysis_units(
    samples: Sequence[DensitySample],
    mode: str,
) -> tuple[AnalysisUnit, ...]:
    """Build deterministic analysis units for the selected mode."""
    units: list[AnalysisUnit] = []
    sample_names = tuple(sample.sample for sample in samples)
    group_members: dict[str, list[str]] = defaultdict(list)
    for sample in samples:
        group_members[sample.group].append(sample.sample)

    if mode in {"single", "separate", "hybrid"}:
        for sample in samples:
            units.append(
                AnalysisUnit(
                    scope="sample",
                    name=sample.sample,
                    group=sample.group,
                    samples=(sample.sample,),
                )
            )
    if mode == "hybrid":
        for group, members in group_members.items():
            if len(members) < 2:
                continue
            units.append(
                AnalysisUnit(
                    scope="group",
                    name=group,
                    group=group,
                    samples=tuple(members),
                )
            )
    if mode in {"pooled", "hybrid"}:
        units.append(
            AnalysisUnit(
                scope="pooled",
                name="all",
                group="all",
                samples=sample_names,
            )
        )
    return tuple(units)


def _density_for_sample(
    chromosome_densities: Mapping[tuple[str, str], ChromDensity | None],
    sample: str,
    strand: str,
) -> ChromDensity | None:
    """Select one sample's unstranded or matching strand density."""
    unstranded_key = (sample, ".")
    if unstranded_key in chromosome_densities:
        return chromosome_densities[unstranded_key]
    return chromosome_densities.get((sample, strand))


def _sample_profiles(
    row: pd.Series,
    samples: Sequence[DensitySample],
    chromosome_densities: Mapping[tuple[str, str], ChromDensity | None],
    genepred_records: Mapping[str, GenePredRecord],
    post_stop_codons: int,
) -> dict[str, tuple[np.ndarray, np.ndarray, str]]:
    """Extract ORF and downstream profiles for all samples."""
    starts = list(row["exon_starts_parsed"])
    ends = list(row["exon_ends_parsed"])
    strand = str(row["strand"])
    transcript = genepred_records.get(str(row["transcript_id"]))
    profiles: dict[str, tuple[np.ndarray, np.ndarray, str]] = {}
    for sample in samples:
        density = _density_for_sample(
            chromosome_densities,
            sample.sample,
            strand,
        )
        nt_profile = extract_transcript_profile(
            density=density,
            starts=starts,
            ends=ends,
            strand=strand,
            orf_id=str(row["orf_id"]),
        )
        if len(nt_profile) != int(row["nt_length"]):
            raise ValueError(
                f"ORF {row['orf_id']} extracted length {len(nt_profile)} "
                f"does not match nt_length {row['nt_length']}."
            )
        downstream, context = extract_transcript_downstream_profile(
            density=density,
            orf_starts=starts,
            orf_ends=ends,
            strand=strand,
            transcript=transcript,
            nt_window=post_stop_codons * 3,
            orf_id=str(row["orf_id"]),
        )
        profiles[sample.sample] = (nt_profile, downstream, context)
    return profiles


def _unit_profile(
    sample_profiles: Mapping[str, tuple[np.ndarray, np.ndarray, str]],
    unit: AnalysisUnit,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Sum sample profiles for one sample/group/pooled unit."""
    nt_profiles = [sample_profiles[sample][0] for sample in unit.samples]
    downstream_profiles = [sample_profiles[sample][1] for sample in unit.samples]
    if not nt_profiles:
        return (
            np.zeros(0, dtype=np.float32),
            np.zeros(0, dtype=np.float32),
            "unavailable",
        )
    nt_profile = np.sum(np.stack(nt_profiles), axis=0, dtype=np.float64)
    available_downstream = [profile for profile in downstream_profiles if len(profile)]
    if available_downstream:
        max_length = max(len(profile) for profile in available_downstream)
        downstream = np.zeros(max_length, dtype=np.float64)
        for profile in available_downstream:
            downstream[: len(profile)] += profile
        context = "pooled_transcript" if len(unit.samples) > 1 else "transcript"
    else:
        downstream = np.zeros(0, dtype=np.float64)
        context = "unavailable"
    return nt_profile, downstream, context


def _score_profile(
    row: pd.Series,
    nt_profile: np.ndarray,
    downstream: np.ndarray,
    release_context: str,
    config: AdaptiveEvidenceConfig,
    post_stop_codons: int,
) -> ProfileScore:
    """Score one representative ORF profile."""
    coding_nt_length = get_coding_nt_length(row, len(nt_profile))
    coding_nt = np.asarray(nt_profile[:coding_nt_length], dtype=np.float64)
    codon_profile = make_codon_profile(nt_profile, coding_nt_length)
    release_thresholds = EvidenceThresholds(
        min_rpf_sum=config.min_rpf_sum,
        min_covered_codon=config.min_covered_codon,
        min_codon_coverage=config.min_codon_coverage,
        strong_periodicity=config.strong_periodicity,
        moderate_periodicity=config.moderate_periodicity,
        pseudocount=config.pseudocount,
    )
    release = calculate_release(
        codon_profile=codon_profile,
        downstream_nt_profile=downstream,
        post_stop_codons=post_stop_codons,
        thresholds=release_thresholds,
        context=release_context,
    )
    adaptive = evaluate_adaptive_evidence(
        nt_profile=coding_nt,
        codon_profile=codon_profile,
        category=str(row.get("category", ".")),
        config=config,
        release_label=str(release.get("release_label", "NA")),
        release_evaluable=bool(release.get("release_evaluable", False)),
    )
    return ProfileScore(
        row=row,
        nt_profile=coding_nt,
        codon_profile=codon_profile,
        downstream_profile=downstream,
        release=release,
        adaptive=adaptive,
    )


def _representative_order(records: pd.DataFrame) -> pd.DataFrame:
    """Sort family representatives from longest to shortest ORF."""
    output = records.copy()
    aa = pd.to_numeric(output.get("aa_length", 0), errors="coerce").fillna(
        (output["nt_length"] - 3) // 3
    )
    output["_aa_sort"] = aa.astype(int)
    output["_structural_sort"] = output["structural_primary"].astype(int)
    return output.sort_values(
        ["_aa_sort", "_structural_sort", "orf_id"],
        ascending=[False, False, True],
        kind="mergesort",
    )


def _extension_metrics(
    longer: ProfileScore,
    next_shorter_codons: int,
    config: AdaptiveEvidenceConfig,
) -> dict[str, object]:
    """Score the member-specific N-terminal extension of one representative."""
    total_codons = len(longer.codon_profile)
    extension_codons = max(0, total_codons - next_shorter_codons)
    if extension_codons <= 0:
        return {
            "extension_codons": 0,
            "extension_rpf": 0.0,
            "extension_covered_codon": 0,
            "extension_frame0_ratio": 0.0,
            "extension_support": False,
        }
    nt_values = longer.nt_profile[: extension_codons * 3]
    codon_values = longer.codon_profile[:extension_codons]
    frame0 = float(nt_values[0::3].sum())
    frame1 = float(nt_values[1::3].sum())
    frame2 = float(nt_values[2::3].sum())
    total = frame0 + frame1 + frame2
    ratio = frame0 / total if total > 0 else 0.0
    rpf_sum = float(nt_values.sum())
    covered = int(np.count_nonzero(codon_values > 0))
    support = (
        rpf_sum >= max(2.0, config.min_window_rpf * 0.5)
        and covered >= min(2, extension_codons)
        and ratio >= config.moderate_periodicity
    )
    return {
        "extension_codons": extension_codons,
        "extension_rpf": rpf_sum,
        "extension_covered_codon": covered,
        "extension_frame0_ratio": ratio,
        "extension_support": support,
    }


def _select_family_start(
    family_records: pd.DataFrame,
    scores: Mapping[str, ProfileScore],
    config: AdaptiveEvidenceConfig,
) -> tuple[str, str, dict[str, dict[str, object]], str]:
    """Choose an evidence-supported start without redefining family structure."""
    ordered = _representative_order(family_records)
    representative_ids = ordered["orf_id"].tolist()
    extension_by_orf: dict[str, dict[str, object]] = {}
    for index, orf_id in enumerate(representative_ids):
        next_codons = (
            len(scores[representative_ids[index + 1]].codon_profile)
            if index + 1 < len(representative_ids)
            else len(scores[orf_id].codon_profile)
        )
        extension_by_orf[orf_id] = _extension_metrics(
            scores[orf_id],
            next_codons,
            config,
        )

    structural_primary = str(
        family_records.loc[family_records["structural_primary"], "orf_id"].iloc[0]
    )
    if len(representative_ids) == 1:
        return structural_primary, "Singleton", extension_by_orf, "singleton_family"

    annotated = [
        orf_id
        for orf_id in representative_ids
        if str(scores[orf_id].row.get("category", "")).lower()
        == "annotated_orf"
        and EVIDENCE_LEVELS.get(
            str(scores[orf_id].adaptive.metrics["translation_evidence"]),
            0,
        )
        >= 1
    ]
    if annotated:
        return annotated[0], "ResolvedAnnotated", extension_by_orf, "annotated_mORF_supported"

    extension_supported = [
        orf_id
        for orf_id in representative_ids
        if bool(extension_by_orf[orf_id]["extension_support"])
    ]
    if extension_supported:
        return (
            extension_supported[0],
            "ResolvedExtension",
            extension_by_orf,
            "upstream_extension_supported",
        )
    return (
        structural_primary,
        "AmbiguousStart",
        extension_by_orf,
        "shared_body_without_start_resolution",
    )


def _family_unit_result(
    unit: AnalysisUnit,
    family_records: pd.DataFrame,
    scores: Mapping[str, ProfileScore],
    config: AdaptiveEvidenceConfig,
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    """Build one family result for one evidence unit."""
    ordered = _representative_order(family_records)
    shortest = ordered.iloc[-1]
    common_orf = str(shortest["orf_id"])
    common_score = scores[common_orf]
    evidence_primary, resolution, extensions, selection_reason = _select_family_start(
        family_records,
        scores,
        config,
    )
    common_metrics = common_score.adaptive.metrics
    representative_evidence = [
        str(score.adaptive.metrics["translation_evidence"])
        for score in scores.values()
    ]
    family_evidence = str(common_metrics["translation_evidence"])
    if EVIDENCE_LEVELS.get(evidence_max(representative_evidence), 0) > EVIDENCE_LEVELS.get(
        family_evidence,
        0,
    ):
        family_evidence = evidence_max(representative_evidence)
        family_reason = "representative_specific_support"
    else:
        family_reason = str(common_metrics["evidence_reason"])

    structural_primary = str(
        family_records.loc[family_records["structural_primary"], "orf_id"].iloc[0]
    )
    row: dict[str, object] = {
        "evidence_scope": unit.scope,
        "evidence_unit": unit.name,
        "group": unit.group,
        "family_id": str(family_records["family_id"].iloc[0]),
        "structural_primary": structural_primary,
        "evidence_primary": evidence_primary,
        "start_resolution": resolution,
        "family_translation_evidence": family_evidence,
        "family_evidence_reason": family_reason + ";" + selection_reason,
        "family_common_body_orf": common_orf,
        "family_common_body_codons": len(common_score.codon_profile),
        "representative_count": len(family_records),
        "supported_representative_count": sum(
            EVIDENCE_LEVELS.get(label, 0) >= 1
            for label in representative_evidence
        ),
    }
    row.update(common_metrics)
    return row, extensions


def _unit_member_rows(
    unit: AnalysisUnit,
    family_records: pd.DataFrame,
    unit_result: Mapping[str, object],
    scores: Mapping[str, ProfileScore],
    extensions: Mapping[str, Mapping[str, object]],
) -> Iterator[dict[str, object]]:
    """Yield representative-level evidence for one sample/group/pooled unit."""
    evidence_primary = str(unit_result["evidence_primary"])
    structural_primary = str(unit_result["structural_primary"])
    start_resolution = str(unit_result["start_resolution"])
    for _, record in _representative_order(family_records).iterrows():
        representative_id = str(record["orf_id"])
        score = scores[representative_id]
        extension = extensions[representative_id]
        label = str(score.adaptive.metrics["translation_evidence"])
        if representative_id == evidence_primary:
            status = "EvidencePrimary"
        elif bool(extension["extension_support"]):
            status = "AlternativeStartSupported"
        elif start_resolution == "AmbiguousStart":
            status = "AmbiguousStart"
        else:
            status = "CollapsedUnsupported"
        yield {
            "evidence_scope": unit.scope,
            "evidence_unit": unit.name,
            "group": unit.group,
            "family_id": unit_result["family_id"],
            "representative_orf_id": representative_id,
            "structural_primary": structural_primary,
            "evidence_primary": evidence_primary,
            "start_resolution": start_resolution,
            "member_category": record.get("category", "."),
            "member_start_codon": record.get("start_codon", "."),
            "member_aa_length": record.get("aa_length", ""),
            "member_translation_evidence": label,
            "member_status": status,
            "evidence_reason": score.adaptive.metrics["evidence_reason"],
            "extension_codons": extension["extension_codons"],
            "extension_rpf": extension["extension_rpf"],
            "extension_covered_codon": extension["extension_covered_codon"],
            "extension_frame0_ratio": extension["extension_frame0_ratio"],
            "extension_support": extension["extension_support"],
            "start_support": score.adaptive.metrics["start_support"],
            "localized_only": score.adaptive.metrics["localized_only"],
            "coverage_span_ratio": score.adaptive.metrics[
                "coverage_span_ratio"
            ],
        }


def _consensus_family_result(
    family_records: pd.DataFrame,
    unit_rows: Sequence[dict[str, object]],
    mode: str,
) -> dict[str, object]:
    """Build one project-level family consensus from sample/group/pooled units."""
    sample_rows = [row for row in unit_rows if row["evidence_scope"] == "sample"]
    group_rows = [row for row in unit_rows if row["evidence_scope"] == "group"]
    pooled_rows = [row for row in unit_rows if row["evidence_scope"] == "pooled"]

    supported_samples = [
        row
        for row in sample_rows
        if EVIDENCE_LEVELS.get(str(row["family_translation_evidence"]), 0) >= 1
    ]
    supported_groups = [
        row
        for row in group_rows
        if EVIDENCE_LEVELS.get(str(row["family_translation_evidence"]), 0) >= 1
    ]
    max_sample = evidence_max(
        str(row["family_translation_evidence"]) for row in sample_rows
    )
    max_group = evidence_max(
        str(row["family_translation_evidence"]) for row in group_rows
    )
    pooled_evidence = (
        str(pooled_rows[0]["family_translation_evidence"])
        if pooled_rows
        else "NoEvidence"
    )

    if mode == "pooled":
        selected = pooled_rows[0]
        reason = "pooled_density_consensus"
    elif mode == "single":
        selected = sample_rows[0]
        reason = "single_sample_consensus"
    elif mode == "separate":
        selected = max(
            sample_rows,
            key=lambda row: (
                EVIDENCE_LEVELS.get(str(row["family_translation_evidence"]), 0),
                float(row.get("coverage_span_ratio", 0.0)),
                float(row.get("rpf_sum", 0.0)),
            ),
        )
        reason = "best_sample_consensus"
    else:
        candidates = group_rows + pooled_rows + sample_rows
        selected = max(
            candidates,
            key=lambda row: (
                EVIDENCE_LEVELS.get(str(row["family_translation_evidence"]), 0),
                2 if row["evidence_scope"] == "group" else 1,
                len(supported_samples),
                float(row.get("coverage_span_ratio", 0.0)),
                float(row.get("rpf_sum", 0.0)),
            ),
        )
        if selected["evidence_scope"] in {"group", "pooled"} and not supported_samples:
            level = EVIDENCE_LEVELS.get(
                str(selected["family_translation_evidence"]),
                0,
            )
            if level > 1:
                selected = dict(selected)
                selected["family_translation_evidence"] = "LowConfidence"
                selected["family_evidence_reason"] = (
                    str(selected["family_evidence_reason"])
                    + ";pooled_only_without_sample_support"
                )
        reason = "hybrid_group_sample_consensus"

    structural_primary_row = family_records.loc[
        family_records["structural_primary"]
    ].iloc[0]
    return {
        "family_id": selected["family_id"],
        "gene_id": structural_primary_row["gene_id"],
        "category": structural_primary_row.get("category", "."),
        "structural_primary": selected["structural_primary"],
        "evidence_primary": selected["evidence_primary"],
        "family_translation_evidence": selected["family_translation_evidence"],
        "start_resolution": selected["start_resolution"],
        "consensus_scope": selected["evidence_scope"],
        "consensus_unit": selected["evidence_unit"],
        "supported_sample_count": len(supported_samples),
        "supported_group_count": len(supported_groups),
        "high_sample_count": sum(
            row["family_translation_evidence"] == "HighConfidence"
            for row in sample_rows
        ),
        "medium_sample_count": sum(
            row["family_translation_evidence"] == "MediumConfidence"
            for row in sample_rows
        ),
        "max_sample_evidence": max_sample,
        "max_group_evidence": max_group,
        "pooled_evidence": pooled_evidence,
        "family_common_body_orf": selected["family_common_body_orf"],
        "family_common_body_codons": selected["family_common_body_codons"],
        "coverage_span_ratio": selected["coverage_span_ratio"],
        "supported_window_count": selected["supported_window_count"],
        "localized_only": selected["localized_only"],
        "selection_reason": reason + ";" + str(selected["family_evidence_reason"]),
    }


def _member_consensus_rows(
    family_records: pd.DataFrame,
    family_members: pd.DataFrame,
    selected_unit: dict[str, object],
    selected_scores: Mapping[str, ProfileScore],
    extensions: Mapping[str, Mapping[str, object]],
) -> Iterator[dict[str, object]]:
    """Yield compact evidence rows for all original family members."""
    evidence_primary = str(selected_unit["evidence_primary"])
    structural_primary = str(selected_unit["structural_primary"])
    representative_status: dict[str, tuple[str, str]] = {}
    for representative_id, score in selected_scores.items():
        label = str(score.adaptive.metrics["translation_evidence"])
        extension_support = bool(extensions[representative_id]["extension_support"])
        if representative_id == evidence_primary:
            status = "EvidencePrimary"
        elif extension_support:
            status = "AlternativeStartSupported"
        elif selected_unit["start_resolution"] == "AmbiguousStart":
            status = "AmbiguousStart"
        else:
            status = "CollapsedUnsupported"
        representative_status[representative_id] = (label, status)

    for member in family_members.to_dict("records"):
        representative_id = str(member["representative_orf_id"])
        label, status = representative_status[representative_id]
        extension = extensions[representative_id]
        score = selected_scores[representative_id]
        yield {
            "family_id": member["family_id"],
            "member_orf_id": member["member_orf_id"],
            "representative_orf_id": representative_id,
            "structural_primary": structural_primary,
            "evidence_primary": evidence_primary,
            "member_category": member.get(
                "member_category",
                score.row.get("category", "."),
            ),
            "member_start_codon": member.get(
                "member_start_codon",
                score.row.get("start_codon", "."),
            ),
            "member_aa_length": member.get(
                "member_aa_length",
                score.row.get("aa_length", ""),
            ),
            "family_role": member.get("family_role", "."),
            "collapse_reason": member.get("collapse_reason", "."),
            "evidence_inherited_from": representative_id,
            "member_translation_evidence": label,
            "member_status": status,
            "extension_codons": extension["extension_codons"],
            "extension_rpf": extension["extension_rpf"],
            "extension_covered_codon": extension["extension_covered_codon"],
            "extension_frame0_ratio": extension["extension_frame0_ratio"],
            "extension_support": extension["extension_support"],
            "start_support": score.adaptive.metrics["start_support"],
            "member_selection_reason": selected_unit["family_evidence_reason"],
        }


def _iter_chromosome_families(
    representatives: pd.DataFrame,
) -> Iterator[tuple[str, list[tuple[str, pd.DataFrame]]]]:
    """Yield chromosome families with linear-time pandas grouping."""
    table = representatives.copy()
    table["_input_order"] = np.arange(len(table), dtype=np.int64)
    for chrom in sorted(table["chrom"].astype(str).unique()):
        chrom_table = table.loc[table["chrom"].eq(chrom)]
        families = [
            (str(family_id), family.copy())
            for family_id, family in chrom_table.groupby(
                "family_id",
                sort=False,
                observed=True,
            )
        ]
        yield str(chrom), families


def _collect_positive_controls(
    family_input: FamilyInput,
    samples: Sequence[DensitySample],
    units: Sequence[AnalysisUnit],
    prepared_set: PreparedDensitySet,
    genepred_records: Mapping[str, GenePredRecord],
    base_config: AdaptiveEvidenceConfig,
    post_stop_codons: int,
    positive_category: str,
) -> dict[tuple[str, str], list[dict[str, object]]]:
    """Collect annotated-mORF metrics for per-unit calibration."""
    controls = family_input.representatives.loc[
        family_input.representatives["category"].str.lower().eq(
            positive_category.lower()
        )
    ].copy()
    output: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    if controls.empty:
        return output
    density_pass = DensityPass(prepared_set)
    for chrom, family_tables in _iter_chromosome_families(controls):
        chromosome_densities = density_pass.chromosome_map(chrom)
        for _family_id, records in family_tables:
            for _, row in records.iterrows():
                sample_profiles = _sample_profiles(
                    row=row,
                    samples=samples,
                    chromosome_densities=chromosome_densities,
                    genepred_records=genepred_records,
                    post_stop_codons=post_stop_codons,
                )
                for unit in units:
                    nt_profile, downstream, context = _unit_profile(
                        sample_profiles,
                        unit,
                    )
                    score = _score_profile(
                        row=row,
                        nt_profile=nt_profile,
                        downstream=downstream,
                        release_context=context,
                        config=base_config,
                        post_stop_codons=post_stop_codons,
                    )
                    output[(unit.scope, unit.name)].append(score.adaptive.metrics)
    return output


def _write_summary(
    path: str | Path,
    values: Mapping[str, object],
    evidence_counts: Counter[str],
) -> None:
    """Write a simple key-value summary table atomically."""
    final = Path(path)
    temporary = final.with_name(final.name + ".tmp")
    final.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in values.items():
            handle.write(f"{key}\t{value}\n")
        for label in sorted(EVIDENCE_LEVELS, key=EVIDENCE_LEVELS.get):
            handle.write(f"family_{label}\t{evidence_counts[label]}\n")
    temporary.replace(final)




@dataclass(frozen=True, slots=True)
class ChromosomeShardResult:
    """Store one completed chromosome shard set."""

    order: int
    chrom: str
    family_count: int
    evidence_counts: Mapping[str, int]
    shards: Mapping[str, str]


class ShardTableWriter:
    """Write header-free buffered TSV shards for later ordered merging."""

    def __init__(self, path: str | Path, columns: Sequence[str]) -> None:
        self.path = Path(path)
        self.columns = tuple(columns)
        self.handle: TextIO | None = None
        self.writer: csv.DictWriter | None = None

    def __enter__(self) -> "ShardTableWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open(
            "w",
            encoding="utf-8",
            newline="",
            buffering=_BUFFER_BYTES,
        )
        self.writer = csv.DictWriter(
            self.handle,
            fieldnames=self.columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        return self

    def write(self, row: Mapping[str, object]) -> None:
        """Write one row without reconstructing unrelated fields."""
        if self.writer is None:
            raise RuntimeError("Shard writer is not active.")
        self.writer.writerow(row)

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        if self.handle is not None:
            self.handle.close()
        self.handle = None
        self.writer = None
        if exception_type is not None:
            self.path.unlink(missing_ok=True)
        return False


_WORKER_CACHE: DensityCacheSet | None = None
_WORKER_REPRESENTATIVES: pd.DataFrame | None = None
_WORKER_CHROM_INDICES: dict[str, np.ndarray] = {}
_WORKER_CONTROL_TABLE: pd.DataFrame | None = None
_WORKER_CONTROL_INDICES: dict[str, np.ndarray] = {}
_WORKER_MEMBERS_INDEXED: pd.DataFrame | None = None
_WORKER_SAMPLES: tuple[DensitySample, ...] = ()
_WORKER_UNITS: tuple[AnalysisUnit, ...] = ()
_WORKER_CONFIGS: dict[tuple[str, str], AdaptiveEvidenceConfig] = {}
_WORKER_BASE_CONFIG: AdaptiveEvidenceConfig | None = None
_WORKER_GENEPRED: Mapping[str, GenePredRecord] = {}
_WORKER_POST_STOP_CODONS = 10
_WORKER_MODE = "single"
_WORKER_SHARD_DIRECTORY = ""
_WORKER_UNIT_MEMBER_MODE = "family"


def _fork_context() -> mp.context.BaseContext | None:
    """Return the POSIX fork context required for copy-on-write tables."""
    try:
        return mp.get_context("fork")
    except ValueError:
        return None


def _available_family_members(family_id: str) -> pd.DataFrame:
    """Return indexed member rows for one family without full-table scans."""
    if _WORKER_MEMBERS_INDEXED is None:
        raise RuntimeError("Family member index is unavailable.")
    try:
        rows = _WORKER_MEMBERS_INDEXED.loc[[family_id]]
    except KeyError:
        return _WORKER_MEMBERS_INDEXED.iloc[0:0].copy()
    return rows.reset_index(drop=True)


def _select_positive_controls(
    representatives: pd.DataFrame,
    category: str,
    maximum: int,
) -> pd.DataFrame:
    """Select a deterministic chromosome-distributed control subset."""
    controls = representatives.loc[
        representatives["category"].astype(str).str.lower().eq(
            str(category).lower()
        )
    ].copy()
    if maximum <= 0 or len(controls) <= maximum:
        return controls

    controls["_control_hash"] = pd.util.hash_pandas_object(
        controls["orf_id"].astype(str),
        index=False,
    ).astype("uint64")
    selected = (
        controls.sort_values("_control_hash", kind="mergesort")
        .head(maximum)
        .drop(columns="_control_hash")
    )
    return selected


def _build_chrom_indices(table: pd.DataFrame) -> dict[str, np.ndarray]:
    """Build chromosome row-index arrays once in the parent process."""
    output: dict[str, np.ndarray] = {}
    chrom_values = table["chrom"].astype(str).to_numpy()
    for chrom in sorted(set(chrom_values)):
        output[str(chrom)] = np.flatnonzero(chrom_values == chrom)
    return output


def _calibration_chromosome_worker(
    payload: tuple[int, str],
) -> tuple[int, str, dict[tuple[str, str], list[dict[str, object]]]]:
    """Collect positive-control metrics for one cached chromosome."""
    order, chrom = payload
    if _WORKER_CACHE is None or _WORKER_CONTROL_TABLE is None:
        raise RuntimeError("Calibration worker state is not initialized.")
    if _WORKER_BASE_CONFIG is None:
        raise RuntimeError("Base evidence configuration is unavailable.")
    indices = _WORKER_CONTROL_INDICES.get(chrom)
    if indices is None or len(indices) == 0:
        return order, chrom, {}

    controls = _WORKER_CONTROL_TABLE.iloc[indices]
    chromosome_densities = load_chromosome_map(_WORKER_CACHE, chrom)
    output: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for _family_id, records in controls.groupby(
        "family_id",
        sort=False,
        observed=True,
    ):
        for _, row in records.iterrows():
            sample_profiles = _sample_profiles(
                row=row,
                samples=_WORKER_SAMPLES,
                chromosome_densities=chromosome_densities,
                genepred_records=_WORKER_GENEPRED,
                post_stop_codons=_WORKER_POST_STOP_CODONS,
            )
            for unit in _WORKER_UNITS:
                nt_profile, downstream, context = _unit_profile(
                    sample_profiles,
                    unit,
                )
                score = _score_profile(
                    row=row,
                    nt_profile=nt_profile,
                    downstream=downstream,
                    release_context=context,
                    config=_WORKER_BASE_CONFIG,
                    post_stop_codons=_WORKER_POST_STOP_CODONS,
                )
                output[(unit.scope, unit.name)].append(
                    score.adaptive.metrics
                )
    return order, chrom, output


def _evaluate_chromosome_worker(
    payload: tuple[int, str],
) -> ChromosomeShardResult:
    """Evaluate all families on one chromosome and write local shards."""
    order, chrom = payload
    if _WORKER_CACHE is None or _WORKER_REPRESENTATIVES is None:
        raise RuntimeError("Evidence worker state is not initialized.")
    indices = _WORKER_CHROM_INDICES.get(chrom)
    if indices is None:
        raise RuntimeError(f"No representative index for chromosome {chrom}.")

    chromosome_densities = load_chromosome_map(_WORKER_CACHE, chrom)
    chrom_table = _WORKER_REPRESENTATIVES.iloc[indices]
    shard_root = Path(_WORKER_SHARD_DIRECTORY)
    shard_paths = {
        "unit": shard_root / f"{order:05d}.unit.tsv",
        "unit_member": shard_root / f"{order:05d}.unit_member.tsv",
        "family": shard_root / f"{order:05d}.family.tsv",
        "member": shard_root / f"{order:05d}.member.tsv",
    }
    evidence_counts: Counter[str] = Counter()
    completed = 0

    with ShardTableWriter(
        shard_paths["unit"], UNIT_COLUMNS
    ) as unit_writer, ShardTableWriter(
        shard_paths["unit_member"], UNIT_MEMBER_COLUMNS
    ) as unit_member_writer, ShardTableWriter(
        shard_paths["family"], FAMILY_COLUMNS
    ) as family_writer, ShardTableWriter(
        shard_paths["member"], MEMBER_COLUMNS
    ) as member_writer:
        for family_id, family_records in chrom_table.groupby(
            "family_id",
            sort=False,
            observed=True,
        ):
            family_id = str(family_id)
            representative_sample_profiles: dict[
                str,
                dict[str, tuple[np.ndarray, np.ndarray, str]],
            ] = {}
            for _, row in family_records.iterrows():
                orf_id = str(row["orf_id"])
                representative_sample_profiles[orf_id] = _sample_profiles(
                    row=row,
                    samples=_WORKER_SAMPLES,
                    chromosome_densities=chromosome_densities,
                    genepred_records=_WORKER_GENEPRED,
                    post_stop_codons=_WORKER_POST_STOP_CODONS,
                )

            unit_rows: list[dict[str, object]] = []
            unit_scores: dict[
                tuple[str, str],
                dict[str, ProfileScore],
            ] = {}
            unit_extensions: dict[
                tuple[str, str],
                dict[str, dict[str, object]],
            ] = {}
            write_unit_members = (
                _WORKER_UNIT_MEMBER_MODE == "all"
                or (
                    _WORKER_UNIT_MEMBER_MODE == "family"
                    and len(family_records) > 1
                )
            )

            for unit in _WORKER_UNITS:
                config = _WORKER_CONFIGS[(unit.scope, unit.name)]
                scores: dict[str, ProfileScore] = {}
                for _, row in family_records.iterrows():
                    orf_id = str(row["orf_id"])
                    nt_profile, downstream, context = _unit_profile(
                        representative_sample_profiles[orf_id],
                        unit,
                    )
                    scores[orf_id] = _score_profile(
                        row=row,
                        nt_profile=nt_profile,
                        downstream=downstream,
                        release_context=context,
                        config=config,
                        post_stop_codons=_WORKER_POST_STOP_CODONS,
                    )
                unit_row, extensions = _family_unit_result(
                    unit=unit,
                    family_records=family_records,
                    scores=scores,
                    config=config,
                )
                unit_writer.write(unit_row)
                if write_unit_members:
                    for member_evidence in _unit_member_rows(
                        unit=unit,
                        family_records=family_records,
                        unit_result=unit_row,
                        scores=scores,
                        extensions=extensions,
                    ):
                        unit_member_writer.write(member_evidence)
                unit_rows.append(unit_row)
                unit_scores[(unit.scope, unit.name)] = scores
                unit_extensions[(unit.scope, unit.name)] = extensions

            consensus = _consensus_family_result(
                family_records=family_records,
                unit_rows=unit_rows,
                mode=_WORKER_MODE,
            )
            family_writer.write(consensus)
            evidence_counts[str(consensus["family_translation_evidence"])] += 1
            selected_key = (
                str(consensus["consensus_scope"]),
                str(consensus["consensus_unit"]),
            )
            selected_unit = next(
                row
                for row in unit_rows
                if (row["evidence_scope"], row["evidence_unit"])
                == selected_key
            )
            for member_row in _member_consensus_rows(
                family_records=family_records,
                family_members=_available_family_members(family_id),
                selected_unit=selected_unit,
                selected_scores=unit_scores[selected_key],
                extensions=unit_extensions[selected_key],
            ):
                member_writer.write(member_row)
            completed += 1

    return ChromosomeShardResult(
        order=order,
        chrom=chrom,
        family_count=completed,
        evidence_counts=dict(evidence_counts),
        shards={key: str(value) for key, value in shard_paths.items()},
    )


def _effective_workers(
    requested: int,
    cache: DensityCacheSet,
    chromosomes: Sequence[str],
) -> int:
    """Choose a memory-aware chromosome worker count."""
    if not chromosomes:
        return 1
    workers = min(max(1, int(requested)), len(chromosomes))
    available = available_memory_bytes()
    if available is None:
        return workers
    largest_cache = max(
        chromosome_cache_bytes(cache, chrom)
        for chrom in chromosomes
    )
    estimated_per_worker = max(
        256 * 1024 * 1024,
        int(largest_cache * 1.5),
    )
    memory_workers = max(1, int(available * 0.65) // estimated_per_worker)
    return max(1, min(workers, memory_workers))


def _merge_shards(
    final_path: str | Path,
    columns: Sequence[str],
    shard_paths: Sequence[str],
) -> None:
    """Merge header-free chromosome shards into one atomic output table."""
    final = Path(final_path)
    temporary = final.with_name(final.name + ".tmp")
    final.parent.mkdir(parents=True, exist_ok=True)
    compressed = final.name.lower().endswith(".gz")
    if compressed:
        output = gzip.open(temporary, "wb", compresslevel=1)
    else:
        output = temporary.open("wb", buffering=_BUFFER_BYTES)
    try:
        output.write(("\t".join(columns) + "\n").encode("utf-8"))
        for shard_path in shard_paths:
            with Path(shard_path).open("rb") as source:
                shutil.copyfileobj(source, output, length=_BUFFER_BYTES)
    finally:
        output.close()
    temporary.replace(final)


def _task_progress_message(
    completed: int,
    total: int,
    progress_label: str,
) -> str:
    """Build one ArgsParser-compatible task progress message.

    Parameters
    ----------
    completed : int
        Number of completed tasks.
    total : int
        Total number of tasks.
    progress_label : str
        Human-readable task unit.

    Returns
    -------
    str
        Progress text accepted by ``progress_print(message)``.
    """
    completed_value = max(0, int(completed))
    total_value = max(0, int(total))
    if total_value == 0:
        return f"{progress_label}: 0/0."
    completed_value = min(completed_value, total_value)
    percentage = completed_value * 100.0 / total_value
    return (
        f"{progress_label}: "
        f"{completed_value:,}/{total_value:,} "
        f"({percentage:.1f}%)."
    )


def _report_task_progress(
    completed: int,
    total: int,
    progress_label: str,
) -> None:
    """Print one task-progress update through the shared output API."""
    progress_print(
        _task_progress_message(
            completed=completed,
            total=total,
            progress_label=progress_label,
        )
    )


def _run_tasks(
    worker,
    payloads: Sequence[tuple[int, str]],
    workers: int,
    progress_label: str,
):
    """Run chromosome tasks sequentially or with forked processes."""
    payload_list = list(payloads)
    total_tasks = len(payload_list)
    if total_tasks == 0:
        _report_task_progress(
            completed=0,
            total=0,
            progress_label=progress_label,
        )
        return []

    context = _fork_context()
    if workers <= 1 or context is None:
        results = []
        for completed, payload in enumerate(payload_list, start=1):
            results.append(worker(payload))
            _report_task_progress(
                completed=completed,
                total=total_tasks,
                progress_label=progress_label,
            )
        return results

    results = []
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
    ) as executor:
        futures = {
            executor.submit(worker, payload): payload
            for payload in payload_list
        }
        completed = 0
        try:
            for future in as_completed(futures):
                results.append(future.result())
                completed += 1
                _report_task_progress(
                    completed=completed,
                    total=total_tasks,
                    progress_label=progress_label,
                )
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    return results


def run_family_riboseq_evidence(args: object) -> EvidencePipelineResult:
    """Run cached, parallel family-aware adaptive evidence analysis."""
    message_print("Load smorf_cluster family and member tables.")
    family_input = read_family_input(
        family_table=getattr(args, "family_table"),
        family_members=getattr(args, "family_members"),
        orf_source=getattr(args, "orf_source", None),
    )
    message_print(
        "Loaded family inputs: "
        f"families={family_input.family_primary.shape[0]:,}, "
        f"representatives={family_input.representatives.shape[0]:,}, "
        f"members={family_input.members.shape[0]:,}."
    )

    message_print("Read density sample and group definitions.")
    samples = read_density_design(args)
    mode = resolve_analysis_mode(
        getattr(args, "analysis_mode", "auto"),
        len(samples),
    )
    units = _analysis_units(samples, mode)
    base_config = _build_base_config(args)
    positive_category = str(
        getattr(args, "positive_category", "annotated_ORF")
    )
    positive_max_controls = int(
        getattr(args, "positive_max_controls", 5000)
    )
    post_stop_codons = int(getattr(args, "post_stop_codons", 10))
    requested_threads = max(1, int(getattr(args, "threads", 1)))
    unit_member_mode = str(
        getattr(args, "unit_member_mode", "family")
    ).lower()

    message_print("Read required transcript models from genePred.")
    genepred_path = getattr(args, "genepred", None)
    coord_mode = str(getattr(args, "coord_mode", "0based-half-open"))
    required_transcripts = set(
        family_input.representatives["transcript_id"].astype(str)
    )
    genepred_records = (
        read_genepred(
            genepred_path,
            coord_mode,
            required_names=required_transcripts,
        )
        if genepred_path
        else {}
    )
    message_print(
        f"Loaded transcript models: {len(genepred_records):,}."
    )

    output_prefix = Path(str(getattr(args, "output")))
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_paths = {
        "unit": Path(f"{output_prefix}.unit_evidence.txt.gz"),
        "unit_member": Path(f"{output_prefix}.unit_member_evidence.txt.gz"),
        "family": Path(f"{output_prefix}.family_evidence.txt"),
        "member": Path(f"{output_prefix}.member_evidence.txt.gz"),
        "calibration": Path(f"{output_prefix}.calibration.txt"),
        "summary": Path(f"{output_prefix}.evidence_summary.txt"),
    }

    required_chromosomes = tuple(
        sorted(
            family_input.representatives["chrom"].astype(str).unique()
        )
    )
    with tempfile.TemporaryDirectory(
        prefix="smorf_evidence.",
        dir=str(output_prefix.parent.resolve()),
    ) as work_directory:
        message_print(
            "Build reusable density caches: "
            f"tracks={sum(len(sample.tracks) for sample in samples):,}, "
            f"workers={requested_threads:,}."
        )

        def cache_progress(completed: int, total: int, track) -> None:
            message_print(
                f"Cached density {completed:,}/{total:,}: "
                f"{track.sample}/{track.strand}, "
                f"chromosomes={len(track.chromosomes):,}."
            )

        cache = build_density_cache(
            samples=samples,
            work_directory=work_directory,
            required_chromosomes=required_chromosomes,
            threads=requested_threads,
            progress_callback=cache_progress,
        )
        message_print(
            f"Density cache ready: {cache.byte_size / (1024 ** 3):.2f} GiB."
        )

        global _WORKER_CACHE
        global _WORKER_REPRESENTATIVES
        global _WORKER_CHROM_INDICES
        global _WORKER_CONTROL_TABLE
        global _WORKER_CONTROL_INDICES
        global _WORKER_MEMBERS_INDEXED
        global _WORKER_SAMPLES
        global _WORKER_UNITS
        global _WORKER_CONFIGS
        global _WORKER_BASE_CONFIG
        global _WORKER_GENEPRED
        global _WORKER_POST_STOP_CODONS
        global _WORKER_MODE
        global _WORKER_SHARD_DIRECTORY
        global _WORKER_UNIT_MEMBER_MODE

        _WORKER_CACHE = cache
        _WORKER_REPRESENTATIVES = family_input.representatives
        _WORKER_CHROM_INDICES = _build_chrom_indices(
            family_input.representatives
        )
        _WORKER_MEMBERS_INDEXED = family_input.members_indexed
        _WORKER_SAMPLES = tuple(samples)
        _WORKER_UNITS = tuple(units)
        _WORKER_BASE_CONFIG = base_config
        _WORKER_GENEPRED = genepred_records
        _WORKER_POST_STOP_CODONS = post_stop_codons
        _WORKER_MODE = mode
        _WORKER_UNIT_MEMBER_MODE = unit_member_mode

        controls_table = _select_positive_controls(
            representatives=family_input.representatives,
            category=positive_category,
            maximum=positive_max_controls,
        )
        _WORKER_CONTROL_TABLE = controls_table
        _WORKER_CONTROL_INDICES = _build_chrom_indices(controls_table)
        control_chromosomes = tuple(
            chrom
            for chrom in required_chromosomes
            if chrom in _WORKER_CONTROL_INDICES
        )
        calibration_workers = _effective_workers(
            requested_threads,
            cache,
            control_chromosomes,
        )
        message_print(
            "Calibrate thresholds with annotated mORFs: "
            f"controls={len(controls_table):,}, "
            f"chromosomes={len(control_chromosomes):,}, "
            f"workers={calibration_workers:,}."
        )
        calibration_results = _run_tasks(
            _calibration_chromosome_worker,
            list(enumerate(control_chromosomes)),
            calibration_workers,
            "calibration chromosomes",
        )
        controls: dict[
            tuple[str, str],
            list[dict[str, object]],
        ] = defaultdict(list)
        for _order, _chrom, result in calibration_results:
            for key, metrics in result.items():
                controls[key].extend(metrics)

        configs: dict[tuple[str, str], AdaptiveEvidenceConfig] = {}
        calibration_rows: list[dict[str, object]] = []
        sample_group = {sample.sample: sample.group for sample in samples}
        for unit in units:
            calibrated, calibration = base_config.with_calibration(
                controls.get((unit.scope, unit.name), [])
            )
            configs[(unit.scope, unit.name)] = calibrated
            calibration_rows.append(
                {
                    "evidence_scope": unit.scope,
                    "evidence_unit": unit.name,
                    "group": (
                        sample_group.get(unit.name, unit.group)
                        if unit.scope == "sample"
                        else unit.group
                    ),
                    **calibration,
                }
            )
        _WORKER_CONFIGS = configs

        with AtomicTableWriter(
            output_paths["calibration"],
            CALIBRATION_COLUMNS,
        ) as calibration_writer:
            for row in calibration_rows:
                calibration_writer.write(row)

        chromosome_workers = _effective_workers(
            requested_threads,
            cache,
            required_chromosomes,
        )
        shard_directory = Path(work_directory) / "evidence_shards"
        shard_directory.mkdir(parents=True, exist_ok=True)
        _WORKER_SHARD_DIRECTORY = str(shard_directory)
        message_print(
            "Evaluate adaptive family evidence: "
            f"chromosomes={len(required_chromosomes):,}, "
            f"workers={chromosome_workers:,}."
        )
        chromosome_results = _run_tasks(
            _evaluate_chromosome_worker,
            list(enumerate(required_chromosomes)),
            chromosome_workers,
            "evidence chromosomes",
        )
        chromosome_results.sort(key=lambda result: result.order)

        _merge_shards(
            output_paths["unit"],
            UNIT_COLUMNS,
            [result.shards["unit"] for result in chromosome_results],
        )
        _merge_shards(
            output_paths["unit_member"],
            UNIT_MEMBER_COLUMNS,
            [result.shards["unit_member"] for result in chromosome_results],
        )
        _merge_shards(
            output_paths["family"],
            FAMILY_COLUMNS,
            [result.shards["family"] for result in chromosome_results],
        )
        _merge_shards(
            output_paths["member"],
            MEMBER_COLUMNS,
            [result.shards["member"] for result in chromosome_results],
        )

    evidence_counts: Counter[str] = Counter()
    for result in chromosome_results:
        evidence_counts.update(result.evidence_counts)
    _write_summary(
        output_paths["summary"],
        {
            "analysis_mode": mode,
            "sample_count": len(samples),
            "group_count": len({sample.group for sample in samples}),
            "analysis_unit_count": len(units),
            "family_count": family_input.family_primary.shape[0],
            "representative_count": family_input.representatives.shape[0],
            "member_count": family_input.members.shape[0],
            "positive_category": positive_category,
            "positive_control_sampled": len(controls_table),
            "density_cache_bytes": cache.byte_size,
            "requested_threads": requested_threads,
            "effective_workers": chromosome_workers,
            "unit_member_mode": unit_member_mode,
            "window_codons": base_config.window_codons,
            "window_step_codons": base_config.window_step_codons,
            "short_max_codons": base_config.short_max_codons,
            "long_min_codons": base_config.long_min_codons,
        },
        evidence_counts,
    )
    return EvidencePipelineResult(
        mode=mode,
        sample_count=len(samples),
        family_count=family_input.family_primary.shape[0],
        representative_count=family_input.representatives.shape[0],
        outputs=tuple(str(path) for path in output_paths.values()),
    )
