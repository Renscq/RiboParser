#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Evaluate sample-level and cross-sample smORF translation evidence.
# Input: Filtered ORF table, P-site density tracks, and optional genePred annotation.
# Output: Sample-level evidence table and ORF-level evidence summary.

"""Main pipeline for smORF Ribo-seq translation-evidence analysis."""

from __future__ import annotations

import csv
import gzip
import math
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from types import TracebackType
from typing import Iterator, TextIO

import numpy as np
import pandas as pd

from utils.ribo.ArgsParser import message_print, progress_print, warning_print
from .smorf_riboseq_constants import (
    DensityTrack,
    EVIDENCE_LEVELS,
    EvidenceThresholds,
)
from .smorf_riboseq_density import (
    ChromDensity,
    iter_density_chromosomes,
    prepare_density,
)
from .smorf_riboseq_io import (
    GenePredRecord,
    get_orf_blocks,
    peek_genepred_names,
    read_density_list,
    read_genepred,
    read_orf_table,
    sample_output_path,
    sanitize_sample_name,
)
from .smorf_riboseq_metrics import (
    calculate_pausing,
    calculate_periodicity,
    calculate_release,
    classify_coverage_shape,
    classify_translation_evidence,
    quantify_profile,
)
from .smorf_riboseq_profile import (
    extract_transcript_downstream_profile,
    extract_transcript_profile,
    get_coding_nt_length,
    make_codon_profile,
)


SAMPLE_OUTPUT_COLUMNS = [
    "sample",
    "density_strand",
    "orf_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "category",
    "genomic_start",
    "genomic_end",
    "nt_length",
    "coding_nt_length",
    "coding_codon_count",
    "block_source",
    "profile_status",
    "profile_reason",
    "rpf_sum",
    "rpf_mean",
    "rpf_per_codon",
    "covered_nt",
    "coverage_ratio",
    "covered_codon",
    "covered_codon_ratio",
    "max_density",
    "frame0_density",
    "frame1_density",
    "frame2_density",
    "frame0_ratio",
    "frame1_ratio",
    "frame2_ratio",
    "frame0_vs_alt_ratio",
    "periodicity_score",
    "periodicity_label",
    "periodicity_evaluable",
    "start_codon_mean",
    "body_codon_mean",
    "pre_stop_codon_mean",
    "start_pause_ratio",
    "stop_pause_ratio",
    "start_pause_label",
    "stop_pause_label",
    "pausing_label",
    "pausing_evaluable",
    "post_stop_mean",
    "release_ratio",
    "release_drop_score",
    "release_label",
    "release_evaluable",
    "release_context",
    "post_stop_nt",
    "codon_gini",
    "max_to_mean_ratio",
    "top10_fraction",
    "coverage_shape",
    "shape_evaluable",
    "evidence_components",
    "evidence_reason",
    "translation_evidence",
]


@dataclass(frozen=True, slots=True)
class SampleEvidenceResult:
    """Store one sample's output information.

    Attributes:
        sample: Density-list sample identifier.
        output: Sample-level evidence path.
        written_rows: Number of evidence rows written.
        sorted_track_count: Number of tracks automatically sorted.
    """

    sample: str
    output: str
    written_rows: int
    sorted_track_count: int


@dataclass(frozen=True, slots=True)
class EvidenceRunResult:
    """Store all sample-specific output information.

    Attributes:
        samples: Per-sample results in density-list order.
        threads: Number of worker processes used.
    """

    samples: tuple[SampleEvidenceResult, ...]
    threads: int


class _EvidenceWriter:
    """Atomically write one sample-level evidence table."""

    def __init__(
        self,
        sample_output: str | Path,
    ) -> None:
        self.sample_final = Path(sample_output)
        self.sample_final.parent.mkdir(parents=True, exist_ok=True)
        self.sample_temp = self.sample_final.with_name(
            self.sample_final.name + ".tmp"
        )
        self.sample_handle: TextIO | None = None
        self.sample_writer: csv.DictWriter | None = None

    @staticmethod
    def _open_text(
        path: Path,
        compressed: bool,
    ) -> TextIO:
        """Open one output file in text mode."""
        if compressed:
            return gzip.open(path, "wt", encoding="utf-8", newline="")
        return path.open("w", encoding="utf-8", newline="")

    def __enter__(self) -> "_EvidenceWriter":
        """Open the temporary sample output."""
        self.sample_handle = self._open_text(
            self.sample_temp,
            str(self.sample_final).lower().endswith(".gz"),
        )
        self.sample_writer = csv.DictWriter(
            self.sample_handle,
            fieldnames=SAMPLE_OUTPUT_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        self.sample_writer.writeheader()
        return self

    def write_sample(self, result: dict[str, object]) -> None:
        """Write one sample-level evidence record.

        Args:
            result: Evidence result dictionary.
        """
        if self.sample_writer is None:
            raise RuntimeError("Evidence writer is not active.")
        self.sample_writer.writerow(
            {
                column: result.get(column, "")
                for column in SAMPLE_OUTPUT_COLUMNS
            }
        )

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        """Commit a successful output or delete the temporary file."""
        if self.sample_handle is not None:
            self.sample_handle.close()
            self.sample_handle = None
            self.sample_writer = None

        if exception_type is None:
            self.sample_temp.replace(self.sample_final)
        else:
            self.sample_temp.unlink(missing_ok=True)
        return False


def build_thresholds(args: object) -> EvidenceThresholds:
    """Build thresholds from a preset and optional overrides.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Validated evidence thresholds.
    """
    thresholds = EvidenceThresholds.from_mode(
        getattr(args, "evidence_mode", "balanced")
    ).with_overrides(
        min_rpf_sum=getattr(args, "min_rpf_sum", None),
        min_covered_codon=getattr(
            args,
            "min_covered_codon",
            None,
        ),
        min_codon_coverage=getattr(
            args,
            "min_codon_coverage",
            None,
        ),
        moderate_periodicity=getattr(
            args,
            "moderate_periodicity",
            None,
        ),
        strong_periodicity=getattr(
            args,
            "strong_periodicity",
            None,
        ),
    )
    return thresholds


def _empty_metrics() -> dict[str, object]:
    """Return zero/NA metrics for an invalid profile."""
    return {
        "rpf_sum": 0.0,
        "rpf_mean": 0.0,
        "rpf_per_codon": 0.0,
        "covered_nt": 0,
        "coverage_ratio": 0.0,
        "covered_codon": 0,
        "covered_codon_ratio": 0.0,
        "max_density": 0.0,
        "frame0_density": 0.0,
        "frame1_density": 0.0,
        "frame2_density": 0.0,
        "frame0_ratio": 0.0,
        "frame1_ratio": 0.0,
        "frame2_ratio": 0.0,
        "frame0_vs_alt_ratio": 0.0,
        "periodicity_score": 0.0,
        "periodicity_label": "NA",
        "periodicity_evaluable": False,
        "start_codon_mean": 0.0,
        "body_codon_mean": 0.0,
        "pre_stop_codon_mean": 0.0,
        "start_pause_ratio": 0.0,
        "stop_pause_ratio": 0.0,
        "start_pause_label": "NA",
        "stop_pause_label": "NA",
        "pausing_label": "NA",
        "pausing_evaluable": False,
        "post_stop_mean": 0.0,
        "release_ratio": 0.0,
        "release_drop_score": 0.0,
        "release_label": "NA",
        "release_evaluable": False,
        "release_context": "unavailable",
        "post_stop_nt": 0,
        "codon_gini": 0.0,
        "max_to_mean_ratio": 0.0,
        "top10_fraction": 0.0,
        "coverage_shape": "NA",
        "shape_evaluable": False,
        "evidence_components": "none",
        "evidence_reason": "invalid_profile",
        "translation_evidence": "NoEvidence",
    }


def _base_result(
    row: dict[str, object],
    track: DensityTrack,
    block_source: str,
) -> dict[str, object]:
    """Build common sample-level output fields."""
    return {
        "sample": track.sample,
        "density_strand": track.strand,
        "orf_id": row["orf_id"],
        "gene_id": row["gene_id"],
        "transcript_id": row["transcript_id"],
        "chrom": row["chrom"],
        "strand": row["strand"],
        "category": row.get("category", "."),
        "genomic_start": int(row["genomic_start"]),
        "genomic_end": int(row["genomic_end"]),
        "nt_length": int(row["nt_length"]),
        "block_source": block_source,
    }


def score_one_orf(
    row: dict[str, object],
    density: ChromDensity | None,
    track: DensityTrack,
    genepred_records: dict[str, GenePredRecord],
    thresholds: EvidenceThresholds,
    post_stop_codons: int,
) -> dict[str, object]:
    """Score one ORF in one sample.

    Args:
        row: ORF metadata record.
        density: Chromosome density or ``None`` when absent.
        track: Density-track metadata.
        genepred_records: Optional genePred mapping.
        thresholds: Evidence thresholds.
        post_stop_codons: Downstream transcript window.

    Returns:
        Sample-level evidence result.
    """
    try:
        starts, ends, block_source = get_orf_blocks(
            row,
            genepred_records,
        )
        base = _base_result(row, track, block_source)

        nt_profile = extract_transcript_profile(
            density=density,
            starts=starts,
            ends=ends,
            strand=str(row["strand"]),
            orf_id=str(row["orf_id"]),
        )
        declared_length = int(row["nt_length"])
        if len(nt_profile) != declared_length:
            raise ValueError(
                f"extracted length {len(nt_profile)} does not match "
                f"declared nt_length {declared_length}"
            )

        coding_nt_length = get_coding_nt_length(
            pd.Series(row),
            len(nt_profile),
        )
        codon_profile = make_codon_profile(
            nt_profile,
            coding_nt_length,
        )

        transcript = genepred_records.get(
            str(row["transcript_id"])
        )
        complete = str(
            row.get("completeness", "")
        ).strip().lower() in {
            "",
            "complete",
            "cmpl",
            "full",
            "true",
            "yes",
        }
        if density is None:
            downstream_profile = np.zeros(0, dtype=np.float32)
            release_context = "no_density"
        elif complete:
            downstream_profile, release_context = (
                extract_transcript_downstream_profile(
                    density=density,
                    orf_starts=starts,
                    orf_ends=ends,
                    strand=str(row["strand"]),
                    transcript=transcript,
                    nt_window=int(post_stop_codons) * 3,
                    orf_id=str(row["orf_id"]),
                )
            )
        else:
            downstream_profile = np.zeros(0, dtype=np.float32)
            release_context = "partial_orf"

        quant = quantify_profile(
            nt_profile,
            codon_profile,
            coding_nt_length,
        )
        periodicity = calculate_periodicity(
            nt_profile,
            coding_nt_length,
            thresholds,
        )
        pausing = calculate_pausing(
            codon_profile,
            thresholds,
        )
        release = calculate_release(
            codon_profile=codon_profile,
            downstream_nt_profile=downstream_profile,
            post_stop_codons=post_stop_codons,
            thresholds=thresholds,
            context=release_context,
        )
        shape = classify_coverage_shape(
            codon_profile,
            thresholds,
        )
        evidence, components, reason = (
            classify_translation_evidence(
                quant=quant,
                periodicity=periodicity,
                release=release,
                coverage_shape=shape,
                thresholds=thresholds,
            )
        )

        base.update(
            {
                "coding_nt_length": int(coding_nt_length),
                "coding_codon_count": int(len(codon_profile)),
                "profile_status": (
                    "NoDensityChrom"
                    if density is None
                    else "OK"
                ),
                "profile_reason": (
                    "chromosome_absent_from_density"
                    if density is None
                    else "PASS"
                ),
            }
        )
        base.update(quant)
        base.update(periodicity)
        base.update(pausing)
        base.update(release)
        base.update(shape)
        base["evidence_components"] = components
        base["evidence_reason"] = reason
        base["translation_evidence"] = evidence
        return base

    except Exception as error:
        base = _base_result(row, track, "invalid")
        base.update(
            {
                "coding_nt_length": 0,
                "coding_codon_count": 0,
                "profile_status": "InvalidProfile",
                "profile_reason": str(error),
            }
        )
        base.update(_empty_metrics())
        return base


def _track_orf_table(
    orf_table: pd.DataFrame,
    track: DensityTrack,
) -> pd.DataFrame:
    """Return ORFs applicable to one density track."""
    if track.strand in {"+", "-"}:
        return orf_table.loc[
            orf_table["strand"].eq(track.strand)
        ]
    return orf_table


def process_one_track(
    orf_table: pd.DataFrame,
    track: DensityTrack,
    genepred_records: dict[str, GenePredRecord],
    thresholds: EvidenceThresholds,
    post_stop_codons: int,
) -> Iterator[tuple[int, dict[str, object]]]:
    """Yield evidence results for one density track.

    The density file is read once. ORFs on chromosomes absent from the track
    receive explicit ``NoDensityChrom`` results.

    Args:
        orf_table: Validated ORF table containing internal ``_orf_index``.
        track: Density track.
        genepred_records: Optional genePred mapping.
        thresholds: Evidence thresholds.
        post_stop_codons: Downstream window.

    Yields:
        ORF index and sample-level result.
    """
    sub_table = _track_orf_table(orf_table, track)
    if sub_table.empty:
        return

    groups = {
        str(chrom): np.asarray(indices, dtype=np.int64)
        for chrom, indices in sub_table.groupby(
            "chrom",
            sort=False,
        ).indices.items()
    }
    seen_chromosomes: set[str] = set()

    message_print(
        "Processing sample={sample}, strand={strand}, "
        "ORFs={orfs:,}, file={path}".format(
            sample=track.sample,
            strand=track.strand,
            orfs=len(sub_table),
            path=track.path,
        )
    )

    for chrom_density in iter_density_chromosomes(
        track.path,
        track.file_format,
    ):
        chrom = chrom_density.chrom
        seen_chromosomes.add(chrom)
        row_indices = groups.get(chrom)
        if row_indices is None:
            continue

        for row in sub_table.iloc[row_indices].to_dict(
            orient="records"
        ):
            result = score_one_orf(
                row=row,
                density=chrom_density,
                track=track,
                genepred_records=genepred_records,
                thresholds=thresholds,
                post_stop_codons=post_stop_codons,
            )
            yield int(row["_orf_index"]), result

    missing_chromosomes = set(groups).difference(seen_chromosomes)
    for chrom in sorted(missing_chromosomes):
        for row in sub_table.iloc[groups[chrom]].to_dict(
            orient="records"
        ):
            result = score_one_orf(
                row=row,
                density=None,
                track=track,
                genepred_records=genepred_records,
                thresholds=thresholds,
                post_stop_codons=post_stop_codons,
            )
            yield int(row["_orf_index"]), result


def _minimum_evidence_rank(value: str) -> int:
    """Convert a minimum evidence label to a numeric rank."""
    mapping = {
        "all": 0,
        "low": 1,
        "medium": 2,
        "high": 3,
    }
    return mapping[str(value).lower()]


_WORKER_ORF_TABLE: pd.DataFrame | None = None
_WORKER_GENEPRED: dict[str, GenePredRecord] | None = None
_WORKER_THRESHOLDS: EvidenceThresholds | None = None
_WORKER_OUTPUT_TEMPLATE: str | None = None
_WORKER_POST_STOP_CODONS = 10
_WORKER_MINIMUM_RANK = 0


def _common_identifier_prefix(values: pd.Series) -> str:
    """Return a stable non-numeric prefix from identifier examples.

    Args:
        values: Identifier values.

    Returns:
        Common prefix with trailing digits and separators removed.
    """
    examples = [
        str(value)
        for value in values.head(200).tolist()
        if str(value)
    ]
    if not examples:
        return ""
    prefix = os.path.commonprefix(examples)
    return prefix.rstrip("0123456789._-")


def _looks_like_scanner_orf_genepred(
    path: str | Path,
    orf_table: pd.DataFrame,
) -> bool:
    """Detect scanner ORF genePred output from sampled record names.

    Args:
        path: genePred path.
        orf_table: Filtered ORF table.

    Returns:
        ``True`` when names strongly resemble ORF IDs and not transcript IDs.
    """
    sampled_names = peek_genepred_names(path, maximum_records=200)
    if not sampled_names:
        return False

    transcript_examples = set(
        orf_table["transcript_id"].astype(str).head(10_000)
    )
    if any(name in transcript_examples for name in sampled_names):
        return False

    orf_prefix = _common_identifier_prefix(
        orf_table["orf_id"].astype(str)
    )
    if len(orf_prefix) < 3:
        return False

    matching = sum(
        name.startswith(orf_prefix)
        for name in sampled_names
    )
    return matching / len(sampled_names) >= 0.80


def _sample_groups(
    tracks: list[DensityTrack],
) -> list[tuple[str, tuple[DensityTrack, ...]]]:
    """Group density tracks by sample while preserving input order.

    Args:
        tracks: Validated density tracks.

    Returns:
        Sample and track tuples.
    """
    grouped: dict[str, list[DensityTrack]] = {}
    for track in tracks:
        grouped.setdefault(track.sample, []).append(track)
    return [
        (sample, tuple(sample_tracks))
        for sample, sample_tracks in grouped.items()
    ]


def _validate_sample_output_names(
    output_template: str | Path,
    sample_groups: list[tuple[str, tuple[DensityTrack, ...]]],
) -> None:
    """Ensure sanitized sample names do not produce output collisions.

    Args:
        output_template: User output template.
        sample_groups: Grouped sample tracks.

    Raises:
        ValueError: If two samples map to one output path.
    """
    observed: dict[str, str] = {}
    for sample, _ in sample_groups:
        output_path = sample_output_path(output_template, sample)
        if output_path in observed:
            raise ValueError(
                "Sample names produce the same output path after "
                f"sanitization: {observed[output_path]} and {sample}"
            )
        observed[output_path] = sample


def _configure_worker(
    orf_table: pd.DataFrame,
    genepred_records: dict[str, GenePredRecord],
    thresholds: EvidenceThresholds,
    output_template: str,
    post_stop_codons: int,
    minimum_rank: int,
) -> None:
    """Configure process-global read-only analysis objects.

    Args:
        orf_table: Validated ORF table.
        genepred_records: Selected transcript annotation.
        thresholds: Evidence thresholds.
        output_template: User output template.
        post_stop_codons: Release window.
        minimum_rank: Minimum sample-level output rank.
    """
    global _WORKER_ORF_TABLE
    global _WORKER_GENEPRED
    global _WORKER_THRESHOLDS
    global _WORKER_OUTPUT_TEMPLATE
    global _WORKER_POST_STOP_CODONS
    global _WORKER_MINIMUM_RANK

    _WORKER_ORF_TABLE = orf_table
    _WORKER_GENEPRED = genepred_records
    _WORKER_THRESHOLDS = thresholds
    _WORKER_OUTPUT_TEMPLATE = str(output_template)
    _WORKER_POST_STOP_CODONS = int(post_stop_codons)
    _WORKER_MINIMUM_RANK = int(minimum_rank)


def _process_sample(
    sample_task: tuple[str, tuple[DensityTrack, ...]],
) -> SampleEvidenceResult:
    """Process all density tracks belonging to one sample.

    Args:
        sample_task: Sample identifier and its one or two tracks.

    Returns:
        Sample-specific output information.

    Raises:
        RuntimeError: If worker configuration is unavailable.
    """
    if (
        _WORKER_ORF_TABLE is None
        or _WORKER_GENEPRED is None
        or _WORKER_THRESHOLDS is None
        or _WORKER_OUTPUT_TEMPLATE is None
    ):
        raise RuntimeError("smORF evidence worker is not configured.")

    sample, tracks = sample_task
    sample_output = sample_output_path(
        _WORKER_OUTPUT_TEMPLATE,
        sample,
    )
    output_directory = Path(sample_output).parent
    sample_rows = 0
    sorted_track_count = 0

    message_print(
        f"Start sample={sample}, "
        f"tracks={len(tracks)}."
    )

    with _EvidenceWriter(
        sample_output=sample_output,
    ) as writer:
        for track in tracks:
            prepared = prepare_density(
                path=track.path,
                file_format=track.file_format,
                work_directory=output_directory,
            )
            if prepared.was_sorted:
                sorted_track_count += 1
                warning_print(
                    "Density file required coordinate sorting: "
                    f"sample={sample}, file={track.path}"
                )

            prepared_track = replace(
                track,
                path=prepared.path,
                file_format=(
                    "bedgraph"
                    if prepared.was_sorted
                    else track.file_format
                ),
            )
            try:
                for orf_index, result in process_one_track(
                    orf_table=_WORKER_ORF_TABLE,
                    track=prepared_track,
                    genepred_records=_WORKER_GENEPRED,
                    thresholds=_WORKER_THRESHOLDS,
                    post_stop_codons=_WORKER_POST_STOP_CODONS,
                ):
                    label_rank = EVIDENCE_LEVELS[
                        str(result["translation_evidence"])
                    ]
                    if label_rank >= _WORKER_MINIMUM_RANK:
                        writer.write_sample(result)
                        sample_rows += 1
            finally:
                prepared.cleanup()

    message_print(
        "Completed sample={sample}, rows={rows:,}, "
        "auto_sorted_tracks={sorted_count}.".format(
            sample=sample,
            rows=sample_rows,
            sorted_count=sorted_track_count,
        )
    )

    return SampleEvidenceResult(
        sample=sample,
        output=sample_output,
        written_rows=sample_rows,
        sorted_track_count=sorted_track_count,
    )


def _multiprocessing_context(
    requested_threads: int,
    sample_count: int,
) -> tuple[multiprocessing.context.BaseContext | None, int]:
    """Resolve a safe multiprocessing context.

    Large ORF and annotation tables are shared copy-on-write through ``fork``.
    Platforms without ``fork`` fall back to one process to avoid serializing
    multi-gigabyte objects into every worker.

    Args:
        requested_threads: User-requested process count.
        sample_count: Number of sample tasks.

    Returns:
        Multiprocessing context and effective worker count.
    """
    workers = max(1, min(int(requested_threads), int(sample_count)))
    if workers <= 1:
        return None, 1

    methods = multiprocessing.get_all_start_methods()
    if sys.platform.startswith("linux") and "fork" in methods:
        return multiprocessing.get_context("fork"), workers

    warning_print(
        "Multiprocessing requires fork for "
        "large shared tables on this platform. Use one process."
    )
    return None, 1


def run_riboseq_evidence(
    args: object,
) -> EvidenceRunResult:
    """Run sample-separated evidence analysis.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Per-sample output paths and effective worker count.
    """
    message_print("Reading ORF table.")
    orf_table = read_orf_table(
        args.orf_table,
        args.coord_mode,
    )
    orf_table["_orf_index"] = np.arange(
        len(orf_table),
        dtype=np.int64,
    )
    message_print(
        f"Loaded PASS ORFs: {len(orf_table):,}"
    )

    tracks = read_density_list(args)
    sample_groups = _sample_groups(tracks)
    _validate_sample_output_names(args.output, sample_groups)

    genepred_records: dict[str, GenePredRecord] = {}
    genepred_path = getattr(args, "genepred", None)
    if genepred_path is not None:
        has_orf_blocks = {
            "exon_starts",
            "exon_ends",
        }.issubset(orf_table.columns)

        if (
            has_orf_blocks
            and _looks_like_scanner_orf_genepred(
                genepred_path,
                orf_table,
            )
        ):
            warning_print(
                "--genepred appears to be the "
                "scanner ORF genePred output. ORF blocks are already in "
                "the message table, so this file is skipped. Use the "
                "original transcript annotation passed to smorf_scanner "
                "-a to enable transcript-aware release analysis."
            )
        else:
            required_names = set(
                orf_table["transcript_id"].astype(str)
            )
            if not has_orf_blocks:
                required_names.update(
                    orf_table["orf_id"].astype(str)
                )

            message_print(
                "Reading selected genePred records. "
                f"Required names: {len(required_names):,}"
            )
            genepred_records = read_genepred(
                genepred_path,
                args.coord_mode,
                required_names=required_names,
            )
            message_print(
                "Retained genePred records: "
                f"{len(genepred_records):,}"
            )
            if not genepred_records:
                warning_print(
                    "No genePred name matched "
                    "transcript_id/required ORF IDs. Release analysis is "
                    "disabled. Use the original transcript genePred passed "
                    "to smorf_scanner -a, or omit --genepred."
                )

    thresholds = build_thresholds(args)
    minimum_rank = _minimum_evidence_rank(
        getattr(args, "min_evidence", "all")
    )
    requested_threads = max(1, int(getattr(args, "threads", 1)))
    context, workers = _multiprocessing_context(
        requested_threads=requested_threads,
        sample_count=len(sample_groups),
    )

    _configure_worker(
        orf_table=orf_table,
        genepred_records=genepred_records,
        thresholds=thresholds,
        output_template=str(args.output),
        post_stop_codons=args.post_stop_codons,
        minimum_rank=minimum_rank,
    )

    if workers == 1:
        results = tuple(
            _process_sample(sample_task)
            for sample_task in sample_groups
        )
    else:
        message_print(
            f"Parallel samples: workers={workers}, "
            f"samples={len(sample_groups)}."
        )
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
        ) as executor:
            results = tuple(
                executor.map(
                    _process_sample,
                    sample_groups,
                    chunksize=1,
                )
            )

    return EvidenceRunResult(
        samples=results,
        threads=workers,
    )
