#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev022
# Function: Resolve candidate-first translation units and audit ORF phase.
# Input: Family tables, scanner ORF table/genePred, and density design table.
# Output: Family evidence, reliable smORFs/genePred, and summary.

"""Bottom-up family-aware smORF evidence engine.

The engine deliberately avoids the previous representative-by-unit workflow.
It builds a SQLite family index, caches each density track once, evaluates one
family scaffold per sample, derives project reliability from independent sample
support across all sample groups, and resolves translated extents without
requiring every adjacent alternative-start segment to pass independently.

Only standard-library modules, NumPy, and the existing RiboParser density reader
are required. A malformed family is recorded as uncertain instead of aborting
an entire genome-wide run.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import shutil
import sqlite3
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Final, Iterable, Iterator, Mapping, Sequence, TextIO

import numpy as np

from utils.ribo.ArgsParser import message_print, progress_print, warning_print
from utils.smorf.smorf_riboseq_density import (
    ChromDensity,
    infer_density_format,
    iter_density_chromosomes,
    prepare_density,
)

SCHEMA_VERSION: Final[int] = 9
STOP_CODONS: Final[frozenset[str]] = frozenset({"TAA", "TAG", "TGA"})
ANNOTATED_CATEGORIES: Final[frozenset[str]] = frozenset(
    {"annotated_orf", "annotated_morf"}
)
FAMILY_BATCH_SIZE: Final[int] = 10_000
SQL_BATCH_SIZE: Final[int] = 50_000
MAX_WORKERS: Final[int] = 8
MEMORY_PER_WORKER_BYTES: Final[int] = 4 * 1024**3
LEVEL_RANK: Final[dict[str, int]] = {
    "NoEvidence": 0,
    "LowConfidence": 1,
    "MediumConfidence": 2,
    "HighConfidence": 3,
}

MANUAL_THRESHOLD_DEFAULTS: Final[dict[str, float | int]] = {
    "min_rpf_sum": 5.0,
    "min_rpf_per_codon": 0.10,
    "min_covered_codon": 3,
    "min_coverage_ratio": 0.10,
    "moderate_periodicity": 0.50,
    "strong_periodicity": 0.60,
    "min_window_rpf": 3.0,
    "min_window_covered": 3,
}
ADVANCED_DEFAULTS: Final[dict[str, float | int]] = {
    "short_max_codons": 30,
    "long_min_codons": 100,
    "window_codons": 20,
    "window_step_codons": 5,
    "min_supported_windows": 2,
    "min_window_gap_codons": 15,
    "min_signal_span": 0.35,
    "localized_span_max": 0.20,
    "localized_top_window_fraction": 0.70,
    "boundary_codons": 5,
    "extent_bins": 5,
    "min_extent_bins": 3,
    "start_resolution_codons": 5,
    "min_frame_margin": 0.10,
    "leading_window_codons": 20,
    "min_exclusion_codons": 10,
    "silent_extension_density_ratio": 0.20,
    "noncanonical_min_exclusive_codons": 8,
    "noncanonical_override_density_ratio": 0.75,
    "noncanonical_override_min_coverage_ratio": 0.50,
    "noncanonical_override_min_frame_margin": 0.25,
    "nested_min_frame_margin": 0.15,
    "nested_min_phase_rpf": 10.0,
    "high_overlap_fraction": 0.70,
    "overlap_min_frame_margin": 0.25,
    "overlap_min_phase_rpf": 15.0,
}
SHORT_MODEL_FLOOR: Final[int] = 20
SHORT_MODEL_CEILING: Final[int] = 30
LONG_MODEL_FLOOR: Final[int] = 90
LONG_MODEL_CEILING: Final[int] = 100
TARGET_LENGTH_CATEGORIES: Final[frozenset[str]] = frozenset(
    {"uorf", "dorf", "lncorf"}
)


@dataclass(frozen=True, slots=True)
class DensityTrack:
    """Describe one P-site density track."""

    sample: str
    group: str
    strand: str
    path: str
    file_format: str


@dataclass(frozen=True, slots=True)
class Thresholds:
    """Store sample-specific evidence thresholds."""

    sample: str
    source: str
    control_count: int
    min_rpf_sum: float
    min_rpf_per_codon: float
    min_covered_codon: int
    min_coverage_ratio: float
    moderate_periodicity: float
    strong_periodicity: float
    min_window_rpf: float
    min_window_covered: int


@dataclass(frozen=True, slots=True)
class EngineConfig:
    """Store biological, calibration, grouping, and runtime configuration."""

    evidence_mode: str = "manual"
    group_column: str = "group"
    reliable_sample: int = 2

    min_rpf_sum: float = 5.0
    min_rpf_per_codon: float = 0.10
    min_covered_codon: int = 3
    min_coverage_ratio: float = 0.10
    moderate_periodicity: float = 0.50
    strong_periodicity: float = 0.60
    min_window_rpf: float = 3.0
    min_window_covered: int = 3

    short_max_codons: int = 30
    long_min_codons: int = 100
    window_codons: int = 20
    window_step_codons: int = 5
    min_supported_windows: int = 2
    min_window_gap_codons: int = 15
    min_signal_span: float = 0.35
    localized_span_max: float = 0.20
    localized_top_window_fraction: float = 0.70
    boundary_codons: int = 5
    extent_bins: int = 5
    min_extent_bins: int = 3
    start_resolution_codons: int = 5
    min_frame_margin: float = 0.10
    leading_window_codons: int = 20
    min_exclusion_codons: int = 10
    silent_extension_density_ratio: float = 0.20
    noncanonical_min_exclusive_codons: int = 8
    noncanonical_override_density_ratio: float = 0.75
    noncanonical_override_min_coverage_ratio: float = 0.50
    noncanonical_override_min_frame_margin: float = 0.25
    nested_min_frame_margin: float = 0.15
    nested_min_phase_rpf: float = 10.0
    high_overlap_fraction: float = 0.70
    overlap_min_frame_margin: float = 0.25
    overlap_min_phase_rpf: float = 15.0
    estimated_short_max_codons: int = 30
    estimated_long_min_codons: int = 100
    hidden_overrides: frozenset[str] = frozenset()

    positive_quantile: float = 0.20
    positive_min_controls: int = 30
    positive_max_controls: int = 5000
    threads: int = 1



@dataclass(frozen=True, slots=True)
class RepresentativeGeometry:
    """Store one family representative in scaffold coordinates."""

    orf_id: str
    transcript_id: str
    category: str
    start_codon: str
    nt_length: int
    coding_nt_length: int
    aa_length: int
    starts: tuple[int, ...]
    ends: tuple[int, ...]
    start_offset: int
    stop_codon: str = "TGA"
    completeness: str = "complete"


@dataclass(frozen=True, slots=True)
class FamilyGeometry:
    """Store one validated family scaffold."""

    family_id: str
    gene_id: str
    chrom: str
    strand: str
    category: str
    family_type: str
    family_size: int
    structural_primary: str
    scaffold_orf_id: str
    scaffold_nt_length: int
    scaffold_coding_nt_length: int
    scaffold_starts: tuple[int, ...]
    scaffold_ends: tuple[int, ...]
    common_body_orf: str
    common_body_start: int
    common_body_end: int
    representatives: tuple[RepresentativeGeometry, ...]


@dataclass(frozen=True, slots=True)
class InvalidFamily:
    """Store a family that cannot be represented by one scaffold."""

    family_id: str
    gene_id: str
    chrom: str
    strand: str
    category: str
    family_type: str
    family_size: int
    structural_primary: str
    reason: str


@dataclass(frozen=True, slots=True)
class SegmentFeatures:
    """Store sufficient statistics for one transcript segment."""

    rpf_sum: float
    rpf_per_codon: float
    covered_codon: int
    coverage_ratio: float
    frame0_density: float
    frame1_density: float
    frame2_density: float
    frame0_ratio: float
    supported_windows: int
    distributed_windows: int
    signal_span: float
    top_window_fraction: float
    start_rpf: float
    body_rpf: float
    end_rpf: float
    localized_only: bool
    score: float


class ScaffoldProfile:
    """Store prefix-indexed codon and frame profiles for one family/sample."""

    __slots__ = (
        "codon_count",
        "codon_prefix",
        "covered_prefix",
        "frame_prefixes",
        "nonzero_codons",
    )

    def __init__(
        self,
        positions: np.ndarray,
        values: np.ndarray,
        coding_nt_length: int,
    ) -> None:
        self.codon_count = max(0, int(coding_nt_length) // 3)
        valid = (
            (positions >= 0)
            & (positions < self.codon_count * 3)
            & np.isfinite(values)
            & (values > 0)
        )
        local_positions = positions[valid].astype(np.int64, copy=False)
        local_values = values[valid].astype(np.float64, copy=False)
        codons = local_positions // 3
        frames = local_positions % 3
        codon_profile = np.bincount(
            codons,
            weights=local_values,
            minlength=self.codon_count,
        ).astype(np.float64, copy=False)
        self.codon_prefix = np.empty(self.codon_count + 1, dtype=np.float64)
        self.codon_prefix[0] = 0.0
        np.cumsum(codon_profile, out=self.codon_prefix[1:])
        self.covered_prefix = np.empty(self.codon_count + 1, dtype=np.int64)
        self.covered_prefix[0] = 0
        np.cumsum(
            (codon_profile > 0).astype(np.int64),
            out=self.covered_prefix[1:],
        )
        self.frame_prefixes: tuple[np.ndarray, ...] = tuple(
            self._prefix(
                np.bincount(
                    codons[frames == frame],
                    weights=local_values[frames == frame],
                    minlength=self.codon_count,
                ).astype(np.float64, copy=False)
            )
            for frame in range(3)
        )
        self.nonzero_codons = np.flatnonzero(codon_profile > 0)

    @staticmethod
    def _prefix(profile: np.ndarray) -> np.ndarray:
        """Return a zero-leading cumulative sum array."""
        output = np.empty(profile.size + 1, dtype=np.float64)
        output[0] = 0.0
        np.cumsum(profile, out=output[1:])
        return output

    def region_features(
        self,
        start: int,
        end: int,
        thresholds: Thresholds,
        config: EngineConfig,
    ) -> SegmentFeatures:
        """Calculate one codon-aligned region from prefix indexes."""
        start_codon = max(0, int(start) // 3)
        end_codon = min(self.codon_count, int(end) // 3)
        codon_count = max(0, end_codon - start_codon)
        if codon_count <= 0:
            return _empty_segment_features()

        rpf_sum = float(
            self.codon_prefix[end_codon] - self.codon_prefix[start_codon]
        )
        if rpf_sum <= 0:
            return _empty_segment_features()
        covered = int(
            self.covered_prefix[end_codon]
            - self.covered_prefix[start_codon]
        )
        coverage_ratio = covered / codon_count
        local_phase = int(start) % 3
        frame_order = (
            local_phase,
            (local_phase + 1) % 3,
            (local_phase + 2) % 3,
        )
        frame_density = [
            float(
                self.frame_prefixes[frame][end_codon]
                - self.frame_prefixes[frame][start_codon]
            )
            for frame in frame_order
        ]
        frame_total = sum(frame_density)
        frame0_ratio = frame_density[0] / frame_total if frame_total else 0.0

        nonzero_left = int(
            np.searchsorted(self.nonzero_codons, start_codon, side="left")
        )
        nonzero_right = int(
            np.searchsorted(self.nonzero_codons, end_codon, side="left")
        )
        if nonzero_right > nonzero_left:
            first_nonzero = int(self.nonzero_codons[nonzero_left])
            last_nonzero = int(self.nonzero_codons[nonzero_right - 1])
            signal_span = (last_nonzero - first_nonzero + 1) / codon_count
        else:
            signal_span = 0.0

        boundary = min(config.boundary_codons, codon_count)
        start_rpf = float(
            self.codon_prefix[start_codon + boundary]
            - self.codon_prefix[start_codon]
        )
        end_rpf = float(
            self.codon_prefix[end_codon]
            - self.codon_prefix[end_codon - boundary]
        )
        body_start = start_codon + boundary
        body_end = max(body_start, end_codon - boundary)
        body_rpf = float(
            self.codon_prefix[body_end] - self.codon_prefix[body_start]
        )

        supported_windows = 0
        supported_starts: list[int] = []
        top_window = 0.0
        if codon_count >= config.window_codons:
            window_starts = list(
                range(
                    0,
                    codon_count - config.window_codons + 1,
                    config.window_step_codons,
                )
            )
            last_start = codon_count - config.window_codons
            if not window_starts or window_starts[-1] != last_start:
                window_starts.append(last_start)
            for relative_start in window_starts:
                left = start_codon + relative_start
                right = left + config.window_codons
                window_sum = float(
                    self.codon_prefix[right] - self.codon_prefix[left]
                )
                top_window = max(top_window, window_sum)
                window_covered = int(
                    self.covered_prefix[right] - self.covered_prefix[left]
                )
                window_frames = [
                    float(
                        self.frame_prefixes[frame][right]
                        - self.frame_prefixes[frame][left]
                    )
                    for frame in frame_order
                ]
                window_total = sum(window_frames)
                window_frame0 = (
                    window_frames[0] / window_total if window_total else 0.0
                )
                if (
                    window_sum >= thresholds.min_window_rpf
                    and window_covered >= thresholds.min_window_covered
                    and window_frame0 >= thresholds.moderate_periodicity
                ):
                    supported_windows += 1
                    supported_starts.append(relative_start)

        top_window_fraction = top_window / rpf_sum if rpf_sum else 0.0
        distributed = 0
        if supported_starts:
            selected = [supported_starts[0]]
            for value in supported_starts[1:]:
                if value - selected[-1] >= config.min_window_gap_codons:
                    selected.append(value)
            distributed = len(selected)
        localized_only = (
            codon_count >= config.long_min_codons
            and signal_span <= config.localized_span_max
            and top_window_fraction >= config.localized_top_window_fraction
        )
        abundance_score = min(
            1.0,
            rpf_sum / max(thresholds.min_rpf_sum, 1e-9),
        )
        coverage_score = min(
            1.0,
            coverage_ratio / max(thresholds.min_coverage_ratio, 1e-9),
        )
        periodicity_score = max(
            0.0,
            min(1.0, (frame0_ratio - 1.0 / 3.0) / (2.0 / 3.0)),
        )
        distribution_score = min(
            1.0,
            max(
                signal_span,
                distributed / max(config.min_supported_windows, 1),
            ),
        )
        score = (
            0.30 * abundance_score
            + 0.25 * coverage_score
            + 0.30 * periodicity_score
            + 0.15 * distribution_score
        )
        return SegmentFeatures(
            rpf_sum=rpf_sum,
            rpf_per_codon=rpf_sum / codon_count,
            covered_codon=covered,
            coverage_ratio=coverage_ratio,
            frame0_density=frame_density[0],
            frame1_density=frame_density[1],
            frame2_density=frame_density[2],
            frame0_ratio=frame0_ratio,
            supported_windows=supported_windows,
            distributed_windows=distributed,
            signal_span=signal_span,
            top_window_fraction=top_window_fraction,
            start_rpf=start_rpf,
            body_rpf=body_rpf,
            end_rpf=end_rpf,
            localized_only=localized_only,
            score=score,
        )


@dataclass(frozen=True, slots=True)
class ExtentResolution:
    """Store the translated-extent and start-site resolution for one family."""

    quant_primary: str
    evidence_primary: str
    translated_extent_status: str
    extent_supported_sample_count: int
    extent_supporting_samples: tuple[str, ...]
    extent_supported_bins: tuple[int, ...]
    extension_rpf_sum: float
    extension_frame0_ratio: float
    start_site_status: str
    start_site_reason: str
    start_interval_orf_ids: tuple[str, ...]
    frame_margin: float
    selection_policy: str
    selection_reason: str
    canonical_anchor_orf: str
    longest_candidate_orf: str
    prior_override_status: str
    leading_support_sample_count: int
    leading_negative_sample_count: int
    noncanonical_extension_support_sample_count: int
    noncanonical_extension_density_ratio: float
    pooled_extent_rpf_sum: float
    pooled_extent_frame0_ratio: float
    pooled_start_rpf_sum: float
    pooled_end_rpf_sum: float
    pooled_extent_sample_count: int


@dataclass(frozen=True, slots=True)
class EngineResult:
    """Store formal output paths and project counts."""

    master_output: str
    reliable_output: str
    reliable_genepred_output: str
    summary_output: str
    total_families: int
    reliable_families: int
    uncertain_families: int
    no_evidence_families: int
    reliable_smorfs: int
    sample_count: int
    effective_workers: int


@dataclass(frozen=True, slots=True)
class ChromosomeTask:
    """Store a resumable chromosome evidence task."""

    chromosome: str
    database_path: str
    work_directory: str
    output_prefix: str
    tracks: tuple[DensityTrack, ...]
    thresholds: tuple[Thresholds, ...]
    config: EngineConfig
    group_names: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ChromosomeResult:
    """Store one chromosome shard and count summary."""

    chromosome: str
    part_path: str
    count_path: str
    total_families: int
    reliable_families: int
    uncertain_families: int
    no_evidence_families: int
    reliable_smorfs: int
    invalid_families: int


@dataclass(slots=True)
class _FamilyAccumulator:
    """Accumulate family translation and candidate-level extent evidence."""

    any_signal: bool
    any_evidence: bool
    support_count: int
    high_count: int
    group_mask: int
    group_support_counts: dict[int, int]
    group_high_counts: dict[int, int]
    supporting_samples: list[str]
    best_level: str
    best_score: float
    best_sample: str
    best_features: SegmentFeatures | None
    candidate_full_support_counts: list[int]
    candidate_extension_support_counts: list[int]
    candidate_support_samples: list[list[str]]
    candidate_bin_support_counts: list[list[int]]
    candidate_bin_support_samples: list[list[str]]
    candidate_cohort_support_samples: list[list[str]]
    candidate_bin_pooled_stats: list[list[list[float]]]
    candidate_boundary_pooled_stats: list[list[list[float]]]
    candidate_best_full_features: list[SegmentFeatures | None]
    candidate_best_extension_features: list[SegmentFeatures | None]
    candidate_pooled_frame_density: list[list[float]]
    candidate_leading_support_counts: list[int]
    candidate_leading_negative_counts: list[int]
    candidate_noncanonical_extension_support_counts: list[int]
    candidate_best_noncanonical_density_ratios: list[float]



class EvidenceEngineError(RuntimeError):
    """Represent a stage-specific evidence-engine failure."""


class _StageLogger:
    """Report stage messages to the active console or nohup stream."""

    def write(self, message: str) -> None:
        message_print(str(message))

    def error(self, stage: str, error: BaseException) -> None:
        warning_print(f"ERROR stage={stage}: {error}")


def _smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain or gzip-compressed text."""
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8", newline="")
    return file_path.open(mode, encoding="utf-8", newline="")


def _split_ints(value: Any) -> tuple[int, ...]:
    """Parse comma-separated integer coordinates."""
    text = str(value).strip().rstrip(",")
    if not text:
        return ()
    return tuple(int(item) for item in text.split(",") if item != "")


def _coding_nt_length(record: Mapping[str, Any]) -> int:
    """Return coding length with a complete terminal stop removed."""
    length = int(record["nt_length"])
    stop = str(record.get("stop_codon", "")).upper().replace("U", "T")
    completeness = str(record.get("completeness", "complete")).lower()
    complete = completeness in {"", "complete", "cmpl", "full", "true", "yes"}
    if complete and stop in STOP_CODONS and length >= 6:
        length -= 3
    return max(0, length - length % 3)


def _file_signature(path: str | Path) -> dict[str, Any]:
    """Return a stable file signature."""
    resolved = Path(path).expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _run_signature(args: object, tracks: Sequence[DensityTrack]) -> str:
    """Hash all material inputs and analysis options."""
    files = {
        "family_table": _file_signature(args.family_table),
        "family_members": _file_signature(args.family_members),
        "orf_source": _file_signature(args.orf_source),
        "orf_genepred": _file_signature(args.orf_genepred),
        "tracks": [_file_signature(track.path) for track in tracks],
    }
    options = {
        key: value
        for key, value in sorted(vars(args).items())
        if key != "threads"
        and isinstance(value, (str, int, float, bool, type(None)))
    }
    payload = json.dumps(
        {
            "schema": SCHEMA_VERSION,
            "files": files,
            "options": options,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _safe_name(value: str) -> str:
    """Return a deterministic filesystem-safe label."""
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:10]
    cleaned = "".join(character if character.isalnum() else "_" for character in value)
    cleaned = cleaned.strip("_")[:80] or "item"
    return f"{cleaned}.{digest}"


def _read_memory_limit() -> int | None:
    """Return the effective cgroup or host memory limit."""
    values: list[int] = []
    for path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            text = path.read_text(encoding="utf-8").strip()
        except OSError:
            continue
        if text and text != "max":
            try:
                value = int(text)
            except ValueError:
                continue
            if 0 < value < 1 << 60:
                values.append(value)
    try:
        with Path("/proc/meminfo").open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    values.append(int(line.split()[1]) * 1024)
                    break
    except OSError:
        pass
    return min(values) if values else None


def _effective_workers(requested: int, tasks: int) -> int:
    """Cap chromosome workers by memory and task count."""
    workers = max(1, min(int(requested), int(tasks), MAX_WORKERS))
    memory = _read_memory_limit()
    if memory is not None:
        workers = min(
            workers,
            max(1, int(memory * 0.70 // MEMORY_PER_WORKER_BYTES)),
        )
    return max(1, workers)


def _parse_density_list(
    path: str | Path,
    group_column: str,
) -> list[DensityTrack]:
    """Parse a precise header-based density design table.

    The design must contain sample, strand, one path column, and the biological
    replicate group column selected by ``--group``. Density format is inferred
    from each data-file path.
    """
    design_path = Path(path).expanduser().resolve()
    with _smart_open(design_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original_header = tuple(reader.fieldnames or ())
        if not original_header:
            raise ValueError(
                "Density list requires a tab-delimited header."
            )
        header_lookup = {
            str(name).strip().lower(): str(name)
            for name in original_header
        }
        path_name = next(
            (
                header_lookup[name]
                for name in ("path", "file", "density")
                if name in header_lookup
            ),
            None,
        )
        sample_name = next(
            (
                header_lookup[name]
                for name in ("sample", "name")
                if name in header_lookup
            ),
            None,
        )
        if sample_name is None or "strand" not in header_lookup:
            raise ValueError(
                "Density list requires sample/name, strand, and "
                "path/file/density columns."
            )
        group_key = str(group_column).strip().lower()
        if group_key not in {"sample", "name"} and group_key not in header_lookup:
            raise ValueError(
                f"Density-list group column was not found: {group_column}"
            )
        if path_name is None:
            raise ValueError(
                "Density list requires one path column: path, file, or density."
            )

        strand_name = header_lookup["strand"]
        group_name = (
            sample_name
            if group_key in {"sample", "name"}
            else header_lookup[group_key]
        )
        tracks: list[DensityTrack] = []
        for line_number, row in enumerate(reader, start=2):
            sample = str(row.get(sample_name, "")).strip()
            group = str(row.get(group_name, "")).strip()
            strand_text = str(row.get(strand_name, "")).strip().lower()
            density_text = str(row.get(path_name, "")).strip()
            if not sample or not group or not strand_text or not density_text:
                raise ValueError(
                    f"Density-list row {line_number} contains an empty required value."
                )
            strand_alias = {
                "+": "+",
                "plus": "+",
                "forward": "+",
                "-": "-",
                "minus": "-",
                "reverse": "-",
                ".": ".",
                "unstranded": ".",
            }
            if strand_text not in strand_alias:
                raise ValueError(
                    f"Density-list row {line_number} has invalid strand: {strand_text}"
                )
            density_path = Path(density_text).expanduser()
            if not density_path.is_absolute():
                density_path = design_path.parent / density_path
            density_path = density_path.resolve()
            if not density_path.is_file():
                raise FileNotFoundError(density_path)
            tracks.append(
                DensityTrack(
                    sample=sample,
                    group=group,
                    strand=strand_alias[strand_text],
                    path=str(density_path),
                    file_format=infer_density_format(density_path),
                )
            )
    if not tracks:
        raise ValueError("Density list contains no data rows.")
    return tracks



def read_density_tracks(args: object) -> list[DensityTrack]:
    """Read the required density design and validate sample/strand uniqueness."""
    density_list = getattr(args, "density_list", None)
    if not density_list:
        raise ValueError("--density-list is required.")
    tracks = _parse_density_list(
        density_list,
        str(getattr(args, "group_column", "group")),
    )

    seen: set[tuple[str, str]] = set()
    sample_groups: dict[str, str] = {}
    for track in tracks:
        key = (track.sample, track.strand)
        if key in seen:
            raise ValueError(
                "Only one density track per sample/strand is supported: "
                f"{track.sample} {track.strand}"
            )
        seen.add(key)
        previous_group = sample_groups.setdefault(track.sample, track.group)
        if previous_group != track.group:
            raise ValueError(
                f"Sample {track.sample} is assigned to multiple groups: "
                f"{previous_group}, {track.group}"
            )
    return tracks



def _connect_database(path: Path) -> sqlite3.Connection:
    """Open the disk index with memory-bounded pragmas."""
    connection = sqlite3.connect(path, timeout=120.0)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-262144")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    return connection


def _create_database_schema(connection: sqlite3.Connection) -> None:
    """Create the family and representative index schema."""
    connection.executescript(
        """
        CREATE TABLE families (
            family_id TEXT PRIMARY KEY,
            gene_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            strand TEXT NOT NULL,
            category TEXT NOT NULL,
            family_type TEXT NOT NULL,
            family_size INTEGER NOT NULL,
            structural_primary TEXT NOT NULL
        );
        CREATE TABLE representatives (
            family_id TEXT NOT NULL,
            orf_id TEXT PRIMARY KEY,
            is_primary INTEGER NOT NULL DEFAULT 0
        );
        CREATE INDEX representatives_family_idx
            ON representatives(family_id);
        CREATE TABLE geometry (
            orf_id TEXT PRIMARY KEY,
            transcript_id TEXT NOT NULL,
            category TEXT NOT NULL,
            start_codon TEXT NOT NULL,
            stop_codon TEXT NOT NULL,
            completeness TEXT NOT NULL,
            chrom TEXT NOT NULL,
            strand TEXT NOT NULL,
            nt_length INTEGER NOT NULL,
            aa_length INTEGER NOT NULL,
            exon_starts TEXT NOT NULL,
            exon_ends TEXT NOT NULL
        );
        CREATE INDEX geometry_chrom_idx ON geometry(chrom);
        """
    )


def _required_columns(header: Sequence[str], required: Sequence[str], name: str) -> None:
    """Validate a tabular header."""
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"{name} is missing column(s): {', '.join(missing)}")


def build_family_database(
    family_table: str | Path,
    family_members: str | Path,
    orf_source: str | Path,
    database_path: Path,
    signature: str,
    logger: _StageLogger,
) -> None:
    """Build or reuse a disk-backed family/representative index."""
    marker = database_path.with_suffix(".complete.json")
    if database_path.is_file() and marker.is_file():
        state = json.loads(marker.read_text(encoding="utf-8"))
        if state.get("signature") == signature:
            logger.write("Reuse family SQLite index.")
            return

    database_path.unlink(missing_ok=True)
    marker.unlink(missing_ok=True)
    connection = _connect_database(database_path)
    try:
        _create_database_schema(connection)

        logger.write("Index family primary records.")
        with _smart_open(family_table, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            _required_columns(
                header,
                (
                    "family_id", "gene_id", "chrom", "strand", "category",
                    "family_type", "family_size", "orf_id", "transcript_id",
                    "start_codon", "stop_codon", "completeness", "nt_length",
                    "aa_length", "exon_starts", "exon_ends",
                ),
                "family table",
            )
            family_rows: list[tuple[Any, ...]] = []
            rep_rows: list[tuple[Any, ...]] = []
            geometry_rows: list[tuple[Any, ...]] = []
            for number, row in enumerate(reader, start=1):
                family_rows.append(
                    (
                        row["family_id"], row["gene_id"], row["chrom"],
                        row["strand"], row["category"], row["family_type"],
                        int(row["family_size"]), row["orf_id"],
                    )
                )
                rep_rows.append((row["family_id"], row["orf_id"], 1))
                geometry_rows.append(
                    (
                        row["orf_id"], row["transcript_id"], row["category"],
                        row["start_codon"], row["stop_codon"],
                        row["completeness"], row["chrom"], row["strand"],
                        int(row["nt_length"]), int(row["aa_length"]),
                        row["exon_starts"], row["exon_ends"],
                    )
                )
                if len(family_rows) >= SQL_BATCH_SIZE:
                    connection.executemany(
                        "INSERT INTO families VALUES (?,?,?,?,?,?,?,?)",
                        family_rows,
                    )
                    connection.executemany(
                        "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                        rep_rows,
                    )
                    connection.executemany(
                        "INSERT OR REPLACE INTO geometry VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                        geometry_rows,
                    )
                    connection.commit()
                    family_rows.clear(); rep_rows.clear(); geometry_rows.clear()
                    progress_print(f"family index: {number:,}")
            if family_rows:
                connection.executemany(
                    "INSERT INTO families VALUES (?,?,?,?,?,?,?,?)",
                    family_rows,
                )
                connection.executemany(
                    "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                    rep_rows,
                )
                connection.executemany(
                    "INSERT OR REPLACE INTO geometry VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    geometry_rows,
                )
                connection.commit()

        logger.write("Index alternative-start representative relationships.")
        with _smart_open(family_members, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            _required_columns(
                header,
                ("family_id", "representative_orf_id"),
                "family members",
            )
            rows: list[tuple[str, str, int]] = []
            for number, row in enumerate(reader, start=1):
                representative = row["representative_orf_id"]
                if not representative:
                    continue
                rows.append(
                    (
                        row["family_id"],
                        representative,
                        int(representative == row.get("primary_orf_id", "")),
                    )
                )
                if len(rows) >= SQL_BATCH_SIZE:
                    connection.executemany(
                        "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                        rows,
                    )
                    connection.commit(); rows.clear()
                    progress_print(f"member index: {number:,}")
            if rows:
                connection.executemany(
                    "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                    rows,
                )
                connection.commit()

        logger.write("Read source ORF table and fill representative geometry.")
        with _smart_open(orf_source, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            _required_columns(
                header,
                (
                    "orf_id", "transcript_id", "category", "start_codon",
                    "stop_codon", "completeness", "chrom", "strand",
                    "nt_length", "aa_length", "exon_starts", "exon_ends",
                ),
                "source ORF table",
            )
            connection.execute(
                """
                CREATE TEMP TABLE source_batch (
                    orf_id TEXT, transcript_id TEXT, category TEXT,
                    start_codon TEXT, stop_codon TEXT, completeness TEXT,
                    chrom TEXT, strand TEXT, nt_length INTEGER,
                    aa_length INTEGER, exon_starts TEXT, exon_ends TEXT
                )
                """
            )
            batch: list[tuple[Any, ...]] = []
            for number, row in enumerate(reader, start=1):
                batch.append(
                    (
                        row["orf_id"], row["transcript_id"], row["category"],
                        row["start_codon"], row["stop_codon"],
                        row["completeness"], row["chrom"], row["strand"],
                        int(row["nt_length"]), int(row["aa_length"]),
                        row["exon_starts"], row["exon_ends"],
                    )
                )
                if len(batch) >= SQL_BATCH_SIZE:
                    connection.executemany(
                        "INSERT INTO source_batch VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                        batch,
                    )
                    connection.execute(
                        """
                        INSERT OR REPLACE INTO geometry
                        SELECT b.* FROM source_batch b
                        INNER JOIN representatives r ON r.orf_id=b.orf_id
                        """
                    )
                    connection.execute("DELETE FROM source_batch")
                    connection.commit(); batch.clear()
                    progress_print(f"source ORFs: {number:,}")
            if batch:
                connection.executemany(
                    "INSERT INTO source_batch VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                connection.execute(
                    """
                    INSERT OR REPLACE INTO geometry
                    SELECT b.* FROM source_batch b
                    INNER JOIN representatives r ON r.orf_id=b.orf_id
                    """
                )
                connection.execute("DELETE FROM source_batch")
                connection.commit()

        missing = connection.execute(
            """
            SELECT COUNT(*) FROM representatives r
            LEFT JOIN geometry g ON g.orf_id=r.orf_id
            WHERE g.orf_id IS NULL
            """
        ).fetchone()[0]
        if missing:
            raise ValueError(
                f"Missing source coordinates for {missing:,} family representatives."
            )
        connection.execute("CREATE INDEX geometry_orf_idx ON geometry(orf_id)")
        connection.execute("ANALYZE")
        connection.commit()
        counts = {
            "families": connection.execute("SELECT COUNT(*) FROM families").fetchone()[0],
            "representatives": connection.execute(
                "SELECT COUNT(*) FROM representatives"
            ).fetchone()[0],
        }
        marker.write_text(
            json.dumps({"signature": signature, **counts}, indent=2),
            encoding="utf-8",
        )
        logger.write(
            "Family index completed: families={families:,}, representatives={representatives:,}.".format(
                **counts
            )
        )
    finally:
        connection.close()


def _density_cache_track_directory(work_directory: Path, track: DensityTrack) -> Path:
    """Return one track cache directory."""
    return work_directory / "density" / _safe_name(
        f"{track.sample}|{track.strand}|{Path(track.path).resolve()}"
    )


def build_density_cache(
    track: DensityTrack,
    work_directory: Path,
    signature: str,
) -> None:
    """Parse one density file once and save chromosome NumPy arrays."""
    directory = _density_cache_track_directory(work_directory, track)
    marker = directory / "complete.json"
    if marker.is_file():
        state = json.loads(marker.read_text(encoding="utf-8"))
        if state.get("signature") == signature:
            return
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True, exist_ok=True)
    sort_directory = directory / "sort"
    prepared = prepare_density(
        path=track.path,
        file_format=track.file_format,
        work_directory=sort_directory,
    )
    chrom_index: dict[str, str] = {}
    try:
        for number, density in enumerate(
            iter_density_chromosomes(
                prepared.path,
                "bedgraph" if prepared.was_sorted else track.file_format,
            ),
            start=1,
        ):
            file_name = _safe_name(density.chrom) + ".npy"
            output_path = directory / file_name
            array = np.empty(
                density.starts.size,
                dtype=[("start", "<i8"), ("end", "<i8"), ("value", "<f4")],
            )
            array["start"] = density.starts
            array["end"] = density.ends
            array["value"] = density.values
            np.save(output_path, array, allow_pickle=False)
            chrom_index[density.chrom] = file_name
            progress_print(
                f"density cache {track.sample}/{track.strand}: chromosome {number}"
            )
    finally:
        prepared.cleanup()
        shutil.rmtree(sort_directory, ignore_errors=True)
    marker.write_text(
        json.dumps(
            {
                "signature": signature,
                "sample": track.sample,
                "group": track.group,
                "strand": track.strand,
                "chromosomes": chrom_index,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def _load_cached_density(
    work_directory: Path,
    track: DensityTrack,
    chrom: str,
) -> ChromDensity | None:
    """Load one chromosome density as memory-mapped arrays."""
    directory = _density_cache_track_directory(work_directory, track)
    marker_path = directory / "complete.json"
    if not marker_path.is_file():
        raise FileNotFoundError(marker_path)
    state = json.loads(marker_path.read_text(encoding="utf-8"))
    file_name = state.get("chromosomes", {}).get(chrom)
    if not file_name:
        return None
    array = np.load(directory / file_name, mmap_mode="r", allow_pickle=False)
    starts = array["start"]
    ends = array["end"]
    values = array["value"]
    return ChromDensity(
        chrom=chrom,
        starts=starts,
        ends=ends,
        values=values,
        max_end=int(ends[-1]) if ends.size else 0,
    )


def _oriented_blocks(
    starts: Sequence[int],
    ends: Sequence[int],
    strand: str,
) -> tuple[tuple[int, int], ...]:
    """Return exon blocks in transcript order."""
    blocks = tuple(sorted(zip(starts, ends), key=lambda item: item[0]))
    return blocks if strand == "+" else tuple(reversed(blocks))


def _is_suffix(
    scaffold_starts: Sequence[int],
    scaffold_ends: Sequence[int],
    rep_starts: Sequence[int],
    rep_ends: Sequence[int],
    strand: str,
) -> bool:
    """Return whether a representative is an exact spliced suffix."""
    scaffold = _oriented_blocks(scaffold_starts, scaffold_ends, strand)
    representative = _oriented_blocks(rep_starts, rep_ends, strand)
    if not representative or len(representative) > len(scaffold):
        return False
    tail = scaffold[-len(representative):]
    if len(representative) > 1 and tail[1:] != representative[1:]:
        return False
    long_first = tail[0]
    short_first = representative[0]
    if strand == "+":
        return short_first[1] == long_first[1] and short_first[0] >= long_first[0]
    return short_first[0] == long_first[0] and short_first[1] <= long_first[1]


def _row_to_mapping(row: sqlite3.Row) -> dict[str, Any]:
    return {key: row[key] for key in row.keys()}


def _build_family_geometry(rows: Sequence[sqlite3.Row]) -> FamilyGeometry | InvalidFamily:
    """Compile one family and convert incompatible geometry to a record."""
    first = rows[0]
    base = {
        "family_id": first["family_id"],
        "gene_id": first["gene_id"],
        "chrom": first["chrom"],
        "strand": first["strand"],
        "category": first["family_category"],
        "family_type": first["family_type"],
        "family_size": int(first["family_size"]),
        "structural_primary": first["structural_primary"],
    }
    try:
        parsed: list[dict[str, Any]] = []
        for row in rows:
            record = _row_to_mapping(row)
            starts = _split_ints(record["exon_starts"])
            ends = _split_ints(record["exon_ends"])
            if not starts or len(starts) != len(ends):
                raise ValueError(f"invalid exon blocks for {record['orf_id']}")
            if any(end <= start for start, end in zip(starts, ends)):
                raise ValueError(f"invalid exon interval for {record['orf_id']}")
            record["starts"] = starts
            record["ends"] = ends
            record["coding_nt_length"] = _coding_nt_length(record)
            parsed.append(record)

        scaffold = max(
            parsed,
            key=lambda record: (
                int(record["nt_length"]),
                record["orf_id"] == base["structural_primary"],
            ),
        )
        representatives: list[RepresentativeGeometry] = []
        scaffold_nt_length = int(scaffold["nt_length"])
        for record in parsed:
            difference = scaffold_nt_length - int(record["nt_length"])
            if difference < 0 or difference % 3 != 0:
                raise ValueError(
                    f"phase-incompatible representative {record['orf_id']}"
                )
            if not _is_suffix(
                scaffold["starts"], scaffold["ends"],
                record["starts"], record["ends"], base["strand"],
            ):
                raise ValueError(
                    f"non-suffix representative {record['orf_id']}"
                )
            representatives.append(
                RepresentativeGeometry(
                    orf_id=record["orf_id"],
                    transcript_id=record["transcript_id"],
                    category=record["category"],
                    start_codon=record["start_codon"],
                    nt_length=int(record["nt_length"]),
                    coding_nt_length=int(record["coding_nt_length"]),
                    aa_length=int(record["aa_length"]),
                    starts=record["starts"],
                    ends=record["ends"],
                    start_offset=difference,
                    stop_codon=record["stop_codon"],
                    completeness=record["completeness"],
                )
            )
        representatives.sort(key=lambda item: (item.start_offset, item.orf_id))
        common = max(representatives, key=lambda item: item.start_offset)
        coding_end = int(scaffold["coding_nt_length"])
        if common.start_offset >= coding_end:
            raise ValueError("empty common-body coding segment")
        return FamilyGeometry(
            **base,
            scaffold_orf_id=scaffold["orf_id"],
            scaffold_nt_length=scaffold_nt_length,
            scaffold_coding_nt_length=coding_end,
            scaffold_starts=scaffold["starts"],
            scaffold_ends=scaffold["ends"],
            common_body_orf=common.orf_id,
            common_body_start=common.start_offset,
            common_body_end=coding_end,
            representatives=tuple(representatives),
        )
    except Exception as error:
        return InvalidFamily(**base, reason=str(error))


def _iter_chromosome_family_batches(
    database_path: str | Path,
    chromosome: str,
    batch_size: int = FAMILY_BATCH_SIZE,
) -> Iterator[list[FamilyGeometry | InvalidFamily]]:
    """Stream compiled family geometry from SQLite in bounded batches."""
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    query = """
        SELECT
            f.family_id, f.gene_id, f.chrom, f.strand,
            f.category AS family_category, f.family_type, f.family_size,
            f.structural_primary, r.orf_id, r.is_primary,
            g.transcript_id, g.category, g.start_codon, g.stop_codon,
            g.completeness, g.nt_length, g.aa_length,
            g.exon_starts, g.exon_ends
        FROM families f
        JOIN representatives r ON r.family_id=f.family_id
        JOIN geometry g ON g.orf_id=r.orf_id
        WHERE f.chrom=?
        ORDER BY f.family_id, g.nt_length DESC, r.orf_id
    """
    try:
        current_id: str | None = None
        current_rows: list[sqlite3.Row] = []
        batch: list[FamilyGeometry | InvalidFamily] = []
        for row in connection.execute(query, (chromosome,)):
            family_id = row["family_id"]
            if current_id is None:
                current_id = family_id
            elif family_id != current_id:
                batch.append(_build_family_geometry(current_rows))
                current_rows = []
                current_id = family_id
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
            current_rows.append(row)
        if current_rows:
            batch.append(_build_family_geometry(current_rows))
        if batch:
            yield batch
    finally:
        connection.close()


def _chromosomes(database_path: str | Path) -> list[str]:
    """Return chromosomes in deterministic lexical-natural order."""
    connection = sqlite3.connect(database_path)
    try:
        values = [row[0] for row in connection.execute(
            "SELECT DISTINCT chrom FROM families"
        )]
    finally:
        connection.close()

    def key(value: str) -> tuple[Any, ...]:
        pieces: list[Any] = []
        token = ""
        numeric = value[:1].isdigit()
        for character in value:
            now_numeric = character.isdigit()
            if token and now_numeric != numeric:
                pieces.append(int(token) if numeric else token)
                token = ""
            token += character
            numeric = now_numeric
        if token:
            pieces.append(int(token) if numeric else token)
        return tuple(pieces)

    return sorted(values, key=key)


def _track_for_family(
    tracks_by_sample: Mapping[str, Mapping[str, DensityTrack]],
    sample: str,
    strand: str,
) -> DensityTrack | None:
    """Select a strand-specific or unstranded track."""
    mapping = tracks_by_sample[sample]
    return mapping.get(strand) or mapping.get(".")


def _extract_sparse_signals(
    families: Sequence[FamilyGeometry],
    density: ChromDensity | None,
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Map one chromosome density to family scaffold transcript offsets."""
    if density is None or density.starts.size == 0 or not families:
        return {}

    block_starts: list[int] = []
    block_ends: list[int] = []
    block_offsets: list[int] = []
    block_families: list[int] = []
    block_strands: list[str] = []
    for family_index, family in enumerate(families):
        blocks = _oriented_blocks(
            family.scaffold_starts,
            family.scaffold_ends,
            family.strand,
        )
        offset = 0
        for start, end in blocks:
            block_starts.append(start)
            block_ends.append(end)
            block_offsets.append(offset)
            block_families.append(family_index)
            block_strands.append(family.strand)
            offset += end - start

    starts_array = np.asarray(block_starts, dtype=np.int64)
    ends_array = np.asarray(block_ends, dtype=np.int64)
    left = np.searchsorted(density.ends, starts_array, side="right")
    right = np.searchsorted(density.starts, ends_array, side="left")

    positions: dict[int, list[int]] = defaultdict(list)
    values: dict[int, list[float]] = defaultdict(list)
    for block_index in np.flatnonzero(right > left):
        family_index = block_families[int(block_index)]
        block_start = block_starts[int(block_index)]
        block_end = block_ends[int(block_index)]
        block_offset = block_offsets[int(block_index)]
        strand = block_strands[int(block_index)]
        for density_index in range(int(left[block_index]), int(right[block_index])):
            value = abs(float(density.values[density_index]))
            if value <= 0 or not math.isfinite(value):
                continue
            overlap_start = max(block_start, int(density.starts[density_index]))
            overlap_end = min(block_end, int(density.ends[density_index]))
            if overlap_end <= overlap_start:
                continue
            if strand == "+":
                tx_start = block_offset + overlap_start - block_start
                tx_end = block_offset + overlap_end - block_start
            else:
                tx_start = block_offset + block_end - overlap_end
                tx_end = block_offset + block_end - overlap_start
            span = tx_end - tx_start
            if span <= 0:
                continue
            positions[family_index].extend(range(tx_start, tx_end))
            values[family_index].extend([value] * span)

    output: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for family_index, family_positions in positions.items():
        order = np.argsort(np.asarray(family_positions, dtype=np.int64), kind="stable")
        output[family_index] = (
            np.asarray(family_positions, dtype=np.int64)[order],
            np.asarray(values[family_index], dtype=np.float64)[order],
        )
    return output


def _empty_segment_features() -> SegmentFeatures:
    """Return an all-zero segment feature record."""
    return SegmentFeatures(
        0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0,
    )


def segment_features(
    positions: np.ndarray,
    values: np.ndarray,
    start: int,
    end: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> SegmentFeatures:
    """Calculate sufficient statistics for one scaffold segment."""
    return ScaffoldProfile(
        positions,
        values,
        max(0, int(end)),
    ).region_features(start, end, thresholds, config)

def classify_sample(
    features: SegmentFeatures,
    codon_count: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> tuple[str, str]:
    """Classify one family common body in one independent sample."""
    if features.rpf_sum <= 0:
        return "NoEvidence", "zero_common_body_rpf"
    core_pass = (
        features.rpf_sum >= thresholds.min_rpf_sum
        and features.rpf_per_codon >= thresholds.min_rpf_per_codon
        and features.covered_codon >= thresholds.min_covered_codon
        and features.coverage_ratio >= thresholds.min_coverage_ratio
    )
    if not core_pass:
        return "LowConfidence", "abundance_or_coverage_below_threshold"
    if features.frame0_ratio < thresholds.moderate_periodicity:
        return "LowConfidence", "periodicity_below_calibrated_threshold"
    if features.localized_only:
        return "LowConfidence", "localized_only_long_orf_signal"

    if codon_count <= config.short_max_codons:
        support = True
    elif codon_count < config.long_min_codons:
        support = (
            features.supported_windows >= 1
            or features.coverage_ratio >= max(0.20, thresholds.min_coverage_ratio * 1.5)
        )
    else:
        support = (
            features.distributed_windows >= config.min_supported_windows
            or features.signal_span >= config.min_signal_span
            or (features.start_rpf > 0 and features.end_rpf > 0)
        )
    if not support:
        return "LowConfidence", "whole_orf_distribution_not_supported"

    high = (
        features.frame0_ratio >= thresholds.strong_periodicity
        and features.coverage_ratio >= max(thresholds.min_coverage_ratio, 0.20)
        and (
            codon_count <= config.short_max_codons
            or features.distributed_windows >= config.min_supported_windows
            or features.signal_span >= config.min_signal_span
        )
    )
    return (
        ("HighConfidence", "strong_replicable_common_body_evidence")
        if high
        else ("MediumConfidence", "moderate_common_body_evidence")
    )


def classify_extension(
    features: SegmentFeatures,
    extension_codons: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> bool:
    """Return whether one cumulative start-specific extension is supported."""
    if extension_codons <= 0:
        return True
    minimum_covered = min(
        max(1, thresholds.min_window_covered),
        max(1, extension_codons),
    )
    core_support = (
        features.rpf_sum >= max(1.0, thresholds.min_window_rpf * 0.5)
        and features.covered_codon >= minimum_covered
        and features.frame0_ratio >= thresholds.moderate_periodicity
    )
    if not core_support or features.localized_only:
        return False
    if extension_codons < config.long_min_codons:
        return True
    return (
        features.distributed_windows >= 1
        or features.signal_span >= max(0.20, config.min_signal_span * 0.75)
    )


def classify_candidate_extent(
    full_features: SegmentFeatures,
    extension_features: SegmentFeatures,
    full_codons: int,
    extension_codons: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> tuple[bool, bool]:
    """Classify one candidate's full range and cumulative extension."""
    full_level, _reason = classify_sample(
        full_features,
        full_codons,
        thresholds,
        config,
    )
    full_supported = LEVEL_RANK[full_level] >= LEVEL_RANK["MediumConfidence"]
    extension_supported = classify_extension(
        extension_features,
        extension_codons,
        thresholds,
        config,
    )
    return full_supported, extension_supported


def _canonical_anchor_index(family: FamilyGeometry) -> int | None:
    """Return the longest canonical ATG representative, if present."""
    canonical = [
        index
        for index, representative in enumerate(family.representatives)
        if representative.start_codon.upper() == "ATG"
    ]
    return min(
        canonical,
        key=lambda index: (
            family.representatives[index].start_offset,
            family.representatives[index].orf_id,
        ),
        default=None,
    )


def _leading_region_end(
    family: FamilyGeometry,
    candidate_index: int,
    canonical_index: int | None,
    config: EngineConfig,
) -> int:
    """Return the diagnostic N-terminal boundary for one candidate."""
    representative = family.representatives[candidate_index]
    boundary = min(
        family.common_body_end,
        representative.start_offset + config.leading_window_codons * 3,
    )
    if canonical_index is not None:
        canonical_start = family.representatives[canonical_index].start_offset
        if representative.start_offset < canonical_start:
            boundary = min(boundary, canonical_start)
    return max(representative.start_offset, boundary)


def _leading_region_is_supported(
    features: SegmentFeatures,
    codon_count: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> bool:
    """Return whether a candidate-specific N terminus contains translation."""
    if codon_count < config.min_exclusion_codons:
        return False
    return (
        features.rpf_sum >= thresholds.min_window_rpf
        and features.covered_codon >= thresholds.min_window_covered
        and features.frame0_ratio >= thresholds.moderate_periodicity
        and not features.localized_only
    )


def _leading_region_is_negative(
    leading: SegmentFeatures,
    common_body: SegmentFeatures,
    codon_count: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> bool:
    """Return whether an observable N terminus is silent relative to its body."""
    if codon_count < config.min_exclusion_codons:
        return False
    if (
        common_body.rpf_sum < thresholds.min_rpf_sum
        or common_body.rpf_per_codon < thresholds.min_rpf_per_codon
        or common_body.frame0_ratio < thresholds.moderate_periodicity
    ):
        return False
    relative_density = leading.rpf_per_codon / max(
        common_body.rpf_per_codon,
        1e-12,
    )
    weak_abundance = (
        leading.rpf_sum < thresholds.min_window_rpf
        or leading.covered_codon < thresholds.min_window_covered
    )
    weak_phase = leading.frame0_ratio < thresholds.moderate_periodicity
    return (
        relative_density <= config.silent_extension_density_ratio
        and (weak_abundance or weak_phase)
    )


def _noncanonical_extension_is_strong(
    extension: SegmentFeatures,
    canonical_body: SegmentFeatures,
    extension_codons: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> bool:
    """Require strong exclusive evidence before replacing a canonical ATG."""
    if extension_codons < config.noncanonical_min_exclusive_codons:
        return False
    if not (
        canonical_body.rpf_sum >= thresholds.min_rpf_sum
        and canonical_body.rpf_per_codon
        >= thresholds.min_rpf_per_codon
        and canonical_body.frame0_ratio
        >= thresholds.moderate_periodicity
    ):
        return False
    density_floor = max(
        thresholds.min_rpf_per_codon,
        canonical_body.rpf_per_codon
        * config.noncanonical_override_density_ratio,
    )
    if not (
        extension.rpf_sum >= max(
            thresholds.min_rpf_sum,
            thresholds.min_window_rpf,
        )
        and extension.rpf_per_codon >= density_floor
        and extension.covered_codon >= thresholds.min_window_covered
        and extension.coverage_ratio >= max(
            thresholds.min_coverage_ratio,
            config.noncanonical_override_min_coverage_ratio,
        )
        and extension.frame0_ratio >= thresholds.strong_periodicity
        and (
            extension.frame0_ratio
            - max(
                extension.frame1_density,
                extension.frame2_density,
            )
            / max(
                extension.frame0_density
                + extension.frame1_density
                + extension.frame2_density,
                1e-12,
            )
            >= config.noncanonical_override_min_frame_margin
        )
        and not extension.localized_only
    ):
        return False
    if extension_codons < config.long_min_codons:
        return True
    return (
        extension.distributed_windows >= 1
        or extension.signal_span >= max(0.20, config.min_signal_span * 0.75)
    )


def _candidate_bin_ranges(
    start: int,
    end: int,
    bin_count: int,
) -> tuple[tuple[int, int], ...]:
    """Split one candidate coding region into codon-aligned extent bins."""
    codons = max(0, (int(end) - int(start)) // 3)
    if codons <= 0:
        return ()
    effective_bins = min(max(1, int(bin_count)), codons)
    boundaries = [
        (index * codons) // effective_bins
        for index in range(effective_bins + 1)
    ]
    return tuple(
        (start + boundaries[index] * 3, start + boundaries[index + 1] * 3)
        for index in range(effective_bins)
        if boundaries[index + 1] > boundaries[index]
    )


def _bin_is_supported(
    features: SegmentFeatures,
    thresholds: Thresholds,
) -> bool:
    """Return whether one longitudinal extent bin contains phased signal."""
    return (
        features.rpf_sum >= max(1.0, thresholds.min_window_rpf * 0.5)
        and features.covered_codon
        >= max(1, math.ceil(thresholds.min_window_covered * 0.5))
        and features.frame0_ratio >= thresholds.moderate_periodicity
    )


def _accumulate_pooled_segment(
    pooled: list[float],
    features: SegmentFeatures,
) -> None:
    """Accumulate raw segment sufficient statistics across samples."""
    pooled[0] += features.rpf_sum
    pooled[1] += features.covered_codon
    pooled[2] += features.frame0_density
    pooled[3] += features.frame1_density
    pooled[4] += features.frame2_density
    if features.rpf_sum > 0:
        pooled[5] += 1


def _pooled_segment_is_supported(
    pooled: Sequence[float],
    config: EngineConfig,
) -> bool:
    """Return whether pooled subthreshold observations form phased evidence."""
    if len(pooled) < 6 or pooled[5] < 1:
        return False
    frames = pooled[2:5]
    total = float(sum(frames))
    if total <= 0:
        return False
    return (
        pooled[0] >= max(1.0, config.min_window_rpf)
        and pooled[1] >= max(2, math.ceil(config.min_window_covered * 0.50))
        and frames[0] / total >= config.moderate_periodicity
        and _frame_margin(frames) >= config.min_frame_margin
    )


def _candidate_is_complete(representative: RepresentativeGeometry) -> bool:
    """Return whether scanner geometry describes a complete start-to-stop ORF."""
    completeness = representative.completeness.strip().lower()
    complete = completeness in {
        "", "complete", "cmpl", "full", "true", "yes",
    }
    return complete and representative.stop_codon.upper() in STOP_CODONS


def _candidate_supported_bin_indices(
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    config: EngineConfig,
) -> list[int]:
    """Return bins supported either directly or by pooled raw statistics."""
    binary = accumulator.candidate_bin_support_counts[candidate_index]
    pooled = accumulator.candidate_bin_pooled_stats[candidate_index]
    return [
        index
        for index in range(max(len(binary), len(pooled)))
        if (
            (index < len(binary) and binary[index] > 0)
            or (
                index < len(pooled)
                and _pooled_segment_is_supported(pooled[index], config)
            )
        )
    ]


def _candidate_boundary_support(
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    config: EngineConfig,
) -> tuple[bool, bool]:
    """Return pooled start- and stop-window support for one candidate."""
    boundary = accumulator.candidate_boundary_pooled_stats[candidate_index]
    start_supported = bool(
        boundary
        and _pooled_segment_is_supported(boundary[0], config)
    )
    end_supported = bool(
        len(boundary) > 1
        and _pooled_segment_is_supported(boundary[1], config)
    )
    binary = accumulator.candidate_bin_support_counts[candidate_index]
    if binary:
        start_supported = start_supported or binary[0] > 0
        end_supported = end_supported or binary[-1] > 0
    return start_supported, end_supported


def _candidate_cohort_sample_count(
    accumulator: _FamilyAccumulator,
    candidate_index: int,
) -> int:
    """Return independent samples contributing candidate-phased evidence."""
    samples = set(accumulator.candidate_cohort_support_samples[candidate_index])
    samples.update(accumulator.candidate_bin_support_samples[candidate_index])
    samples.update(accumulator.candidate_support_samples[candidate_index])
    return len(samples)


def _candidate_cohort_extent_status(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    config: EngineConfig,
) -> str:
    """Classify long-ORF extent assembled across independent samples.

    A long ORF may be sparsely covered in each sample while distinct samples
    support different parts of the same coding range. The cohort model accepts
    that distributed evidence, but it never treats a middle-only signal as
    proof of the complete parent ORF.
    """
    representative = family.representatives[candidate_index]
    codon_count = representative.coding_nt_length // 3
    bin_counts = accumulator.candidate_bin_support_counts[candidate_index]
    if (
        codon_count < config.long_min_codons
        or not bin_counts
        or not _candidate_is_complete(representative)
    ):
        return ""
    sample_count = _candidate_cohort_sample_count(
        accumulator,
        candidate_index,
    )
    if sample_count < config.reliable_sample:
        return ""
    supported_bins = _candidate_supported_bin_indices(
        accumulator,
        candidate_index,
        config,
    )
    minimum_bins = min(config.min_extent_bins, len(bin_counts))
    if len(supported_bins) < minimum_bins:
        return ""
    segmented_frames = _candidate_segmented_frame_density(
        accumulator,
        candidate_index,
        supported_bins,
    )
    if _frame_margin(segmented_frames) < config.min_frame_margin:
        return ""
    if segmented_frames[0] < max(
        config.min_window_rpf * minimum_bins,
        config.min_rpf_sum * config.reliable_sample,
    ):
        return ""
    start_supported, end_supported = _candidate_boundary_support(
        accumulator,
        candidate_index,
        config,
    )
    if supported_bins[0] != 0 or not start_supported:
        return ""

    span_fraction = (
        supported_bins[-1] - supported_bins[0] + 1
    ) / len(bin_counts)
    if (
        supported_bins[-1] == len(bin_counts) - 1
        and end_supported
    ):
        return "Supported"
    downstream_start = max(1, math.floor(len(bin_counts) * 0.60))
    if (
        span_fraction >= max(0.60, config.min_signal_span)
        and any(index >= downstream_start for index in supported_bins)
    ):
        return "Compatible"
    return ""


def _candidate_segmented_frame_density(
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    supported_bins: Sequence[int] | None = None,
) -> tuple[float, float, float]:
    """Return pooled frame counts from candidate-phased supported bins.

    Strong signal from an independently translated overlapping frame may
    dominate whole-ORF counts. Restricting the decision statistic to bins that
    independently support the candidate prevents that unrelated signal from
    vetoing a distributed complete ORF while retaining an explicit whole-ORF
    phase audit in the output.
    """
    pooled = accumulator.candidate_bin_pooled_stats[candidate_index]
    indices = (
        range(len(pooled))
        if supported_bins is None
        else supported_bins
    )
    frames = [0.0, 0.0, 0.0]
    for index in indices:
        if index < 0 or index >= len(pooled):
            continue
        frames[0] += float(pooled[index][2])
        frames[1] += float(pooled[index][3])
        frames[2] += float(pooled[index][4])
    if sum(frames) <= 0:
        fallback = accumulator.candidate_pooled_frame_density[candidate_index]
        return float(fallback[0]), float(fallback[1]), float(fallback[2])
    return frames[0], frames[1], frames[2]


def _manual_thresholds(
    sample: str,
    config: EngineConfig,
) -> Thresholds:
    """Return the eight user-defined manual evidence thresholds."""
    return Thresholds(
        sample=sample,
        source="manual",
        control_count=0,
        min_rpf_sum=config.min_rpf_sum,
        min_rpf_per_codon=config.min_rpf_per_codon,
        min_covered_codon=config.min_covered_codon,
        min_coverage_ratio=config.min_coverage_ratio,
        moderate_periodicity=config.moderate_periodicity,
        strong_periodicity=config.strong_periodicity,
        min_window_rpf=config.min_window_rpf,
        min_window_covered=config.min_window_covered,
    )




def _normalized_length_category(category: str) -> str | None:
    """Map target smORF categories to uORF, dORF, or lncORF."""
    value = str(category).strip().lower()
    for target in TARGET_LENGTH_CATEGORIES:
        if target in value:
            return target
    return None


def resolve_length_models(
    database_path: Path,
    config: EngineConfig,
    logger: _StageLogger,
) -> EngineConfig:
    """Resolve robust short and long cutoffs from target-smORF distributions.

    Equal-tertile cutoffs are inappropriate for smORFs because their length
    distribution is strongly right-skewed. The raw short boundary is estimated
    from the median category-specific lower quartile, and the raw long boundary
    from the median category-specific upper decile. Conservative floors prevent
    project-specific distributions from switching evidence models at implausibly
    small values. These thresholds classify evidence models; they do not filter
    ORFs by length.
    """
    by_category: dict[str, list[int]] = defaultdict(list)
    connection = sqlite3.connect(database_path)
    try:
        query = """
            SELECT f.category, g.nt_length, g.stop_codon, g.completeness
            FROM families f
            JOIN geometry g ON g.orf_id=f.structural_primary
        """
        for category, nt_length, stop_codon, completeness in connection.execute(query):
            normalized = _normalized_length_category(category)
            if normalized is None:
                continue
            coding = _coding_nt_length(
                {
                    "nt_length": nt_length,
                    "stop_codon": stop_codon,
                    "completeness": completeness,
                }
            )
            codons = coding // 3
            if codons > 0:
                by_category[normalized].append(codons)
    finally:
        connection.close()

    short_values: list[float] = []
    long_values: list[float] = []
    for category in sorted(TARGET_LENGTH_CATEGORIES):
        lengths = by_category.get(category, [])
        if not lengths:
            continue
        short_values.append(float(np.quantile(lengths, 0.25)))
        long_values.append(float(np.quantile(lengths, 0.90)))

    estimated_short = (
        int(round(float(np.median(short_values))))
        if short_values
        else config.short_max_codons
    )
    estimated_long = (
        int(round(float(np.median(long_values))))
        if long_values
        else config.long_min_codons
    )
    short_cutoff = config.short_max_codons
    long_cutoff = config.long_min_codons
    if "short_max_codons" not in config.hidden_overrides and short_values:
        short_cutoff = int(
            np.clip(
                estimated_short,
                SHORT_MODEL_FLOOR,
                SHORT_MODEL_CEILING,
            )
        )
    if "long_min_codons" not in config.hidden_overrides and long_values:
        long_cutoff = int(
            np.clip(
                estimated_long,
                LONG_MODEL_FLOOR,
                LONG_MODEL_CEILING,
            )
        )
    if long_cutoff <= short_cutoff:
        long_cutoff = max(LONG_MODEL_FLOOR, short_cutoff + 1)

    logger.write(
        "Length models: estimated short<=%d and long>=%d codons; "
        "applied short<=%d, medium=%d-%d, long>=%d codons."
        % (
            estimated_short,
            estimated_long,
            short_cutoff,
            short_cutoff + 1,
            long_cutoff - 1,
            long_cutoff,
        )
    )
    return replace(
        config,
        short_max_codons=short_cutoff,
        long_min_codons=long_cutoff,
        estimated_short_max_codons=estimated_short,
        estimated_long_min_codons=estimated_long,
    )



def _longest_covered_run(
    positions: np.ndarray,
    coding_nt_length: int,
) -> int:
    """Return the longest consecutive covered-codon run."""
    valid = positions[(positions >= 0) & (positions < coding_nt_length)]
    if valid.size == 0:
        return 0
    codons = np.unique(valid // 3)
    if codons.size == 0:
        return 0
    breaks = np.flatnonzero(np.diff(codons) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks + 1, [codons.size]))
    return int(np.max(ends - starts))


def _select_expression_controls(
    records: list[dict[str, float | int]],
    config: EngineConfig,
) -> list[dict[str, float | int]]:
    """Select annotated_ORFs using the lower expression quantile."""
    if len(records) < config.positive_min_controls:
        return []
    ordered = sorted(
        records,
        key=lambda item: float(item["rpf_per_codon"]),
    )
    expression = np.asarray(
        [float(item["rpf_per_codon"]) for item in ordered],
        dtype=np.float64,
    )
    cutoff = float(np.quantile(expression, config.positive_quantile))
    selected = [
        item
        for item in ordered
        if float(item["rpf_per_codon"]) >= cutoff
    ]
    if len(selected) < config.positive_min_controls:
        selected = ordered[-config.positive_min_controls:]
    if len(selected) > config.positive_max_controls:
        indices = np.linspace(
            0,
            len(selected) - 1,
            config.positive_max_controls,
            dtype=np.int64,
        )
        selected = [selected[int(index)] for index in indices]
    return selected


def _derive_canonical_model(
    config: EngineConfig,
    records: Sequence[Mapping[str, float | int]],
    logger: _StageLogger,
) -> EngineConfig:
    """Derive non-stricter hidden parameters from annotated_ORF controls.

    Canonical controls may relax the balanced hidden defaults when library
    quality is poor, but must never make the smORF whole-ORF model stricter.
    """
    if not records:
        return config

    quantile = config.positive_quantile
    longest_runs = np.asarray(
        [float(item["longest_run"]) for item in records],
        dtype=np.float64,
    )
    signal_spans = np.asarray(
        [float(item["signal_span"]) for item in records],
        dtype=np.float64,
    )
    top_fractions = np.asarray(
        [float(item["top_window_fraction"]) for item in records],
        dtype=np.float64,
    )

    default_window = int(ADVANCED_DEFAULTS["window_codons"])
    window_codons = config.window_codons
    if "window_codons" not in config.hidden_overrides:
        window_codons = int(
            np.clip(
                round(float(np.quantile(longest_runs, quantile))),
                12,
                default_window,
            )
        )

    window_step = config.window_step_codons
    if "window_step_codons" not in config.hidden_overrides:
        window_step = max(
            1,
            min(
                int(ADVANCED_DEFAULTS["window_step_codons"]),
                int(round(window_codons * 0.25)),
            ),
        )

    window_gap = config.min_window_gap_codons
    if "min_window_gap_codons" not in config.hidden_overrides:
        window_gap = max(
            window_step,
            min(
                int(ADVANCED_DEFAULTS["min_window_gap_codons"]),
                int(round(window_codons * 0.75)),
            ),
        )

    supported_windows = config.min_supported_windows
    if "min_supported_windows" not in config.hidden_overrides:
        # Keep the validated balanced requirement. Annotated mORFs are too
        # long to provide a transferable absolute window-count threshold.
        supported_windows = int(
            ADVANCED_DEFAULTS["min_supported_windows"]
        )

    signal_span = config.min_signal_span
    if "min_signal_span" not in config.hidden_overrides:
        signal_span = float(
            np.clip(
                min(
                    float(ADVANCED_DEFAULTS["min_signal_span"]),
                    float(np.quantile(signal_spans, quantile)),
                ),
                0.20,
                float(ADVANCED_DEFAULTS["min_signal_span"]),
            )
        )

    localized_span = config.localized_span_max
    if "localized_span_max" not in config.hidden_overrides:
        localized_span = float(
            np.clip(
                min(
                    float(ADVANCED_DEFAULTS["localized_span_max"]),
                    float(
                        np.quantile(
                            signal_spans,
                            max(0.05, quantile / 2.0),
                        )
                    ),
                ),
                0.08,
                float(ADVANCED_DEFAULTS["localized_span_max"]),
            )
        )

    top_fraction = config.localized_top_window_fraction
    if "localized_top_window_fraction" not in config.hidden_overrides:
        # A larger threshold makes localized-only rejection more conservative.
        top_fraction = float(
            np.clip(
                max(
                    float(
                        ADVANCED_DEFAULTS[
                            "localized_top_window_fraction"
                        ]
                    ),
                    float(np.quantile(top_fractions, 1.0 - quantile)),
                ),
                float(
                    ADVANCED_DEFAULTS[
                        "localized_top_window_fraction"
                    ]
                ),
                0.90,
            )
        )

    boundary = config.boundary_codons
    if "boundary_codons" not in config.hidden_overrides:
        boundary = max(
            3,
            min(
                int(ADVANCED_DEFAULTS["boundary_codons"]),
                int(round(window_codons * 0.25)),
            ),
        )

    resolved = replace(
        config,
        window_codons=window_codons,
        window_step_codons=window_step,
        min_supported_windows=supported_windows,
        min_window_gap_codons=window_gap,
        min_signal_span=signal_span,
        localized_span_max=min(
            localized_span,
            signal_span * 0.90,
        ),
        localized_top_window_fraction=top_fraction,
        boundary_codons=boundary,
    )
    logger.write(
        "Canonical model: window=%d, step=%d, supported_windows=%d, gap=%d, "
        "signal_span=%.3f, localized_span=%.3f, top_window_fraction=%.3f, "
        "boundary=%d."
        % (
            resolved.window_codons,
            resolved.window_step_codons,
            resolved.min_supported_windows,
            resolved.min_window_gap_codons,
            resolved.min_signal_span,
            resolved.localized_span_max,
            resolved.localized_top_window_fraction,
            resolved.boundary_codons,
        )
    )
    return resolved


def _control_families(
    database_path: str | Path,
    category: str,
    maximum: int,
) -> list[tuple[str, FamilyGeometry]]:
    """Select deterministic annotated-ORF controls."""
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    query = """
        SELECT
            f.family_id, f.gene_id, f.chrom, f.strand,
            f.category AS family_category, f.family_type, f.family_size,
            f.structural_primary, r.orf_id, r.is_primary,
            g.transcript_id, g.category, g.start_codon, g.stop_codon,
            g.completeness, g.nt_length, g.aa_length,
            g.exon_starts, g.exon_ends
        FROM families f
        JOIN representatives r ON r.family_id=f.family_id
        JOIN geometry g ON g.orf_id=r.orf_id
        WHERE f.category=?
        ORDER BY f.family_id, g.nt_length DESC
        LIMIT ?
    """
    output: list[tuple[str, FamilyGeometry]] = []
    try:
        current: str | None = None
        rows: list[sqlite3.Row] = []
        for row in connection.execute(query, (category, maximum * 4)):
            family_id = row["family_id"]
            if current is None:
                current = family_id
            elif family_id != current:
                geometry = _build_family_geometry(rows)
                if isinstance(geometry, FamilyGeometry):
                    output.append((geometry.chrom, geometry))
                    if len(output) >= maximum:
                        break
                rows = []
                current = family_id
            rows.append(row)
        if rows and len(output) < maximum:
            geometry = _build_family_geometry(rows)
            if isinstance(geometry, FamilyGeometry):
                output.append((geometry.chrom, geometry))
    finally:
        connection.close()
    return output


def calibrate_samples(
    database_path: Path,
    work_directory: Path,
    tracks: Sequence[DensityTrack],
    config: EngineConfig,
    logger: _StageLogger,
) -> tuple[list[Thresholds], EngineConfig]:
    """Resolve manual thresholds or derive canonical thresholds from annotated_ORFs."""
    samples = sorted({track.sample for track in tracks})
    output_path = work_directory / "calibration.json"

    if config.evidence_mode == "manual":
        thresholds = [
            _manual_thresholds(sample, config)
            for sample in samples
        ]
        output_path.write_text(
            json.dumps(
                {
                    "mode": "manual",
                    "thresholds": [asdict(value) for value in thresholds],
                    "model": asdict(config),
                },
                indent=2,
                default=list,
            ),
            encoding="utf-8",
        )
        logger.write(
            "Manual evidence mode: use the eight user-defined thresholds."
        )
        return thresholds, config

    controls = _control_families(
        database_path,
        "annotated_ORF",
        max(
            config.positive_min_controls,
            config.positive_max_controls * 2,
        ),
    )
    if len(controls) < config.positive_min_controls:
        raise ValueError(
            "Canonical mode requires annotated_ORF controls, but only "
            f"{len(controls):,} valid controls were found. Re-run "
            "smorf_cluster with annotated_ORF retained in --keep-categories."
        )

    tracks_by_sample: dict[str, dict[str, DensityTrack]] = defaultdict(dict)
    for track in tracks:
        tracks_by_sample[track.sample][track.strand] = track
    controls_by_chrom: dict[str, list[FamilyGeometry]] = defaultdict(list)
    for chrom, geometry in controls:
        controls_by_chrom[chrom].append(geometry)

    provisional = _manual_thresholds("canonical_provisional", config)
    records_by_sample: dict[str, list[dict[str, float | int]]] = {}
    all_selected: list[dict[str, float | int]] = []

    for sample_number, sample in enumerate(samples, start=1):
        records: list[dict[str, float | int]] = []
        mapping = tracks_by_sample[sample]
        for chrom, families in controls_by_chrom.items():
            plus_track = mapping.get("+") or mapping.get(".")
            minus_track = mapping.get("-") or mapping.get(".")
            density_by_strand = {
                "+": (
                    _load_cached_density(
                        work_directory,
                        plus_track,
                        chrom,
                    )
                    if plus_track
                    else None
                ),
                "-": (
                    _load_cached_density(
                        work_directory,
                        minus_track,
                        chrom,
                    )
                    if minus_track
                    else None
                ),
            }
            for strand in ("+", "-"):
                subset = [
                    family
                    for family in families
                    if family.strand == strand
                ]
                if not subset:
                    continue
                signals = _extract_sparse_signals(
                    subset,
                    density_by_strand[strand],
                )
                for index, family in enumerate(subset):
                    signal = signals.get(index)
                    if signal is None:
                        continue
                    features = segment_features(
                        signal[0],
                        signal[1],
                        0,
                        family.scaffold_coding_nt_length,
                        provisional,
                        config,
                    )
                    if features.rpf_sum <= 0:
                        continue
                    codon_count = max(
                        1,
                        family.scaffold_coding_nt_length // 3,
                    )
                    records.append(
                        {
                            "rpf_sum": features.rpf_sum,
                            "rpf_per_codon": features.rpf_per_codon,
                            "covered_codon": features.covered_codon,
                            "coverage_ratio": features.coverage_ratio,
                            "frame0_ratio": features.frame0_ratio,
                            "signal_span": features.signal_span,
                            "top_window_fraction": (
                                features.top_window_fraction
                            ),
                            "codon_count": codon_count,
                            "longest_run": _longest_covered_run(
                                signal[0],
                                family.scaffold_coding_nt_length,
                            ),
                        }
                    )

        selected = _select_expression_controls(records, config)
        if len(selected) < config.positive_min_controls:
            raise ValueError(
                "Canonical mode requires at least "
                f"{config.positive_min_controls} expressed annotated_ORFs "
                f"per sample. Sample {sample} provided only "
                f"{len(selected)} after expression-quantile selection. "
                "Retain annotated_ORF in smorf_cluster or use "
                "--evidence-mode manual."
            )
        records_by_sample[sample] = selected
        all_selected.extend(selected)
        logger.write(
            "Canonical controls sample=%s: raw=%d, selected=%d."
            % (sample, len(records), len(selected))
        )
        progress_print(
            f"canonical controls: {sample_number}/{len(samples)}"
        )

    config = _derive_canonical_model(
        config,
        all_selected,
        logger,
    )
    thresholds: list[Thresholds] = []
    q = config.positive_quantile
    for sample in samples:
        selected = records_by_sample[sample]

        def array(name: str) -> np.ndarray:
            return np.asarray(
                [float(item[name]) for item in selected],
                dtype=np.float64,
            )

        rpf_sum = array("rpf_sum")
        rpf_per_codon = array("rpf_per_codon")
        covered = array("covered_codon")
        coverage = array("coverage_ratio")
        frame0 = array("frame0_ratio")
        expected_window_rpf = rpf_per_codon * config.window_codons
        expected_window_covered = coverage * config.window_codons

        moderate = float(
            np.clip(np.quantile(frame0, q), 0.40, 0.85)
        )
        strong = float(
            np.clip(
                np.median(frame0),
                moderate + 0.05,
                0.95,
            )
        )
        # Canonical controls calibrate the lower data-quality boundary.
        # Thresholds may become more permissive than the balanced defaults,
        # but never stricter, because annotated mORFs are generally much more
        # abundant and longer than smORFs.
        moderate = max(
            0.40,
            min(
                float(MANUAL_THRESHOLD_DEFAULTS["moderate_periodicity"]),
                moderate,
            ),
        )
        strong = max(
            moderate + 0.05,
            min(
                float(MANUAL_THRESHOLD_DEFAULTS["strong_periodicity"]),
                strong,
            ),
        )
        threshold = Thresholds(
            sample=sample,
            source="canonical:annotated_ORF",
            control_count=len(selected),
            min_rpf_sum=max(
                1.0,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_rpf_sum"]),
                    float(np.quantile(rpf_sum, q)),
                ),
            ),
            min_rpf_per_codon=max(
                0.01,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_rpf_per_codon"]),
                    float(np.quantile(rpf_per_codon, q)),
                ),
            ),
            min_covered_codon=max(
                1,
                min(
                    int(MANUAL_THRESHOLD_DEFAULTS["min_covered_codon"]),
                    int(math.floor(float(np.quantile(covered, q)))),
                ),
            ),
            min_coverage_ratio=max(
                0.01,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_coverage_ratio"]),
                    float(np.quantile(coverage, q)),
                ),
            ),
            moderate_periodicity=moderate,
            strong_periodicity=min(0.95, strong),
            min_window_rpf=max(
                1.0,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_window_rpf"]),
                    float(np.quantile(expected_window_rpf, q)),
                ),
            ),
            min_window_covered=max(
                1,
                min(
                    int(MANUAL_THRESHOLD_DEFAULTS["min_window_covered"]),
                    int(
                        math.floor(
                            float(
                                np.quantile(
                                    expected_window_covered,
                                    q,
                                )
                            )
                        )
                    ),
                ),
            ),
        )
        thresholds.append(threshold)
        logger.write(
            "Canonical thresholds sample=%s: controls=%d, RPF=%.3f, "
            "RPF/codon=%.4f, covered=%d, coverage=%.3f, "
            "moderate=%.3f, strong=%.3f, window_RPF=%.3f, "
            "window_covered=%d."
            % (
                sample,
                threshold.control_count,
                threshold.min_rpf_sum,
                threshold.min_rpf_per_codon,
                threshold.min_covered_codon,
                threshold.min_coverage_ratio,
                threshold.moderate_periodicity,
                threshold.strong_periodicity,
                threshold.min_window_rpf,
                threshold.min_window_covered,
            )
        )

    output_path.write_text(
        json.dumps(
            {
                "mode": "canonical",
                "thresholds": [asdict(value) for value in thresholds],
                "model": asdict(config),
            },
            indent=2,
            default=list,
        ),
        encoding="utf-8",
    )
    return thresholds, config



def _new_accumulator(
    family: FamilyGeometry,
    config: EngineConfig,
) -> _FamilyAccumulator:
    """Create a bounded accumulator for one family."""
    candidate_count = len(family.representatives)
    return _FamilyAccumulator(
        any_signal=False,
        any_evidence=False,
        support_count=0,
        high_count=0,
        group_mask=0,
        group_support_counts={},
        group_high_counts={},
        supporting_samples=[],
        best_level="NoEvidence",
        best_score=-1.0,
        best_sample="",
        best_features=None,
        candidate_full_support_counts=[0] * candidate_count,
        candidate_extension_support_counts=[0] * candidate_count,
        candidate_support_samples=[[] for _ in range(candidate_count)],
        candidate_bin_support_counts=[
            (
                [0] * min(
                    family.representatives[index].coding_nt_length // 3,
                    config.extent_bins,
                )
                if family.representatives[index].coding_nt_length // 3
                >= config.long_min_codons
                else []
            )
            for index in range(candidate_count)
        ],
        candidate_bin_support_samples=[
            [] for _ in range(candidate_count)
        ],
        candidate_cohort_support_samples=[
            [] for _ in range(candidate_count)
        ],
        candidate_bin_pooled_stats=[
            (
                [
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    for _ in range(
                        min(
                            family.representatives[index].coding_nt_length // 3,
                            config.extent_bins,
                        )
                    )
                ]
                if family.representatives[index].coding_nt_length // 3
                >= config.long_min_codons
                else []
            )
            for index in range(candidate_count)
        ],
        candidate_boundary_pooled_stats=[
            (
                [
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                ]
                if family.representatives[index].coding_nt_length // 3
                >= config.long_min_codons
                else []
            )
            for index in range(candidate_count)
        ],
        candidate_best_full_features=[None] * candidate_count,
        candidate_best_extension_features=[None] * candidate_count,
        candidate_pooled_frame_density=[
            [0.0, 0.0, 0.0] for _ in range(candidate_count)
        ],
        candidate_leading_support_counts=[0] * candidate_count,
        candidate_leading_negative_counts=[0] * candidate_count,
        candidate_noncanonical_extension_support_counts=[0]
        * candidate_count,
        candidate_best_noncanonical_density_ratios=[0.0]
        * candidate_count,
    )



def _output_columns() -> tuple[str, ...]:
    """Return focused family output columns."""
    return (
        "family_id", "gene_id", "chrom", "strand", "category",
        "family_type", "family_size", "structural_primary",
        "common_body_orf", "provisional_primary", "quant_primary",
        "evidence_primary",
        "translation_unit_id", "translation_unit_index",
        "translation_unit_count", "translation_unit_status",
        "evidence_status", "family_translation_status",
        "translated_extent_status", "extent_supported_sample_count",
        "extent_supporting_samples", "extent_supported_bins",
        "extension_rpf_sum", "extension_frame0_ratio",
        "family_translation_evidence", "reliability_reason",
        "start_site_status", "start_site_reason", "start_interval_orf_ids",
        "selection_policy", "selection_reason", "canonical_anchor_orf",
        "longest_candidate_orf", "prior_override_status",
        "leading_support_sample_count", "leading_negative_sample_count",
        "noncanonical_extension_support_sample_count",
        "noncanonical_extension_density_ratio",
        "pooled_extent_rpf_sum", "pooled_extent_frame0_ratio",
        "pooled_start_rpf_sum", "pooled_end_rpf_sum",
        "pooled_extent_sample_count",
        "nested_competition_status", "dominant_parent_orf",
        "competition_overlap_fraction", "competition_frame_relation",
        "competition_reason", "frame_margin",
        "reliable_group",
        "supported_sample_count", "max_group_sample_support",
        "high_confidence_sample_count", "supporting_group_count",
        "supporting_groups", "supporting_samples", "best_sample",
        "best_sample_evidence", "best_rpf_sum", "best_rpf_per_codon",
        "best_covered_codon", "best_coverage_ratio",
        "best_frame0_ratio", "best_supported_windows",
        "best_distributed_windows", "best_signal_span",
        "best_top_window_fraction", "best_localized_only",
        "representative_count", "representative_orf_ids",
        "representative_start_codons", "invalid_family_reason",
        "complete_candidate_count", "complete_candidate_orf_ids",
        "candidate_audit_status", "candidate_failure_reasons",
        "pooled_frame0_density", "pooled_frame1_density",
        "pooled_frame2_density", "pooled_dominant_frame",
        "suggested_psite_shift_nt", "phase_audit_status",
    )



def _candidate_extent_supported(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    config: EngineConfig,
) -> bool:
    """Return whether one candidate has direct or cohort-wide extent support."""
    representative = family.representatives[candidate_index]
    codon_count = representative.coding_nt_length // 3
    cohort_status = _candidate_cohort_extent_status(
        family,
        accumulator,
        candidate_index,
        config,
    )
    if cohort_status in {"Supported", "Compatible"}:
        return True
    if (
        accumulator.candidate_full_support_counts[candidate_index]
        >= config.reliable_sample
    ):
        return True
    if (
        codon_count < config.long_min_codons
        and _frame_margin(
            accumulator.candidate_pooled_frame_density[candidate_index]
        )
        < config.min_frame_margin
    ):
        return False
    extension_codons = max(
        0,
        (family.common_body_start - representative.start_offset) // 3,
    )
    extension_supported = (
        extension_codons <= config.start_resolution_codons
        or accumulator.candidate_extension_support_counts[candidate_index]
        >= config.reliable_sample
    )
    bin_counts = accumulator.candidate_bin_support_counts[candidate_index]
    supported_bins = _candidate_supported_bin_indices(
        accumulator,
        candidate_index,
        config,
    )
    minimum_bins = min(config.min_extent_bins, len(bin_counts))
    upstream_limit = max(1, math.ceil(len(bin_counts) * 0.40))
    upstream_supported = any(index < upstream_limit for index in supported_bins)
    return (
        accumulator.support_count >= config.reliable_sample
        and extension_supported
        and len(supported_bins) >= minimum_bins
        and upstream_supported
    )


def _prior_candidate_index(family: FamilyGeometry) -> int:
    """Return the structural prior used for candidate-level auditing."""
    canonical_index = _canonical_anchor_index(family)
    if canonical_index is not None:
        return canonical_index
    return min(
        range(len(family.representatives)),
        key=lambda index: (
            family.representatives[index].start_offset,
            family.representatives[index].orf_id,
        ),
    )


def _candidate_failure_reasons(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    config: EngineConfig,
) -> tuple[str, ...]:
    """Return explicit unmet candidate-first evidence requirements."""
    representative = family.representatives[candidate_index]
    reasons: list[str] = []
    if not _candidate_is_complete(representative):
        reasons.append("incomplete_start_to_stop_candidate")
    pooled_frames = accumulator.candidate_pooled_frame_density[candidate_index]
    if sum(pooled_frames) <= 0:
        reasons.append("no_candidate_p_site_signal")
    sample_count = _candidate_cohort_sample_count(accumulator, candidate_index)
    if sample_count < config.reliable_sample:
        reasons.append(
            f"candidate_sample_replication_{sample_count}_below_"
            f"{config.reliable_sample}"
        )
    bin_counts = accumulator.candidate_bin_support_counts[candidate_index]
    supported_bins = _candidate_supported_bin_indices(
        accumulator,
        candidate_index,
        config,
    )
    if bin_counts:
        minimum_bins = min(config.min_extent_bins, len(bin_counts))
        if len(supported_bins) < minimum_bins:
            reasons.append(
                f"supported_extent_bins_{len(supported_bins)}_below_"
                f"{minimum_bins}"
            )
        start_supported, end_supported = _candidate_boundary_support(
            accumulator,
            candidate_index,
            config,
        )
        if not supported_bins or supported_bins[0] != 0 or not start_supported:
            reasons.append("missing_candidate_N_terminal_phase_support")
        downstream_start = max(1, math.floor(len(bin_counts) * 0.60))
        if not any(index >= downstream_start for index in supported_bins):
            reasons.append("missing_candidate_downstream_extent_support")
        if (
            supported_bins
            and supported_bins[-1] == len(bin_counts) - 1
            and not end_supported
        ):
            reasons.append("missing_candidate_terminal_phase_support")
        segmented_frames = _candidate_segmented_frame_density(
            accumulator,
            candidate_index,
            supported_bins,
        )
        if (
            supported_bins
            and _frame_margin(segmented_frames) < config.min_frame_margin
        ):
            reasons.append("candidate_segmented_frame_margin_below_threshold")
    elif (
        sum(pooled_frames) > 0
        and _frame_margin(pooled_frames) < config.min_frame_margin
    ):
        reasons.append("candidate_frame_margin_below_threshold")
    if not reasons and not _candidate_extent_supported(
        family,
        accumulator,
        candidate_index,
        config,
    ):
        reasons.append("candidate_extent_not_reliably_resolved")
    return tuple(reasons)


def _phase_audit(
    frame_density: Sequence[float],
    min_margin: float,
) -> tuple[int, int, str]:
    """Return dominant frame, suggested coordinate shift, and audit status."""
    frames = [float(value) for value in frame_density]
    total = sum(frames)
    if total <= 0:
        return -1, 0, "NoSignal"
    order = sorted(range(3), key=lambda index: frames[index], reverse=True)
    best = order[0]
    advantage = (frames[best] - frames[order[1]]) / total
    if advantage < min_margin:
        return best, 0, "AmbiguousFrame"
    if best == 0:
        return 0, 0, "Frame0Consistent"
    suggested_shift = -1 if best == 1 else 1
    return best, suggested_shift, "AlternativeFrameDominant"


def _frame_margin(frame_density: Sequence[float]) -> float:
    """Calculate the candidate frame-0 advantage over the alternative frames."""
    total = float(sum(frame_density))
    if total <= 0:
        return 0.0
    ratios = [float(value) / total for value in frame_density]
    return ratios[0] - max(ratios[1], ratios[2])


def _resolve_candidate_extent(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    config: EngineConfig,
) -> ExtentResolution:
    """Select an ORF using canonical/length priors plus replicated evidence."""
    selection_policy = (
        "canonical_ATG_then_longest_unless_replicated_RPF_override"
    )
    supported_indices = [
        index
        for index in range(len(family.representatives))
        if _candidate_extent_supported(family, accumulator, index, config)
    ]
    canonical_index = _canonical_anchor_index(family)
    longest_index = min(
        range(len(family.representatives)),
        key=lambda index: (
            family.representatives[index].start_offset,
            family.representatives[index].orf_id,
        ),
    )
    prior_index = canonical_index if canonical_index is not None else longest_index

    def contradicted(index: int) -> bool:
        return (
            accumulator.candidate_leading_negative_counts[index]
            >= config.reliable_sample
            and accumulator.candidate_leading_support_counts[index]
            < config.reliable_sample
        )

    selected_index = prior_index
    if canonical_index is not None:
        strong_upstream_noncanonical = [
            index
            for index, representative in enumerate(family.representatives)
            if representative.start_offset
            < family.representatives[canonical_index].start_offset
            and representative.start_codon.upper() != "ATG"
            and accumulator.candidate_noncanonical_extension_support_counts[
                index
            ]
            >= config.reliable_sample
            and (
                accumulator.candidate_leading_support_counts[index]
                >= config.reliable_sample
                or (
                    family.representatives[canonical_index].start_offset
                    - representative.start_offset
                )
                // 3
                < config.min_exclusion_codons
            )
            and not contradicted(index)
        ]
        if strong_upstream_noncanonical:
            selected_index = min(
                strong_upstream_noncanonical,
                key=lambda index: (
                    family.representatives[index].start_offset,
                    family.representatives[index].orf_id,
                ),
            )
            selection_reason = (
                "upstream_noncanonical_start_has_strong_exclusive_"
                "replicated_translation"
            )
            override_status = "OverriddenByStrongNoncanonicalExtension"
        else:
            selection_reason = "longest_canonical_ATG_retained"
            override_status = "CanonicalATGRetained"
    else:
        selection_reason = "longest_candidate_retained"
        override_status = "LongestCandidateRetained"

    if selected_index == prior_index and contradicted(prior_index):
        downstream_supported = [
            index
            for index in supported_indices
            if family.representatives[index].start_offset
            > family.representatives[prior_index].start_offset
            and not contradicted(index)
        ]
        if downstream_supported:
            selected_index = min(
                downstream_supported,
                key=lambda index: (
                    family.representatives[index].start_codon.upper() != "ATG",
                    family.representatives[index].start_offset,
                    family.representatives[index].orf_id,
                ),
            )
            selection_reason = (
                "longer_prior_rejected_by_replicated_N_terminal_silence_"
                "and_internal_extent_support"
            )
            override_status = "ShortenedByReplicatedNterminalSilence"
        else:
            diagnostic_frames = (
                accumulator.candidate_pooled_frame_density[prior_index]
            )
            diagnostic_total = float(sum(diagnostic_frames))
            diagnostic_boundary = (
                accumulator.candidate_boundary_pooled_stats[prior_index]
            )
            diagnostic_extension = (
                accumulator.candidate_best_extension_features[prior_index]
            )
            diagnostic_bins = tuple(
                index + 1
                for index in _candidate_supported_bin_indices(
                    accumulator,
                    prior_index,
                    config,
                )
            )
            diagnostic_samples = tuple(
                dict.fromkeys(
                    accumulator.candidate_support_samples[prior_index]
                    + accumulator.candidate_bin_support_samples[prior_index]
                    + accumulator.candidate_cohort_support_samples[prior_index]
                )
            )
            return ExtentResolution(
                quant_primary=family.representatives[prior_index].orf_id,
                evidence_primary="",
                translated_extent_status="Partial",
                extent_supported_sample_count=_candidate_cohort_sample_count(
                    accumulator,
                    prior_index,
                ),
                extent_supporting_samples=diagnostic_samples,
                extent_supported_bins=diagnostic_bins,
                extension_rpf_sum=(
                    0.0
                    if diagnostic_extension is None
                    else diagnostic_extension.rpf_sum
                ),
                extension_frame0_ratio=(
                    0.0
                    if diagnostic_extension is None
                    else diagnostic_extension.frame0_ratio
                ),
                start_site_status="Unresolved",
                start_site_reason=(
                    "longer_prior_contradicted_without_supported_internal_"
                    "replacement"
                ),
                start_interval_orf_ids=tuple(
                    family.representatives[index].orf_id
                    for index in supported_indices
                ),
                frame_margin=_frame_margin(diagnostic_frames),
                selection_policy=selection_policy,
                selection_reason=(
                    "replicated_N_terminal_silence_without_replacement"
                ),
                canonical_anchor_orf=(
                    ""
                    if canonical_index is None
                    else family.representatives[canonical_index].orf_id
                ),
                longest_candidate_orf=(
                    family.representatives[longest_index].orf_id
                ),
                prior_override_status="NoDefensibleExtent",
                leading_support_sample_count=(
                    accumulator.candidate_leading_support_counts[prior_index]
                ),
                leading_negative_sample_count=(
                    accumulator.candidate_leading_negative_counts[prior_index]
                ),
                noncanonical_extension_support_sample_count=0,
                noncanonical_extension_density_ratio=0.0,
                pooled_extent_rpf_sum=diagnostic_total,
                pooled_extent_frame0_ratio=(
                    diagnostic_frames[0] / diagnostic_total
                    if diagnostic_total > 0
                    else 0.0
                ),
                pooled_start_rpf_sum=(
                    float(diagnostic_boundary[0][0])
                    if diagnostic_boundary
                    else 0.0
                ),
                pooled_end_rpf_sum=(
                    float(diagnostic_boundary[1][0])
                    if len(diagnostic_boundary) > 1
                    else 0.0
                ),
                pooled_extent_sample_count=_candidate_cohort_sample_count(
                    accumulator,
                    prior_index,
                ),
            )

    selected = family.representatives[selected_index]
    selected_is_directly_supported = selected_index in supported_indices
    cohort_extent_status = _candidate_cohort_extent_status(
        family,
        accumulator,
        selected_index,
        config,
    )
    interval_ids = tuple(
        dict.fromkeys(
            [selected.orf_id]
            + [
                family.representatives[index].orf_id
                for index in supported_indices
            ]
        )
    )
    nearby = [
        item.orf_id
        for item in family.representatives
        if abs(item.start_offset - selected.start_offset)
        <= config.start_resolution_codons * 3
    ]
    if len(family.representatives) == 1 and selected_is_directly_supported:
        start_status = "UniqueCandidate"
        start_reason = "singleton_family"
        evidence_primary = selected.orf_id
    elif not selected_is_directly_supported:
        start_status = "Compatible"
        start_reason = selection_reason + "_without_negative_counterevidence"
        evidence_primary = ""
    elif len(nearby) > 1:
        start_status = "Compatible"
        start_reason = selection_reason + "_nearby_starts_not_separable"
        evidence_primary = ""
    elif len(supported_indices) == 1:
        start_status = "Exact"
        start_reason = selection_reason + "_unique_supported_boundary"
        evidence_primary = selected.orf_id
    else:
        start_status = "Interval"
        start_reason = selection_reason + "_multiple_supported_boundaries"
        evidence_primary = ""

    full_features = accumulator.candidate_best_full_features[selected_index]
    extension_features = accumulator.candidate_best_extension_features[selected_index]
    pooled_frames = accumulator.candidate_pooled_frame_density[selected_index]
    supported_bins = tuple(
        index + 1
        for index in _candidate_supported_bin_indices(
            accumulator,
            selected_index,
            config,
        )
    )
    support_samples = tuple(
        dict.fromkeys(
            accumulator.candidate_support_samples[selected_index]
            + accumulator.candidate_bin_support_samples[selected_index]
            + accumulator.candidate_cohort_support_samples[selected_index]
            + accumulator.supporting_samples
        )
    )
    extension_total = (
        0.0 if extension_features is None else extension_features.rpf_sum
    )
    extension_frame0 = (
        0.0 if extension_features is None else extension_features.frame0_ratio
    )
    pooled_total = float(sum(pooled_frames))
    boundary = accumulator.candidate_boundary_pooled_stats[selected_index]
    return ExtentResolution(
        quant_primary=selected.orf_id,
        evidence_primary=evidence_primary,
        translated_extent_status=(
            "Supported"
            if selected_is_directly_supported
            else (
                cohort_extent_status
                if cohort_extent_status
                else "Compatible"
            )
        ),
        extent_supported_sample_count=max(
            accumulator.candidate_full_support_counts[selected_index],
            accumulator.candidate_extension_support_counts[selected_index],
            accumulator.candidate_leading_support_counts[selected_index],
            len(
                set(
                    accumulator.candidate_bin_support_samples[
                        selected_index
                    ]
                )
            ),
            _candidate_cohort_sample_count(
                accumulator,
                selected_index,
            ),
        ),
        extent_supporting_samples=support_samples,
        extent_supported_bins=supported_bins,
        extension_rpf_sum=extension_total,
        extension_frame0_ratio=extension_frame0,
        start_site_status=start_status,
        start_site_reason=start_reason,
        start_interval_orf_ids=interval_ids,
        frame_margin=_frame_margin(pooled_frames),
        selection_policy=selection_policy,
        selection_reason=selection_reason,
        canonical_anchor_orf=(
            ""
            if canonical_index is None
            else family.representatives[canonical_index].orf_id
        ),
        longest_candidate_orf=family.representatives[longest_index].orf_id,
        prior_override_status=override_status,
        leading_support_sample_count=(
            accumulator.candidate_leading_support_counts[selected_index]
        ),
        leading_negative_sample_count=(
            accumulator.candidate_leading_negative_counts[selected_index]
        ),
        noncanonical_extension_support_sample_count=(
            accumulator.candidate_noncanonical_extension_support_counts[
                selected_index
            ]
        ),
        noncanonical_extension_density_ratio=(
            accumulator.candidate_best_noncanonical_density_ratios[
                selected_index
            ]
        ),
        pooled_extent_rpf_sum=pooled_total,
        pooled_extent_frame0_ratio=(
            pooled_frames[0] / pooled_total if pooled_total > 0 else 0.0
        ),
        pooled_start_rpf_sum=(
            float(boundary[0][0]) if boundary else 0.0
        ),
        pooled_end_rpf_sum=(
            float(boundary[1][0]) if len(boundary) > 1 else 0.0
        ),
        pooled_extent_sample_count=_candidate_cohort_sample_count(
            accumulator,
            selected_index,
        ),
    )




def _family_output_row(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    group_names: Sequence[str],
    config: EngineConfig,
) -> dict[str, Any]:
    """Convert one family accumulator into a sample-replicated final row.

    ``--group`` is annotation metadata. Reliability is determined from
    independent sample replication across the whole project, not restricted to
    one group.
    """
    supporting_group_indices = sorted(
        accumulator.group_support_counts
    )
    supporting_groups = ",".join(
        str(group_names[index])
        for index in supporting_group_indices
    )
    max_group_support = max(
        accumulator.group_support_counts.values(),
        default=0,
    )

    resolution = _resolve_candidate_extent(
        family,
        accumulator,
        config,
    )
    provisional_primary = resolution.quant_primary
    selected_index = next(
        (
            index
            for index, representative in enumerate(family.representatives)
            if representative.orf_id == provisional_primary
        ),
        _prior_candidate_index(family),
    )
    candidate_extent_indices = [
        index
        for index in range(len(family.representatives))
        if _candidate_extent_supported(
            family,
            accumulator,
            index,
            config,
        )
    ]
    cohort_extent_candidates = [
        (index, status)
        for index in range(len(family.representatives))
        if (
            status := _candidate_cohort_extent_status(
                family,
                accumulator,
                index,
                config,
            )
        )
        in {"Supported", "Compatible"}
    ]
    common_body_reliable = (
        accumulator.support_count >= config.reliable_sample
    )
    candidate_extent_reliable = bool(candidate_extent_indices)
    cohort_extent_reliable = bool(cohort_extent_candidates)

    if (
        common_body_reliable
        or candidate_extent_reliable
        or cohort_extent_reliable
    ):
        status = "Reliable"
        final_evidence = (
            "HighConfidence"
            if common_body_reliable
            and accumulator.high_count >= config.reliable_sample
            else "MediumConfidence"
        )
        if common_body_reliable:
            reason = (
                f"common_body_supported_by_{accumulator.support_count}_"
                f"independent_samples_across_"
                f"{len(supporting_group_indices)}_groups"
            )
        else:
            extent_candidates = (
                cohort_extent_candidates
                if cohort_extent_candidates
                else [
                    (index, "Supported")
                    for index in candidate_extent_indices
                ]
            )
            best_index, best_status = max(
                extent_candidates,
                key=lambda item: (
                    item[1] == "Supported",
                    _candidate_cohort_sample_count(
                        accumulator,
                        item[0],
                    ),
                    -family.representatives[item[0]].start_offset,
                ),
            )
            reason = (
                "distributed_long_ORF_extent_"
                f"{best_status.lower()}_by_"
                f"{_candidate_cohort_sample_count(accumulator, best_index)}_"
                "independent_samples"
            )
        if len(supporting_group_indices) == 1:
            reliable_group = str(
                group_names[supporting_group_indices[0]]
            )
        elif len(supporting_group_indices) > 1:
            reliable_group = "cross_group"
        else:
            reliable_group = "cohort_distributed"
    elif not accumulator.any_signal:
        status = "NoEvidence"
        final_evidence = "NoEvidence"
        reliable_group = ""
        reason = "no_common_body_p_site_signal_in_any_sample"
        resolution = replace(
            resolution,
            quant_primary="",
            evidence_primary="",
            translated_extent_status="NotEvaluated",
            start_site_status="NotEvaluated",
            start_site_reason="family_translation_not_detected",
            selection_reason="family_translation_not_detected",
            prior_override_status="NotEvaluated",
        )
    else:
        status = "Uncertain"
        final_evidence = accumulator.best_level
        reliable_group = ""
        reason = (
            "insufficient_independent_sample_replication"
            if accumulator.any_evidence
            else "signal_present_but_no_sample_passed_translation_thresholds"
        )
        failure_reasons = _candidate_failure_reasons(
            family,
            accumulator,
            selected_index,
            config,
        )
        resolution = replace(
            resolution,
            quant_primary="",
            evidence_primary="",
            translated_extent_status="Unsupported",
            start_site_status="Unresolved",
            start_site_reason=(
                ";".join(failure_reasons)
                if failure_reasons
                else "family_translation_not_reliably_replicated"
            ),
            selection_reason="candidate_first_resolution_unresolved",
        )

    features = accumulator.best_features
    complete_candidates = tuple(
        item.orf_id
        for item in family.representatives
        if _candidate_is_complete(item)
    )
    failure_reasons = _candidate_failure_reasons(
        family,
        accumulator,
        selected_index,
        config,
    )
    selected_frames = accumulator.candidate_pooled_frame_density[
        selected_index
    ]
    dominant_frame, suggested_shift, phase_status = _phase_audit(
        selected_frames,
        config.min_frame_margin,
    )
    if resolution.quant_primary:
        candidate_audit_status = "InputCompleteCandidateSelected"
    elif complete_candidates:
        candidate_audit_status = "InputCompleteCandidateUnresolved"
    else:
        candidate_audit_status = "NoCompleteInputCandidate"
    return {
        "family_id": family.family_id,
        "gene_id": family.gene_id,
        "chrom": family.chrom,
        "strand": family.strand,
        "category": family.category,
        "family_type": family.family_type,
        "family_size": family.family_size,
        "structural_primary": family.structural_primary,
        "common_body_orf": family.common_body_orf,
        "provisional_primary": provisional_primary,
        "quant_primary": resolution.quant_primary,
        "evidence_primary": resolution.evidence_primary,
        "translation_unit_id": "",
        "translation_unit_index": 0,
        "translation_unit_count": 0,
        "translation_unit_status": "NotAssigned",
        "evidence_status": status,
        "family_translation_status": status,
        "translated_extent_status": resolution.translated_extent_status,
        "extent_supported_sample_count": (
            resolution.extent_supported_sample_count
        ),
        "extent_supporting_samples": ",".join(
            resolution.extent_supporting_samples
        ),
        "extent_supported_bins": ",".join(
            str(value) for value in resolution.extent_supported_bins
        ),
        "extension_rpf_sum": resolution.extension_rpf_sum,
        "extension_frame0_ratio": resolution.extension_frame0_ratio,
        "family_translation_evidence": final_evidence,
        "reliability_reason": reason,
        "start_site_status": resolution.start_site_status,
        "start_site_reason": resolution.start_site_reason,
        "start_interval_orf_ids": ",".join(
            resolution.start_interval_orf_ids
        ),
        "selection_policy": resolution.selection_policy,
        "selection_reason": resolution.selection_reason,
        "canonical_anchor_orf": resolution.canonical_anchor_orf,
        "longest_candidate_orf": resolution.longest_candidate_orf,
        "prior_override_status": resolution.prior_override_status,
        "leading_support_sample_count": (
            resolution.leading_support_sample_count
        ),
        "leading_negative_sample_count": (
            resolution.leading_negative_sample_count
        ),
        "noncanonical_extension_support_sample_count": (
            resolution.noncanonical_extension_support_sample_count
        ),
        "noncanonical_extension_density_ratio": (
            resolution.noncanonical_extension_density_ratio
        ),
        "pooled_extent_rpf_sum": resolution.pooled_extent_rpf_sum,
        "pooled_extent_frame0_ratio": resolution.pooled_extent_frame0_ratio,
        "pooled_start_rpf_sum": resolution.pooled_start_rpf_sum,
        "pooled_end_rpf_sum": resolution.pooled_end_rpf_sum,
        "pooled_extent_sample_count": resolution.pooled_extent_sample_count,
        "nested_competition_status": (
            "Pending"
            if status == "Reliable"
            and resolution.translated_extent_status
            in {"Supported", "Compatible"}
            else "NotEvaluated"
        ),
        "dominant_parent_orf": "",
        "frame_margin": resolution.frame_margin,
        "reliable_group": reliable_group,
        "supported_sample_count": accumulator.support_count,
        "max_group_sample_support": max_group_support,
        "high_confidence_sample_count": accumulator.high_count,
        "supporting_group_count": len(supporting_group_indices),
        "supporting_groups": supporting_groups,
        "supporting_samples": ",".join(
            accumulator.supporting_samples
        ),
        "best_sample": accumulator.best_sample,
        "best_sample_evidence": accumulator.best_level,
        "best_rpf_sum": (
            0.0 if features is None else features.rpf_sum
        ),
        "best_rpf_per_codon": (
            0.0 if features is None else features.rpf_per_codon
        ),
        "best_covered_codon": (
            0 if features is None else features.covered_codon
        ),
        "best_coverage_ratio": (
            0.0 if features is None else features.coverage_ratio
        ),
        "best_frame0_ratio": (
            0.0 if features is None else features.frame0_ratio
        ),
        "best_supported_windows": (
            0 if features is None else features.supported_windows
        ),
        "best_distributed_windows": (
            0 if features is None else features.distributed_windows
        ),
        "best_signal_span": (
            0.0 if features is None else features.signal_span
        ),
        "best_top_window_fraction": (
            0.0
            if features is None
            else features.top_window_fraction
        ),
        "best_localized_only": (
            False if features is None else features.localized_only
        ),
        "representative_count": len(family.representatives),
        "representative_orf_ids": ",".join(
            item.orf_id for item in family.representatives
        ),
        "representative_start_codons": ",".join(
            item.start_codon for item in family.representatives
        ),
        "invalid_family_reason": "",
        "complete_candidate_count": len(complete_candidates),
        "complete_candidate_orf_ids": ",".join(complete_candidates),
        "candidate_audit_status": candidate_audit_status,
        "candidate_failure_reasons": ";".join(failure_reasons),
        "pooled_frame0_density": selected_frames[0],
        "pooled_frame1_density": selected_frames[1],
        "pooled_frame2_density": selected_frames[2],
        "pooled_dominant_frame": dominant_frame,
        "suggested_psite_shift_nt": suggested_shift,
        "phase_audit_status": phase_status,
    }




def _invalid_output_row(family: InvalidFamily) -> dict[str, Any]:
    """Convert an invalid family into an explicit uncertain row."""
    return {
        "family_id": family.family_id,
        "gene_id": family.gene_id,
        "chrom": family.chrom,
        "strand": family.strand,
        "category": family.category,
        "family_type": family.family_type,
        "family_size": family.family_size,
        "structural_primary": family.structural_primary,
        "common_body_orf": "",
        "provisional_primary": "",
        "quant_primary": "",
        "evidence_primary": "",
        "translation_unit_id": "",
        "translation_unit_index": 0,
        "translation_unit_count": 0,
        "translation_unit_status": "NotAssigned",
        "evidence_status": "Uncertain",
        "family_translation_status": "Uncertain",
        "translated_extent_status": "NotEvaluated",
        "extent_supported_sample_count": 0,
        "extent_supporting_samples": "",
        "extent_supported_bins": "",
        "extension_rpf_sum": 0.0,
        "extension_frame0_ratio": 0.0,
        "family_translation_evidence": "InvalidFamilyGeometry",
        "reliability_reason": "family_geometry_validation_failed",
        "start_site_status": "Unresolved",
        "start_site_reason": "invalid_family_geometry",
        "start_interval_orf_ids": "",
        "selection_policy": (
            "canonical_ATG_then_longest_unless_replicated_RPF_override"
        ),
        "selection_reason": "invalid_family_geometry",
        "canonical_anchor_orf": "",
        "longest_candidate_orf": "",
        "prior_override_status": "NotEvaluated",
        "leading_support_sample_count": 0,
        "leading_negative_sample_count": 0,
        "noncanonical_extension_support_sample_count": 0,
        "noncanonical_extension_density_ratio": 0.0,
        "pooled_extent_rpf_sum": 0.0,
        "pooled_extent_frame0_ratio": 0.0,
        "pooled_start_rpf_sum": 0.0,
        "pooled_end_rpf_sum": 0.0,
        "pooled_extent_sample_count": 0,
        "nested_competition_status": "NotEvaluated",
        "dominant_parent_orf": "",
        "frame_margin": 0.0,
        "reliable_group": "",
        "supported_sample_count": 0,
        "max_group_sample_support": 0,
        "high_confidence_sample_count": 0,
        "supporting_group_count": 0,
        "supporting_groups": "",
        "supporting_samples": "",
        "best_sample": "",
        "best_sample_evidence": "InvalidFamilyGeometry",
        "best_rpf_sum": 0.0,
        "best_rpf_per_codon": 0.0,
        "best_covered_codon": 0,
        "best_coverage_ratio": 0.0,
        "best_frame0_ratio": 0.0,
        "best_supported_windows": 0,
        "best_distributed_windows": 0,
        "best_signal_span": 0.0,
        "best_top_window_fraction": 0.0,
        "best_localized_only": False,
        "representative_count": 0,
        "representative_orf_ids": "",
        "representative_start_codons": "",
        "invalid_family_reason": family.reason,
        "complete_candidate_count": 0,
        "complete_candidate_orf_ids": "",
        "candidate_audit_status": "InvalidFamilyGeometry",
        "candidate_failure_reasons": "invalid_family_geometry",
        "pooled_frame0_density": 0.0,
        "pooled_frame1_density": 0.0,
        "pooled_frame2_density": 0.0,
        "pooled_dominant_frame": -1,
        "suggested_psite_shift_nt": 0,
        "phase_audit_status": "NotEvaluated",
    }


def _format_value(value: Any) -> str:
    if isinstance(value, float):
        return format(value, ".10g") if math.isfinite(value) else "NA"
    if isinstance(value, bool):
        return "True" if value else "False"
    return str(value)


def _is_annotated_category(value: Any) -> bool:
    """Return whether a category is reserved for canonical calibration."""
    return str(value).strip().lower() in ANNOTATED_CATEGORIES


def _write_rows(handle: TextIO, rows: Sequence[Mapping[str, Any]]) -> None:
    columns = _output_columns()
    handle.write(
        "".join(
            "\t".join(_format_value(row.get(column, "")) for column in columns)
            + "\n"
            for row in rows
        )
    )


def _process_chromosome(task: ChromosomeTask) -> ChromosomeResult:
    """Process one chromosome in bounded family batches."""
    work_directory = Path(task.work_directory)
    part_directory = work_directory / "parts"
    part_directory.mkdir(parents=True, exist_ok=True)
    safe_chrom = _safe_name(task.chromosome)
    part_path = part_directory / f"{safe_chrom}.part.tsv.gz"
    count_path = part_directory / f"{safe_chrom}.counts.json"
    done_path = part_directory / f"{safe_chrom}.done.json"
    if done_path.is_file() and part_path.is_file() and count_path.is_file():
        values = json.loads(count_path.read_text(encoding="utf-8"))
        return ChromosomeResult(
            chromosome=task.chromosome,
            part_path=str(part_path),
            count_path=str(count_path),
            **values,
        )

    thresholds_by_sample = {value.sample: value for value in task.thresholds}
    tracks_by_sample: dict[str, dict[str, DensityTrack]] = defaultdict(dict)
    sample_groups: dict[str, str] = {}
    for track in task.tracks:
        tracks_by_sample[track.sample][track.strand] = track
        sample_groups[track.sample] = track.group
    samples = sorted(tracks_by_sample)
    group_index = {group: index for index, group in enumerate(task.group_names)}
    if len(group_index) > 63:
        raise ValueError("At most 63 biological groups are supported.")

    # Open each sample/strand chromosome cache once per worker task. The arrays
    # are memory-mapped, so keeping the handles does not load all density data
    # into resident memory and avoids reopening thousands of files per batch.
    density_by_sample_strand: dict[tuple[str, str], ChromDensity | None] = {}
    for sample in samples:
        for strand in ("+", "-"):
            track = _track_for_family(tracks_by_sample, sample, strand)
            density_by_sample_strand[(sample, strand)] = (
                _load_cached_density(work_directory, track, task.chromosome)
                if track is not None
                else None
            )

    temporary = part_path.with_name(part_path.name + ".tmp")
    counts = {
        "total_families": 0,
        "reliable_families": 0,
        "uncertain_families": 0,
        "no_evidence_families": 0,
        "reliable_smorfs": 0,
        "invalid_families": 0,
    }
    with gzip.open(
        temporary,
        "wt",
        encoding="utf-8",
        compresslevel=1,
        newline="",
    ) as output_handle:
        for batch_number, batch in enumerate(
            _iter_chromosome_family_batches(
                task.database_path,
                task.chromosome,
            ),
            start=1,
        ):
            valid_families = [item for item in batch if isinstance(item, FamilyGeometry)]
            invalid_families = [item for item in batch if isinstance(item, InvalidFamily)]
            accumulators = [
                _new_accumulator(family, task.config)
                for family in valid_families
            ]
            canonical_indices = [
                _canonical_anchor_index(family)
                for family in valid_families
            ]
            diagnostic_index_sets = []
            for family, canonical_index in zip(
                valid_families,
                canonical_indices,
            ):
                if canonical_index is None:
                    diagnostic_index_sets.append({0})
                    continue
                canonical_start = family.representatives[
                    canonical_index
                ].start_offset
                diagnostic_index_sets.append(
                    {canonical_index}
                    | {
                        index
                        for index, representative in enumerate(
                            family.representatives
                        )
                        if representative.start_offset < canonical_start
                        and representative.start_codon.upper() != "ATG"
                    }
                )

            for sample in samples:
                threshold = thresholds_by_sample[sample]
                for strand in ("+", "-"):
                    indices = [
                        index for index, family in enumerate(valid_families)
                        if family.strand == strand
                    ]
                    if not indices:
                        continue
                    subset = [valid_families[index] for index in indices]
                    density = density_by_sample_strand[(sample, strand)]
                    signals = _extract_sparse_signals(subset, density)
                    for subset_index, signal in signals.items():
                        family_index = indices[subset_index]
                        family = valid_families[family_index]
                        accumulator = accumulators[family_index]
                        profile = ScaffoldProfile(
                            signal[0],
                            signal[1],
                            family.scaffold_coding_nt_length,
                        )
                        features = profile.region_features(
                            family.common_body_start,
                            family.common_body_end,
                            threshold,
                            task.config,
                        )
                        if features.rpf_sum > 0:
                            accumulator.any_signal = True
                        codon_count = (
                            family.common_body_end - family.common_body_start
                        ) // 3
                        level, _ = classify_sample(
                            features,
                            codon_count,
                            threshold,
                            task.config,
                        )
                        if LEVEL_RANK[level] >= LEVEL_RANK["LowConfidence"]:
                            accumulator.any_evidence = True
                        if (
                            LEVEL_RANK[level] > LEVEL_RANK[accumulator.best_level]
                            or (
                                LEVEL_RANK[level] == LEVEL_RANK[accumulator.best_level]
                                and features.score > accumulator.best_score
                            )
                        ):
                            accumulator.best_level = level
                            accumulator.best_score = features.score
                            accumulator.best_sample = sample
                            accumulator.best_features = features
                        group = sample_groups[sample]
                        group_id = group_index[group]
                        if LEVEL_RANK[level] >= LEVEL_RANK["MediumConfidence"]:
                            accumulator.support_count += 1
                            accumulator.supporting_samples.append(sample)
                            accumulator.group_mask |= 1 << group_id
                            accumulator.group_support_counts[group_id] = (
                                accumulator.group_support_counts.get(group_id, 0)
                                + 1
                            )
                            if level == "HighConfidence":
                                accumulator.high_count += 1
                                accumulator.group_high_counts[group_id] = (
                                    accumulator.group_high_counts.get(group_id, 0)
                                    + 1
                                )

                        canonical_index = canonical_indices[family_index]
                        diagnostic_indices = diagnostic_index_sets[family_index]
                        canonical_features = features
                        canonical_start: int | None = None
                        if canonical_index is not None:
                            canonical_start = family.representatives[
                                canonical_index
                            ].start_offset
                            if canonical_start != family.common_body_start:
                                canonical_features = profile.region_features(
                                    canonical_start,
                                    family.common_body_end,
                                    threshold,
                                    task.config,
                                )

                        for candidate_index, representative in enumerate(
                            family.representatives
                        ):
                            candidate_start = representative.start_offset
                            candidate_end = family.common_body_end
                            extension_end = family.common_body_start
                            full = (
                                features
                                if candidate_start == family.common_body_start
                                else profile.region_features(
                                    candidate_start, candidate_end,
                                    threshold, task.config,
                                )
                            )
                            if full.rpf_sum > 0:
                                accumulator.any_signal = True
                            extension = profile.region_features(
                                candidate_start, extension_end,
                                threshold, task.config,
                            )
                            full_codons = max(
                                0, (candidate_end - candidate_start) // 3
                            )
                            extension_codons = max(
                                0, (extension_end - candidate_start) // 3
                            )
                            full_supported, extension_supported = (
                                classify_candidate_extent(
                                    full,
                                    extension,
                                    full_codons,
                                    extension_codons,
                                    threshold,
                                    task.config,
                                )
                            )
                            if full_supported:
                                accumulator.candidate_full_support_counts[
                                    candidate_index
                                ] += 1
                            if extension_supported:
                                accumulator.candidate_extension_support_counts[
                                    candidate_index
                                ] += 1
                            if full_supported or (
                                extension_supported and extension.rpf_sum > 0
                            ):
                                accumulator.candidate_support_samples[
                                    candidate_index
                                ].append(sample)

                            if candidate_index in diagnostic_indices:
                                leading_end = _leading_region_end(
                                    family,
                                    candidate_index,
                                    canonical_index,
                                    task.config,
                                )
                                leading_codons = max(
                                    0,
                                    (leading_end - candidate_start) // 3,
                                )
                                leading = profile.region_features(
                                    candidate_start,
                                    leading_end,
                                    threshold,
                                    task.config,
                                )
                                if _leading_region_is_supported(
                                    leading,
                                    leading_codons,
                                    threshold,
                                    task.config,
                                ):
                                    accumulator.candidate_leading_support_counts[
                                        candidate_index
                                    ] += 1
                                if (
                                    LEVEL_RANK[level]
                                    >= LEVEL_RANK["MediumConfidence"]
                                    and _leading_region_is_negative(
                                        leading,
                                        features,
                                        leading_codons,
                                        threshold,
                                        task.config,
                                    )
                                ):
                                    accumulator.candidate_leading_negative_counts[
                                        candidate_index
                                    ] += 1

                            if (
                                canonical_index is not None
                                and canonical_start is not None
                                and candidate_start < canonical_start
                                and representative.start_codon.upper() != "ATG"
                            ):
                                noncanonical_extension = (
                                    profile.region_features(
                                        candidate_start,
                                        canonical_start,
                                        threshold,
                                        task.config,
                                    )
                                )
                                noncanonical_codons = (
                                    canonical_start - candidate_start
                                ) // 3
                                density_ratio = (
                                    noncanonical_extension.rpf_per_codon
                                    / canonical_features.rpf_per_codon
                                    if canonical_features.rpf_per_codon > 0
                                    else 0.0
                                )
                                accumulator.candidate_best_noncanonical_density_ratios[
                                    candidate_index
                                ] = max(
                                    accumulator.candidate_best_noncanonical_density_ratios[
                                        candidate_index
                                    ],
                                    density_ratio,
                                )
                                if _noncanonical_extension_is_strong(
                                    noncanonical_extension,
                                    canonical_features,
                                    noncanonical_codons,
                                    threshold,
                                    task.config,
                                ):
                                    accumulator.candidate_noncanonical_extension_support_counts[
                                        candidate_index
                                    ] += 1

                            previous_full = (
                                accumulator.candidate_best_full_features[
                                    candidate_index
                                ]
                            )
                            if (
                                previous_full is None
                                or full.score > previous_full.score
                            ):
                                accumulator.candidate_best_full_features[
                                    candidate_index
                                ] = full
                            previous_extension = (
                                accumulator.candidate_best_extension_features[
                                    candidate_index
                                ]
                            )
                            if (
                                previous_extension is None
                                or extension.score > previous_extension.score
                            ):
                                accumulator.candidate_best_extension_features[
                                    candidate_index
                                ] = extension
                            pooled = accumulator.candidate_pooled_frame_density[
                                candidate_index
                            ]
                            pooled[0] += full.frame0_density
                            pooled[1] += full.frame1_density
                            pooled[2] += full.frame2_density
                            if (
                                full.rpf_sum > 0
                                and full.frame0_density >= 1.0
                                and _frame_margin(
                                    (
                                        full.frame0_density,
                                        full.frame1_density,
                                        full.frame2_density,
                                    )
                                )
                                >= task.config.min_frame_margin
                            ):
                                accumulator.candidate_cohort_support_samples[
                                    candidate_index
                                ].append(sample)

                            bin_ranges = (
                                _candidate_bin_ranges(
                                    candidate_start,
                                    candidate_end,
                                    task.config.extent_bins,
                                )
                                if full_codons >= task.config.long_min_codons
                                else ()
                            )
                            sample_has_supported_bin = False
                            for bin_index, (bin_start, bin_end) in enumerate(
                                bin_ranges
                            ):
                                bin_features = profile.region_features(
                                    bin_start, bin_end,
                                    threshold, task.config,
                                )
                                _accumulate_pooled_segment(
                                    accumulator.candidate_bin_pooled_stats[
                                        candidate_index
                                    ][bin_index],
                                    bin_features,
                                )
                                if _bin_is_supported(
                                    bin_features,
                                    threshold,
                                ):
                                    accumulator.candidate_bin_support_counts[
                                        candidate_index
                                    ][bin_index] += 1
                                    sample_has_supported_bin = True
                            if sample_has_supported_bin:
                                accumulator.candidate_bin_support_samples[
                                    candidate_index
                                ].append(sample)
                            if bin_ranges:
                                boundary_codons = min(
                                    task.config.window_codons,
                                    full_codons,
                                )
                                start_boundary = profile.region_features(
                                    candidate_start,
                                    candidate_start + boundary_codons * 3,
                                    threshold,
                                    task.config,
                                )
                                end_boundary = profile.region_features(
                                    candidate_end - boundary_codons * 3,
                                    candidate_end,
                                    threshold,
                                    task.config,
                                )
                                _accumulate_pooled_segment(
                                    accumulator.candidate_boundary_pooled_stats[
                                        candidate_index
                                    ][0],
                                    start_boundary,
                                )
                                _accumulate_pooled_segment(
                                    accumulator.candidate_boundary_pooled_stats[
                                        candidate_index
                                    ][1],
                                    end_boundary,
                                )

            rows: list[Mapping[str, Any]] = []
            for family, accumulator in zip(valid_families, accumulators):
                row = _family_output_row(
                    family, accumulator, task.group_names, task.config
                )
                rows.append(row)
                counts["total_families"] += 1
                if row["evidence_status"] == "Reliable":
                    counts["reliable_families"] += 1
                    if (
                        row["translated_extent_status"]
                        in {"Supported", "Compatible"}
                        and row["quant_primary"]
                        and not _is_annotated_category(row["category"])
                    ):
                        counts["reliable_smorfs"] += 1
                elif row["evidence_status"] == "NoEvidence":
                    counts["no_evidence_families"] += 1
                else:
                    counts["uncertain_families"] += 1
            for invalid in invalid_families:
                rows.append(_invalid_output_row(invalid))
                counts["total_families"] += 1
                counts["uncertain_families"] += 1
                counts["invalid_families"] += 1
            _write_rows(output_handle, rows)
            progress_print(
                f"{task.chromosome}: family batches {batch_number:,}, "
                f"families {counts['total_families']:,}"
            )

    os.replace(temporary, part_path)
    count_path.write_text(json.dumps(counts, indent=2), encoding="utf-8")
    done_path.write_text(
        json.dumps({"chromosome": task.chromosome}),
        encoding="utf-8",
    )
    return ChromosomeResult(
        chromosome=task.chromosome,
        part_path=str(part_path),
        count_path=str(count_path),
        **counts,
    )


def _reliable_matches(
    connection: sqlite3.Connection,
    orf_ids: Sequence[str],
) -> dict[str, int]:
    """Return reliable-index match states for one genePred input batch."""
    unique_ids = tuple(dict.fromkeys(orf_ids))
    matches: dict[str, int] = {}
    query_size = 800
    for start in range(0, len(unique_ids), query_size):
        chunk = unique_ids[start : start + query_size]
        placeholders = ",".join("?" for _ in chunk)
        query = (
            "SELECT orf_id, matched FROM reliable_export "
            f"WHERE orf_id IN ({placeholders})"
        )
        for row in connection.execute(query, chunk):
            matches[str(row["orf_id"])] = int(row["matched"])
    return matches


def _export_reliable_genepred(
    source_path: str | Path,
    temporary_path: Path,
    connection: sqlite3.Connection,
) -> int:
    """Filter scanner genePred records using the validated Reliable ID index."""
    expected = int(
        connection.execute(
            "SELECT COUNT(*) FROM reliable_export"
        ).fetchone()[0]
    )
    if expected == 0:
        temporary_path.write_text("", encoding="utf-8")
        return 0

    duplicate_count = 0
    duplicate_examples: list[str] = []
    batch: list[tuple[str, str, int]] = []
    batch_size = 25_000

    def flush_batch(output_handle: TextIO) -> None:
        nonlocal duplicate_count
        if not batch:
            return
        matches = _reliable_matches(
            connection,
            [orf_id for orf_id, _line, _number in batch],
        )
        newly_matched: set[str] = set()
        for orf_id, raw_line, line_number in batch:
            matched = matches.get(orf_id)
            if matched is None:
                continue
            if matched or orf_id in newly_matched:
                duplicate_count += 1
                if len(duplicate_examples) < 10:
                    duplicate_examples.append(
                        f"{orf_id}@line{line_number}"
                    )
                continue
            if len(raw_line.rstrip("\r\n").split("\t")) < 10:
                raise EvidenceEngineError(
                    "Reliable genePred record has fewer than 10 columns at "
                    f"line {line_number}: {orf_id}"
                )
            output_handle.write(raw_line)
            if not raw_line.endswith(("\n", "\r")):
                output_handle.write("\n")
            newly_matched.add(orf_id)
        if newly_matched:
            connection.executemany(
                "UPDATE reliable_export SET matched=1 WHERE orf_id=?",
                ((orf_id,) for orf_id in newly_matched),
            )
        batch.clear()

    with _smart_open(source_path, "rt") as source, temporary_path.open(
        "w",
        encoding="utf-8",
        buffering=8 * 1024 * 1024,
    ) as target:
        for line_number, raw_line in enumerate(source, start=1):
            if raw_line.startswith("#"):
                flush_batch(target)
                target.write(raw_line)
                continue
            text = raw_line.rstrip("\r\n")
            if not text:
                continue
            orf_id = text.split("\t", 1)[0]
            batch.append((orf_id, raw_line, line_number))
            if len(batch) >= batch_size:
                flush_batch(target)
        flush_batch(target)

    if duplicate_count:
        raise EvidenceEngineError(
            "Scanner genePred contains duplicate records for "
            f"{duplicate_count:,} Reliable ORF occurrence(s). Examples: "
            + ", ".join(duplicate_examples)
        )

    missing_count = int(
        connection.execute(
            "SELECT COUNT(*) FROM reliable_export WHERE matched=0"
        ).fetchone()[0]
    )
    if missing_count:
        examples = [
            str(row[0])
            for row in connection.execute(
                "SELECT orf_id FROM reliable_export "
                "WHERE matched=0 ORDER BY orf_id LIMIT 10"
            )
        ]
        raise EvidenceEngineError(
            f"Scanner genePred is missing {missing_count:,} Reliable ORF(s). "
            "Examples: " + ", ".join(examples)
        )
    return expected


def _geometry_blocks(row: Mapping[str, Any]) -> tuple[tuple[int, int], ...]:
    """Return sorted genomic blocks from one geometry row."""
    starts = _split_ints(row["exon_starts"])
    ends = _split_ints(row["exon_ends"])
    return tuple(sorted(zip(starts, ends), key=lambda item: item[0]))


def _blocks_contain(
    parent: Sequence[tuple[int, int]],
    child: Sequence[tuple[int, int]],
) -> bool:
    """Return whether every child coding block is covered by parent blocks."""
    parent_index = 0
    for child_start, child_end in child:
        while (
            parent_index < len(parent)
            and parent[parent_index][1] <= child_start
        ):
            parent_index += 1
        if parent_index >= len(parent):
            return False
        parent_start, parent_end = parent[parent_index]
        if parent_start > child_start or parent_end < child_end:
            return False
    return True


def _block_length(blocks: Sequence[tuple[int, int]]) -> int:
    """Return the total spliced length of genomic blocks."""
    return sum(max(0, end - start) for start, end in blocks)


def _block_overlap_length(
    first: Sequence[tuple[int, int]],
    second: Sequence[tuple[int, int]],
) -> int:
    """Return the genomic overlap length between two sorted block sets."""
    first_index = 0
    second_index = 0
    overlap = 0
    while first_index < len(first) and second_index < len(second):
        first_start, first_end = first[first_index]
        second_start, second_end = second[second_index]
        overlap += max(
            0,
            min(first_end, second_end) - max(first_start, second_start),
        )
        if first_end <= second_end:
            first_index += 1
        else:
            second_index += 1
    return overlap


def _shorter_overlap_fraction(
    first: Sequence[tuple[int, int]],
    second: Sequence[tuple[int, int]],
) -> float:
    """Return overlap as a fraction of the shorter spliced ORF."""
    denominator = min(_block_length(first), _block_length(second))
    if denominator <= 0:
        return 0.0
    return _block_overlap_length(first, second) / denominator


def _first_overlap_position(
    first: Sequence[tuple[int, int]],
    second: Sequence[tuple[int, int]],
    strand: str,
) -> int | None:
    """Return the first shared genomic base in transcript orientation."""
    overlaps: list[tuple[int, int]] = []
    first_index = 0
    second_index = 0
    while first_index < len(first) and second_index < len(second):
        start = max(first[first_index][0], second[second_index][0])
        end = min(first[first_index][1], second[second_index][1])
        if end > start:
            overlaps.append((start, end))
        if first[first_index][1] <= second[second_index][1]:
            first_index += 1
        else:
            second_index += 1
    if not overlaps:
        return None
    return overlaps[0][0] if strand == "+" else overlaps[-1][1] - 1


def _same_overlap_frame(
    first_blocks: Sequence[tuple[int, int]],
    second_blocks: Sequence[tuple[int, int]],
    strand: str,
) -> bool:
    """Return whether two overlapping ORFs use the same coding frame."""
    position = _first_overlap_position(first_blocks, second_blocks, strand)
    if position is None:
        return False
    first_offset = _transcript_offset(first_blocks, strand, position)
    second_offset = _transcript_offset(second_blocks, strand, position)
    return (
        first_offset is not None
        and second_offset is not None
        and first_offset % 3 == second_offset % 3
    )


def _transcript_offset(
    blocks: Sequence[tuple[int, int]],
    strand: str,
    genomic_position: int,
) -> int | None:
    """Map one genomic base to a spliced offset from the parent ORF start."""
    offset = 0
    for start, end in _oriented_blocks(
        [item[0] for item in blocks],
        [item[1] for item in blocks],
        strand,
    ):
        if start <= genomic_position < end:
            return (
                offset + genomic_position - start
                if strand == "+"
                else offset + end - 1 - genomic_position
            )
        offset += end - start
    return None


def _same_coding_frame(
    parent_blocks: Sequence[tuple[int, int]],
    child_blocks: Sequence[tuple[int, int]],
    strand: str,
) -> bool:
    """Return whether the child start is frame-compatible with the parent."""
    oriented_child = _oriented_blocks(
        [item[0] for item in child_blocks],
        [item[1] for item in child_blocks],
        strand,
    )
    if not oriented_child:
        return False
    child_start = (
        oriented_child[0][0]
        if strand == "+"
        else oriented_child[0][1] - 1
    )
    offset = _transcript_offset(parent_blocks, strand, child_start)
    return offset is not None and offset % 3 == 0


def _assign_translation_units(
    parsed: Sequence[
        tuple[sqlite3.Row, tuple[tuple[int, int], ...]]
    ],
) -> dict[str, tuple[str, int, int, str]]:
    """Assign overlap-connected ORFs to independent transcript units."""
    if not parsed:
        return {}
    parents = list(range(len(parsed)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first: int, second: int) -> None:
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parents[second_root] = first_root

    interval_order = sorted(
        range(len(parsed)),
        key=lambda index: (
            min(start for start, _ in parsed[index][1]),
            max(end for _, end in parsed[index][1]),
        ),
    )
    active: list[int] = []
    for index in interval_order:
        blocks = parsed[index][1]
        left = min(start for start, _ in blocks)
        active = [
            other
            for other in active
            if max(end for _, end in parsed[other][1]) > left
        ]
        for other in active:
            if _block_overlap_length(blocks, parsed[other][1]) > 0:
                union(index, other)
        active.append(index)

    components: dict[int, list[int]] = defaultdict(list)
    for index in range(len(parsed)):
        components[find(index)].append(index)
    first_row = parsed[0][0]
    strand = str(first_row["strand"])
    ordered_components = sorted(
        components.values(),
        key=lambda indices: (
            min(
                start
                for index in indices
                for start, _ in parsed[index][1]
            )
            if strand == "+"
            else -max(
                end
                for index in indices
                for _, end in parsed[index][1]
            )
        ),
    )
    unit_count = len(ordered_components)
    gene_id = str(first_row["gene_id"])
    transcript_id = str(first_row["transcript_id"])
    assignments: dict[str, tuple[str, int, int, str]] = {}
    for unit_index, indices in enumerate(ordered_components, start=1):
        unit_id = f"{gene_id}|{transcript_id}|TU{unit_index:03d}"
        unit_status = (
            "OverlappingCompetitionUnit"
            if len(indices) > 1
            else "IndependentNonOverlappingUnit"
        )
        for index in indices:
            assignments[str(parsed[index][0]["family_id"])] = (
                unit_id,
                unit_index,
                unit_count,
                unit_status,
            )
    return assignments


def _apply_nested_competition(
    raw_master_path: Path,
    final_master_path: Path,
    connection: sqlite3.Connection,
    config: EngineConfig,
) -> tuple[dict[str, int], dict[str, int]]:
    """Assign transcript units and suppress explainable nested ORFs."""
    connection.execute("DROP TABLE IF EXISTS temp.candidate_calls")
    connection.execute(
        """
        CREATE TEMP TABLE candidate_calls (
            family_id TEXT PRIMARY KEY,
            orf_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            strand TEXT NOT NULL,
            category TEXT NOT NULL,
            evidence_status TEXT NOT NULL,
            frame_margin REAL NOT NULL,
            phase_rpf REAL NOT NULL,
            extent_sample_count INTEGER NOT NULL,
            extent_status TEXT NOT NULL,
            start_site_status TEXT NOT NULL,
            prior_override_status TEXT NOT NULL,
            leading_support_count INTEGER NOT NULL,
            noncanonical_support_count INTEGER NOT NULL,
            competition_status TEXT NOT NULL DEFAULT 'Independent',
            dominant_parent_orf TEXT NOT NULL DEFAULT '',
            overlap_fraction REAL NOT NULL DEFAULT 0.0,
            frame_relation TEXT NOT NULL DEFAULT 'NotEvaluated',
            competition_reason TEXT NOT NULL DEFAULT ''
        ) WITHOUT ROWID
        """
    )
    batch: list[tuple[Any, ...]] = []
    with gzip.open(raw_master_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            provisional_primary = (
                row.get("provisional_primary")
                or row.get("quant_primary")
                or ""
            )
            if not provisional_primary:
                continue
            if (
                int(row.get("complete_candidate_count", "0") or 0) < 1
                and not row.get("quant_primary")
            ):
                continue
            batch.append(
                (
                    row["family_id"],
                    provisional_primary,
                    row["gene_id"],
                    row["chrom"],
                    row["strand"],
                    row["category"],
                    row.get("evidence_status", "Uncertain"),
                    float(row.get("frame_margin", "0") or 0.0),
                    float(
                        row.get("pooled_frame0_density", "0")
                        or row.get("best_rpf_sum", "0")
                        or 0.0
                    ),
                    int(
                        row.get("pooled_extent_sample_count", "0")
                        or row.get("extent_supported_sample_count", "0")
                        or 0
                    ),
                    row.get("translated_extent_status", ""),
                    row.get("start_site_status", ""),
                    row.get("prior_override_status", ""),
                    int(row.get("leading_support_sample_count", "0") or 0),
                    int(
                        row.get(
                            "noncanonical_extension_support_sample_count",
                            "0",
                        )
                        or 0
                    ),
                )
            )
            if len(batch) >= SQL_BATCH_SIZE:
                connection.executemany(
                    "INSERT INTO candidate_calls "
                    "(family_id,orf_id,gene_id,chrom,strand,category,"
                    "evidence_status,frame_margin,phase_rpf,"
                    "extent_sample_count,extent_status,start_site_status,"
                    "prior_override_status,leading_support_count,"
                    "noncanonical_support_count) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                batch.clear()
        if batch:
            connection.executemany(
                "INSERT INTO candidate_calls "
                "(family_id,orf_id,gene_id,chrom,strand,category,"
                "evidence_status,frame_margin,phase_rpf,"
                "extent_sample_count,extent_status,start_site_status,"
                "prior_override_status,leading_support_count,"
                "noncanonical_support_count) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
    connection.execute(
        "CREATE INDEX candidate_calls_group_idx ON "
        "candidate_calls(gene_id, chrom, strand)"
    )

    query = """
        SELECT c.*, g.transcript_id, g.start_codon, g.nt_length,
               g.exon_starts, g.exon_ends
        FROM candidate_calls c
        JOIN geometry g ON g.orf_id=c.orf_id
        ORDER BY c.gene_id, g.transcript_id, c.chrom, c.strand,
                 g.nt_length DESC, c.orf_id
    """
    updates: list[tuple[str, str, float, str, str, str]] = []
    group_rows: list[sqlite3.Row] = []
    group_key: tuple[str, str, str, str] | None = None
    unit_assignments: dict[str, tuple[str, int, int, str]] = {}

    def resolve_group(rows: Sequence[sqlite3.Row]) -> None:
        parsed: list[tuple[sqlite3.Row, tuple[tuple[int, int], ...]]] = [
            (row, _geometry_blocks(row)) for row in rows
        ]
        current_units = _assign_translation_units(parsed)
        unit_assignments.update(current_units)
        suppressed_parent_families: set[str] = set()
        for child_index, (child, child_blocks) in enumerate(parsed):
            competitors: list[
                tuple[
                    sqlite3.Row,
                    tuple[tuple[int, int], ...],
                    float,
                    bool,
                ]
            ] = []
            for parent, parent_blocks in parsed[:child_index]:
                if (
                    current_units[str(parent["family_id"])][0]
                    != current_units[str(child["family_id"])][0]
                ):
                    continue
                if str(parent["family_id"]) in suppressed_parent_families:
                    continue
                if str(parent["evidence_status"]) != "Reliable":
                    continue
                if int(parent["nt_length"]) <= int(child["nt_length"]):
                    continue
                contained = _blocks_contain(parent_blocks, child_blocks)
                overlap_fraction = _shorter_overlap_fraction(
                    parent_blocks,
                    child_blocks,
                )
                if (
                    contained
                    or overlap_fraction >= config.high_overlap_fraction
                ):
                    competitors.append(
                        (
                            parent,
                            parent_blocks,
                            overlap_fraction,
                            contained,
                        )
                    )
            if not competitors:
                continue
            parent, parent_blocks, overlap_fraction, contained = max(
                competitors,
                key=lambda value: (
                    value[2],
                    str(value[0]["start_codon"]).upper() == "ATG",
                    int(value[0]["nt_length"]),
                    str(value[0]["orf_id"]),
                ),
            )
            same_frame = _same_overlap_frame(
                parent_blocks,
                child_blocks,
                str(child["strand"]),
            )
            frame_relation = "SameFrame" if same_frame else "DifferentFrame"
            parent_noncanonical_supported = (
                str(parent["start_codon"]).upper() != "ATG"
                and str(parent["extent_status"]) == "Supported"
                and int(parent["leading_support_count"])
                >= config.reliable_sample
                and int(parent["noncanonical_support_count"])
                >= config.reliable_sample
                and str(parent["prior_override_status"])
                == "OverriddenByStrongNoncanonicalExtension"
            )
            canonical_child_over_parent = (
                same_frame
                and str(child["start_codon"]).upper() == "ATG"
                and str(parent["start_codon"]).upper() != "ATG"
                and not parent_noncanonical_supported
            )
            if canonical_child_over_parent:
                updates.append(
                    (
                        "ShadowedByCanonicalATG",
                        str(child["orf_id"]),
                        overlap_fraction,
                        frame_relation,
                        "canonical_ATG_preferred_without_replicated_"
                        "noncanonical_exclusive_translation",
                        str(parent["family_id"]),
                    )
                )
                suppressed_parent_families.add(str(parent["family_id"]))
                continue
            strong_overlap = (
                overlap_fraction >= config.high_overlap_fraction
            )
            required_margin = (
                config.overlap_min_frame_margin
                if strong_overlap
                else config.nested_min_frame_margin
            )
            required_rpf = (
                config.overlap_min_phase_rpf
                if strong_overlap
                else config.nested_min_phase_rpf
            )
            if (
                str(parent["start_codon"]).upper() == "ATG"
                and str(child["start_codon"]).upper() != "ATG"
            ):
                required_margin += 0.05
                required_rpf *= 1.25
            required_rpf = max(
                required_rpf,
                float(parent["phase_rpf"]) * 0.25,
            )
            independent_phase = (
                str(child["extent_status"]) == "Supported"
                and
                float(child["frame_margin"])
                >= required_margin
                and float(child["phase_rpf"])
                >= required_rpf
                and int(child["extent_sample_count"])
                >= config.reliable_sample
            )
            start_site_resolved = str(child["start_site_status"]) in {
                "Exact",
                "UniqueCandidate",
            }
            noncanonical_start_resolved = (
                str(child["start_codon"]).upper() == "ATG"
                or (
                    str(child["prior_override_status"])
                    == "OverriddenByStrongNoncanonicalExtension"
                    and int(child["leading_support_count"])
                    >= config.reliable_sample
                    and int(child["noncanonical_support_count"])
                    >= config.reliable_sample
                )
            )
            independent_start = (
                start_site_resolved and noncanonical_start_resolved
            )
            if not same_frame and independent_phase and independent_start:
                updates.append(
                    (
                        "IndependentOverlappingORF",
                        str(parent["orf_id"]),
                        overlap_fraction,
                        frame_relation,
                        "independent_frame_advantage_and_replicated_"
                        "translation_support",
                        str(child["family_id"]),
                    )
                )
            else:
                status = (
                    "ShadowedByLongORF"
                    if contained
                    else "ShadowedByHighOverlapORF"
                )
                reason = (
                    "same_frame_overlap_explained_by_longer_ORF"
                    if same_frame
                    else "insufficient_independent_start_frame_or_"
                    "abundance_evidence_in_high_overlap_region"
                )
                updates.append(
                    (
                        status,
                        str(parent["orf_id"]),
                        overlap_fraction,
                        frame_relation,
                        reason,
                        str(child["family_id"]),
                    )
                )

    for row in connection.execute(query):
        current_key = (
            str(row["gene_id"]),
            str(row["transcript_id"]),
            str(row["chrom"]),
            str(row["strand"]),
        )
        if group_key is None:
            group_key = current_key
        elif current_key != group_key:
            resolve_group(group_rows)
            group_rows = []
            group_key = current_key
        group_rows.append(row)
    if group_rows:
        resolve_group(group_rows)
    if updates:
        connection.executemany(
            "UPDATE candidate_calls SET competition_status=?, "
            "dominant_parent_orf=?, overlap_fraction=?, frame_relation=?, "
            "competition_reason=? WHERE family_id=?",
            updates,
        )

    competition = {
        str(row["family_id"]): (
            str(row["competition_status"]),
            str(row["dominant_parent_orf"]),
            float(row["overlap_fraction"]),
            str(row["frame_relation"]),
            str(row["competition_reason"]),
        )
        for row in connection.execute(
            "SELECT family_id, competition_status, dominant_parent_orf, "
            "overlap_fraction, frame_relation, competition_reason "
            "FROM candidate_calls"
        )
    }
    counts: dict[str, int] = defaultdict(int)
    unit_counts: dict[str, int] = defaultdict(int)
    for unit_id, _, _, unit_status in set(unit_assignments.values()):
        unit_counts[unit_status] += 1
    with gzip.open(raw_master_path, "rt", encoding="utf-8") as source, final_master_path.open(
        "w",
        encoding="utf-8",
        buffering=8 * 1024 * 1024,
        newline="",
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(
            target,
            fieldnames=_output_columns(),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in reader:
            status = competition.get(row["family_id"])
            if status is not None:
                row["nested_competition_status"] = status[0]
                row["dominant_parent_orf"] = status[1]
                row["competition_overlap_fraction"] = status[2]
                row["competition_frame_relation"] = status[3]
                row["competition_reason"] = status[4]
                counts[status[0]] += 1
            unit = unit_assignments.get(row["family_id"])
            if unit is not None:
                row["translation_unit_id"] = unit[0]
                row["translation_unit_index"] = unit[1]
                row["translation_unit_count"] = unit[2]
                row["translation_unit_status"] = unit[3]
            writer.writerow(
                {
                    column: _format_value(row.get(column, ""))
                    for column in _output_columns()
                }
            )
    return dict(counts), dict(unit_counts)


def _merge_outputs(
    output_prefix: Path,
    chromosome_results: Sequence[ChromosomeResult],
    database_path: Path,
    config: EngineConfig,
    sample_count: int,
    effective_workers: int,
    orf_genepred: str | Path,
    logger: _StageLogger | None = None,
) -> EngineResult:
    """Merge evidence and atomically export the reliable ORF annotation."""
    master_path = Path(str(output_prefix) + ".smorf_evidence.txt")
    reliable_path = Path(str(output_prefix) + ".reliable_smorf.txt")
    genepred_path = Path(str(output_prefix) + ".reliable_smorf.genepred")
    summary_path = Path(str(output_prefix) + ".evidence_summary.txt")
    master_raw_tmp = master_path.with_name(master_path.name + ".raw.tmp")
    master_tmp = master_path.with_name(master_path.name + ".tmp")
    reliable_tmp = reliable_path.with_name(reliable_path.name + ".tmp")
    genepred_tmp = genepred_path.with_name(genepred_path.name + ".tmp")
    summary_tmp = summary_path.with_name(summary_path.name + ".tmp")
    temporary_paths = (
        master_raw_tmp,
        master_tmp,
        reliable_tmp,
        genepred_tmp,
        summary_tmp,
    )
    reliable_columns = (
        "family_id", "gene_id", "chrom", "strand", "category",
        "family_type", "family_size", "structural_primary",
        "provisional_primary", "quant_primary", "evidence_primary",
        "transcript_id",
        "translation_unit_id", "translation_unit_index",
        "translation_unit_count", "translation_unit_status",
        "start_codon", "stop_codon",
        "nt_length", "aa_length", "exon_starts", "exon_ends",
        "family_translation_status", "translated_extent_status",
        "extent_supported_sample_count", "extent_supporting_samples",
        "extent_supported_bins", "start_site_status",
        "start_interval_orf_ids", "selection_policy", "selection_reason",
        "canonical_anchor_orf", "longest_candidate_orf",
        "prior_override_status", "leading_support_sample_count",
        "leading_negative_sample_count",
        "noncanonical_extension_support_sample_count",
        "noncanonical_extension_density_ratio",
        "pooled_extent_rpf_sum", "pooled_extent_frame0_ratio",
        "pooled_start_rpf_sum", "pooled_end_rpf_sum",
        "pooled_extent_sample_count",
        "nested_competition_status",
        "dominant_parent_orf", "competition_overlap_fraction",
        "competition_frame_relation", "competition_reason", "frame_margin",
        "supported_sample_count", "supporting_samples", "best_sample",
        "best_sample_evidence", "best_rpf_sum", "best_coverage_ratio",
        "best_frame0_ratio", "reliability_reason", "start_site_reason",
        "complete_candidate_count", "complete_candidate_orf_ids",
        "candidate_audit_status", "candidate_failure_reasons",
        "pooled_frame0_density", "pooled_frame1_density",
        "pooled_frame2_density", "pooled_dominant_frame",
        "suggested_psite_shift_nt", "phase_audit_status",
    )
    reliable_smorfs_written = 0
    annotated_controls_excluded = 0
    reliable_genepred_records = 0
    connection: sqlite3.Connection | None = None
    try:
        with gzip.open(
            master_raw_tmp,
            "wt",
            encoding="utf-8",
            compresslevel=1,
        ) as handle:
            handle.write("\t".join(_output_columns()) + "\n")
        with master_raw_tmp.open("ab") as target:
            for result in chromosome_results:
                with Path(result.part_path).open("rb") as source:
                    shutil.copyfileobj(
                        source,
                        target,
                        length=16 * 1024 * 1024,
                    )

        connection = sqlite3.connect(database_path)
        connection.row_factory = sqlite3.Row
        if logger is not None:
            logger.write(
                "Stage5b: Resolve transcript translation units and nested "
                "ORF competition."
            )
        competition_counts, translation_unit_counts = _apply_nested_competition(
            master_raw_tmp,
            master_tmp,
            connection,
            config,
        )
        connection.execute(
            """
            CREATE TEMP TABLE reliable_export (
                orf_id TEXT PRIMARY KEY,
                matched INTEGER NOT NULL DEFAULT 0
            ) WITHOUT ROWID
            """
        )
        with reliable_tmp.open(
            "w", encoding="utf-8", buffering=8 * 1024 * 1024
        ) as output_handle:
            output_handle.write("\t".join(reliable_columns) + "\n")
            with master_tmp.open("r", encoding="utf-8") as input_handle:
                reader = csv.DictReader(input_handle, delimiter="\t")
                for row in reader:
                    if not (
                        row["evidence_status"] == "Reliable"
                        and row["translated_extent_status"]
                        in {"Supported", "Compatible"}
                        and row["quant_primary"]
                        and row["nested_competition_status"]
                        not in {
                            "ShadowedByLongORF",
                            "ShadowedByHighOverlapORF",
                            "ShadowedByCanonicalATG",
                        }
                    ):
                        continue
                    if _is_annotated_category(row.get("category", "")):
                        annotated_controls_excluded += 1
                        continue
                    geometry = connection.execute(
                        """
                        SELECT transcript_id, category, start_codon, stop_codon,
                               nt_length, aa_length, exon_starts, exon_ends
                        FROM geometry WHERE orf_id=?
                        """,
                        (row["quant_primary"],),
                    ).fetchone()
                    if geometry is None:
                        raise EvidenceEngineError(
                            "Reliable quantification primary is absent from "
                            f"the geometry index: {row['quant_primary']}"
                        )
                    if _is_annotated_category(geometry["category"]):
                        annotated_controls_excluded += 1
                        continue
                    try:
                        connection.execute(
                            "INSERT INTO reliable_export (orf_id) VALUES (?)",
                            (row["quant_primary"],),
                        )
                    except sqlite3.IntegrityError as error:
                        raise EvidenceEngineError(
                            "One ORF is the quantification primary of multiple "
                            "reliable families: "
                            f"{row['quant_primary']}"
                        ) from error
                    output = {
                        **row,
                        "transcript_id": geometry["transcript_id"],
                        "start_codon": geometry["start_codon"],
                        "stop_codon": geometry["stop_codon"],
                        "nt_length": geometry["nt_length"],
                        "aa_length": geometry["aa_length"],
                        "exon_starts": geometry["exon_starts"],
                        "exon_ends": geometry["exon_ends"],
                    }
                    output_handle.write(
                        "\t".join(
                            _format_value(output.get(column, ""))
                            for column in reliable_columns
                        )
                        + "\n"
                    )
                    reliable_smorfs_written += 1
        connection.commit()

        if logger is not None:
            logger.write(
                "Stage6: Export reliable smORF genePred annotation."
            )
        reliable_genepred_records = _export_reliable_genepred(
            source_path=orf_genepred,
            temporary_path=genepred_tmp,
            connection=connection,
        )
        if reliable_genepred_records != reliable_smorfs_written:
            raise EvidenceEngineError(
                "Reliable output and genePred record counts differ: "
                f"table={reliable_smorfs_written:,}, "
                f"genePred={reliable_genepred_records:,}."
            )

        totals = {
            "total_families": sum(
                item.total_families for item in chromosome_results
            ),
            "reliable_families": sum(
                item.reliable_families for item in chromosome_results
            ),
            "uncertain_families": sum(
                item.uncertain_families for item in chromosome_results
            ),
            "no_evidence_families": sum(
                item.no_evidence_families for item in chromosome_results
            ),
            "reliable_smorfs": reliable_smorfs_written,
            "invalid_families": sum(
                item.invalid_families for item in chromosome_results
            ),
        }
        evidence_counts: dict[str, int] = defaultdict(int)
        reason_counts: dict[str, int] = defaultdict(int)
        phase_audit_counts: dict[str, int] = defaultdict(int)
        annotated_phase_shift_counts: dict[int, int] = defaultdict(int)
        replicated_signal_families = 0
        with master_tmp.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                evidence_counts[row.get("best_sample_evidence", "")] += 1
                reason_counts[row.get("reliability_reason", "")] += 1
                phase_audit_counts[row.get("phase_audit_status", "")] += 1
                if (
                    _is_annotated_category(row.get("category", ""))
                    and row.get("phase_audit_status")
                    == "AlternativeFrameDominant"
                ):
                    annotated_phase_shift_counts[
                        int(row.get("suggested_psite_shift_nt", "0") or 0)
                    ] += 1
                try:
                    supported = int(
                        row.get("supported_sample_count", "0") or 0
                    )
                    if supported >= config.reliable_sample:
                        replicated_signal_families += 1
                except ValueError:
                    pass

        with summary_tmp.open("w", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write(f"sample_count\t{sample_count}\n")
            handle.write(f"effective_workers\t{effective_workers}\n")
            for key, value in totals.items():
                handle.write(f"{key}\t{value}\n")
            handle.write(
                "reliable_genepred_records\t"
                f"{reliable_genepred_records}\n"
            )
            handle.write(
                "annotated_controls_excluded_from_reliable_smorf\t"
                f"{annotated_controls_excluded}\n"
            )
            handle.write(
                "reliable_family_definition\t"
                "Candidate-first distributed extent evidence or Medium/High "
                "common-body evidence in at least "
                f"{config.reliable_sample} independent samples\n"
            )
            handle.write(
                "reliable_smorf_definition\t"
                "Reliable family with a supported or prior-compatible "
                "translated extent, a quant_primary, and no dominant ORF "
                "within the same transcript translation unit\n"
            )
            handle.write(
                "uncertain_definition\t"
                "Signal without sufficient independent-sample replication, "
                "localized long-ORF signal, or invalid family geometry\n"
            )
            handle.write(f"evidence_mode\t{config.evidence_mode}\n")
            handle.write(f"group_column\t{config.group_column}\n")
            handle.write(f"reliable_sample\t{config.reliable_sample}\n")
            handle.write(
                "estimated_short_max_codons\t"
                f"{config.estimated_short_max_codons}\n"
            )
            handle.write(
                "estimated_long_min_codons\t"
                f"{config.estimated_long_min_codons}\n"
            )
            handle.write(f"short_max_codons\t{config.short_max_codons}\n")
            handle.write(f"long_min_codons\t{config.long_min_codons}\n")
            handle.write(f"window_codons\t{config.window_codons}\n")
            handle.write(f"window_step_codons\t{config.window_step_codons}\n")
            handle.write(f"min_supported_windows\t{config.min_supported_windows}\n")
            handle.write(f"min_window_gap_codons\t{config.min_window_gap_codons}\n")
            handle.write(f"min_signal_span\t{config.min_signal_span}\n")
            handle.write(f"localized_span_max\t{config.localized_span_max}\n")
            handle.write(
                "localized_top_window_fraction\t"
                f"{config.localized_top_window_fraction}\n"
            )
            handle.write(f"boundary_codons\t{config.boundary_codons}\n")
            handle.write(f"extent_bins\t{config.extent_bins}\n")
            handle.write(f"min_extent_bins\t{config.min_extent_bins}\n")
            handle.write(
                f"start_resolution_codons\t"
                f"{config.start_resolution_codons}\n"
            )
            handle.write(
                f"leading_window_codons\t{config.leading_window_codons}\n"
            )
            handle.write(
                f"min_exclusion_codons\t{config.min_exclusion_codons}\n"
            )
            handle.write(
                "silent_extension_density_ratio\t"
                f"{config.silent_extension_density_ratio}\n"
            )
            handle.write(
                "noncanonical_override_density_ratio\t"
                f"{config.noncanonical_override_density_ratio}\n"
            )
            handle.write(
                "noncanonical_override_min_coverage_ratio\t"
                f"{config.noncanonical_override_min_coverage_ratio}\n"
            )
            handle.write(
                "noncanonical_override_min_frame_margin\t"
                f"{config.noncanonical_override_min_frame_margin}\n"
            )
            handle.write(
                "noncanonical_min_exclusive_codons\t"
                f"{config.noncanonical_min_exclusive_codons}\n"
            )
            handle.write(
                f"nested_min_frame_margin\t"
                f"{config.nested_min_frame_margin}\n"
            )
            handle.write(
                f"nested_min_phase_rpf\t{config.nested_min_phase_rpf}\n"
            )
            handle.write(
                f"high_overlap_fraction\t{config.high_overlap_fraction}\n"
            )
            handle.write(
                "overlap_min_frame_margin\t"
                f"{config.overlap_min_frame_margin}\n"
            )
            handle.write(
                f"overlap_min_phase_rpf\t{config.overlap_min_phase_rpf}\n"
            )
            for label, count in sorted(competition_counts.items()):
                handle.write(f"nested_competition_{label}\t{count}\n")
            handle.write(
                "translation_unit_total\t"
                f"{sum(translation_unit_counts.values())}\n"
            )
            for label, count in sorted(translation_unit_counts.items()):
                handle.write(f"translation_unit_{label}\t{count}\n")
            handle.write(
                f"families_with_replicated_medium_high_support\t"
                f"{replicated_signal_families}\n"
            )
            for label, count in sorted(evidence_counts.items()):
                handle.write(
                    f"best_sample_evidence_{label or 'NA'}\t{count}\n"
                )
            for label, count in sorted(phase_audit_counts.items()):
                handle.write(
                    f"phase_audit_{label or 'NA'}\t{count}\n"
                )
            for shift, count in sorted(annotated_phase_shift_counts.items()):
                handle.write(
                    f"annotated_control_suggested_psite_shift_{shift:+d}_nt\t"
                    f"{count}\n"
                )
            for reason, count in sorted(reason_counts.items()):
                handle.write(
                    f"reliability_reason_{reason or 'NA'}\t{count}\n"
                )
        if logger is not None and annotated_phase_shift_counts:
            dominant_shift, dominant_count = max(
                annotated_phase_shift_counts.items(),
                key=lambda item: item[1],
            )
            if dominant_count >= config.positive_min_controls:
                warning_print(
                    "Annotated ORF phase audit detected a systematic "
                    f"alternative frame: suggested P-site coordinate shift "
                    f"{dominant_shift:+d} nt in {dominant_count:,} controls. "
                    "Verify the upstream P-site offset before interpreting "
                    "candidate reading frames."
                )
    except BaseException:
        for path in temporary_paths:
            path.unlink(missing_ok=True)
        raise
    finally:
        if connection is not None:
            connection.close()

    for temporary, final in (
        (master_tmp, master_path),
        (reliable_tmp, reliable_path),
        (genepred_tmp, genepred_path),
        (summary_tmp, summary_path),
    ):
        os.replace(temporary, final)
    master_raw_tmp.unlink(missing_ok=True)

    return EngineResult(
        master_output=str(master_path),
        reliable_output=str(reliable_path),
        reliable_genepred_output=str(genepred_path),
        summary_output=str(summary_path),
        sample_count=sample_count,
        effective_workers=effective_workers,
        **{key: totals[key] for key in (
            "total_families", "reliable_families", "uncertain_families",
            "no_evidence_families", "reliable_smorfs",
        )},
    )


def _build_config(args: object) -> EngineConfig:
    """Build an engine configuration from public and hidden arguments."""
    hidden_fields = (
        "short_max_codons",
        "long_min_codons",
        "window_codons",
        "window_step_codons",
        "min_supported_windows",
        "min_window_gap_codons",
        "min_signal_span",
        "localized_span_max",
        "localized_top_window_fraction",
        "boundary_codons",
        "extent_bins",
        "min_extent_bins",
        "start_resolution_codons",
        "min_frame_margin",
        "leading_window_codons",
        "min_exclusion_codons",
        "silent_extension_density_ratio",
        "noncanonical_min_exclusive_codons",
        "noncanonical_override_density_ratio",
        "noncanonical_override_min_coverage_ratio",
        "noncanonical_override_min_frame_margin",
        "nested_min_frame_margin",
        "nested_min_phase_rpf",
        "high_overlap_fraction",
        "overlap_min_frame_margin",
        "overlap_min_phase_rpf",
    )
    overrides = frozenset(
        field
        for field in hidden_fields
        if getattr(args, field, None) is not None
    )

    def advanced(name: str) -> float | int:
        value = getattr(args, name, None)
        return ADVANCED_DEFAULTS[name] if value is None else value

    return EngineConfig(
        evidence_mode=str(args.evidence_mode),
        group_column=str(args.group_column),
        reliable_sample=int(args.reliable_sample),
        min_rpf_sum=float(args.min_rpf_sum),
        min_rpf_per_codon=float(args.min_rpf_per_codon),
        min_covered_codon=int(args.min_covered_codon),
        min_coverage_ratio=float(args.min_codon_coverage),
        moderate_periodicity=float(args.moderate_periodicity),
        strong_periodicity=float(args.strong_periodicity),
        min_window_rpf=float(args.min_window_rpf),
        min_window_covered=int(
            args.min_window_covered_codon
        ),
        short_max_codons=int(advanced("short_max_codons")),
        long_min_codons=int(advanced("long_min_codons")),
        window_codons=int(advanced("window_codons")),
        window_step_codons=int(advanced("window_step_codons")),
        min_supported_windows=int(
            advanced("min_supported_windows")
        ),
        min_window_gap_codons=int(
            advanced("min_window_gap_codons")
        ),
        min_signal_span=float(advanced("min_signal_span")),
        localized_span_max=float(
            advanced("localized_span_max")
        ),
        localized_top_window_fraction=float(
            advanced("localized_top_window_fraction")
        ),
        boundary_codons=int(advanced("boundary_codons")),
        extent_bins=int(advanced("extent_bins")),
        min_extent_bins=int(advanced("min_extent_bins")),
        start_resolution_codons=int(
            advanced("start_resolution_codons")
        ),
        min_frame_margin=float(advanced("min_frame_margin")),
        leading_window_codons=int(advanced("leading_window_codons")),
        min_exclusion_codons=int(advanced("min_exclusion_codons")),
        silent_extension_density_ratio=float(
            advanced("silent_extension_density_ratio")
        ),
        noncanonical_min_exclusive_codons=int(
            advanced("noncanonical_min_exclusive_codons")
        ),
        noncanonical_override_density_ratio=float(
            advanced("noncanonical_override_density_ratio")
        ),
        noncanonical_override_min_coverage_ratio=float(
            advanced("noncanonical_override_min_coverage_ratio")
        ),
        noncanonical_override_min_frame_margin=float(
            advanced("noncanonical_override_min_frame_margin")
        ),
        nested_min_frame_margin=float(
            advanced("nested_min_frame_margin")
        ),
        nested_min_phase_rpf=float(
            advanced("nested_min_phase_rpf")
        ),
        high_overlap_fraction=float(
            advanced("high_overlap_fraction")
        ),
        overlap_min_frame_margin=float(
            advanced("overlap_min_frame_margin")
        ),
        overlap_min_phase_rpf=float(
            advanced("overlap_min_phase_rpf")
        ),
        hidden_overrides=overrides,
        positive_quantile=float(args.positive_quantile),
        positive_min_controls=int(args.positive_min_controls),
        positive_max_controls=int(args.positive_max_controls),
        threads=int(args.threads),
    )



def run_family_evidence_engine(args: object) -> EngineResult:
    """Run the complete bottom-up family evidence engine."""
    tracks = read_density_tracks(args)
    config = _build_config(args)
    output_prefix = Path(args.output).expanduser().resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    work_directory = Path(str(output_prefix) + ".evidence_work")
    logger = _StageLogger()
    signature = _run_signature(args, tracks)

    if work_directory.exists():
        shutil.rmtree(work_directory)
    work_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = work_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {"signature": signature},
            indent=2,
        ),
        encoding="utf-8",
    )

    database_path = work_directory / "family_index.sqlite"
    try:
        logger.write("Stage1: Build or reuse the family SQLite index.")
        build_family_database(
            args.family_table,
            args.family_members,
            args.orf_source,
            database_path,
            signature,
            logger,
        )

        config = resolve_length_models(
            database_path,
            config,
            logger,
        )

        logger.write("Stage2: Parse every density track once and build chromosome caches.")
        for number, track in enumerate(tracks, start=1):
            logger.write(
                f"Density track {number}/{len(tracks)}: sample={track.sample}, strand={track.strand}."
            )
            build_density_cache(track, work_directory, signature)

        logger.write("Stage3: Calibrate every sample with annotated ORFs.")
        thresholds, config = calibrate_samples(
            database_path,
            work_directory,
            tracks,
            config,
            logger,
        )

        logger.write("Stage4: Evaluate chromosome family scaffolds.")
        chromosomes = _chromosomes(database_path)
        groups = tuple(sorted({track.group for track in tracks}))
        tasks = [
            ChromosomeTask(
                chromosome=chrom,
                database_path=str(database_path),
                work_directory=str(work_directory),
                output_prefix=str(output_prefix),
                tracks=tuple(tracks),
                thresholds=tuple(thresholds),
                config=config,
                group_names=groups,
            )
            for chrom in chromosomes
        ]
        workers = _effective_workers(config.threads, len(tasks))
        logger.write(
            f"Chromosomes={len(tasks)}, requested_workers={config.threads}, effective_workers={workers}."
        )
        results_by_chrom: dict[str, ChromosomeResult] = {}
        if workers == 1:
            for number, task in enumerate(tasks, start=1):
                result = _process_chromosome(task)
                results_by_chrom[result.chromosome] = result
                progress_print(f"evidence chromosomes: {number}/{len(tasks)}")
        else:
            methods = mp.get_all_start_methods()
            context = mp.get_context("fork") if "fork" in methods else mp.get_context()
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
                futures = {executor.submit(_process_chromosome, task): task.chromosome for task in tasks}
                completed = 0
                for future in as_completed(futures):
                    chrom = futures[future]
                    try:
                        result = future.result()
                    except Exception as error:
                        raise EvidenceEngineError(
                            f"Chromosome task failed: {chrom}: {error}"
                        ) from error
                    results_by_chrom[result.chromosome] = result
                    completed += 1
                    progress_print(f"evidence chromosomes: {completed}/{len(tasks)}")
        ordered_results = [results_by_chrom[chrom] for chrom in chromosomes]

        logger.write("Stage5: Merge focused family and reliable-smORF outputs.")
        result = _merge_outputs(
            output_prefix,
            ordered_results,
            database_path,
            config,
            sample_count=len({track.sample for track in tracks}),
            effective_workers=workers,
            orf_genepred=args.orf_genepred,
            logger=logger,
        )
        logger.write(
            "Completed: total={total:,}, reliable={reliable:,}, uncertain={uncertain:,}, no_evidence={none:,}.".format(
                total=result.total_families,
                reliable=result.reliable_families,
                uncertain=result.uncertain_families,
                none=result.no_evidence_families,
            )
        )
        shutil.rmtree(work_directory, ignore_errors=True)
        return result
    except Exception as error:
        logger.error("run_family_evidence_engine", error)
        shutil.rmtree(work_directory, ignore_errors=True)
        raise
