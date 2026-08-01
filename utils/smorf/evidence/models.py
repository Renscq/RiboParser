#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Define family-evidence data models.
# Input: Evidence configuration and density arrays.
# Output: Typed geometry, feature, and result records.

"""Define family-evidence data models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from utils.ribo.ArgsParser import message_print, warning_print


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

        rpf_sum = float(self.codon_prefix[end_codon] - self.codon_prefix[start_codon])
        if rpf_sum <= 0:
            return _empty_segment_features()
        covered = int(self.covered_prefix[end_codon] - self.covered_prefix[start_codon])
        coverage_ratio = covered / codon_count
        local_phase = int(start) % 3
        frame_order = (
            local_phase,
            (local_phase + 1) % 3,
            (local_phase + 2) % 3,
        )
        frame_density = [
            float(self.frame_prefixes[frame][end_codon] - self.frame_prefixes[frame][start_codon])
            for frame in frame_order
        ]
        frame_total = sum(frame_density)
        frame0_ratio = frame_density[0] / frame_total if frame_total else 0.0

        nonzero_left = int(np.searchsorted(self.nonzero_codons, start_codon, side="left"))
        nonzero_right = int(np.searchsorted(self.nonzero_codons, end_codon, side="left"))
        if nonzero_right > nonzero_left:
            first_nonzero = int(self.nonzero_codons[nonzero_left])
            last_nonzero = int(self.nonzero_codons[nonzero_right - 1])
            signal_span = (last_nonzero - first_nonzero + 1) / codon_count
        else:
            signal_span = 0.0

        boundary = min(config.boundary_codons, codon_count)
        start_rpf = float(
            self.codon_prefix[start_codon + boundary] - self.codon_prefix[start_codon]
        )
        end_rpf = float(self.codon_prefix[end_codon] - self.codon_prefix[end_codon - boundary])
        body_start = start_codon + boundary
        body_end = max(body_start, end_codon - boundary)
        body_rpf = float(self.codon_prefix[body_end] - self.codon_prefix[body_start])

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
                window_sum = float(self.codon_prefix[right] - self.codon_prefix[left])
                top_window = max(top_window, window_sum)
                window_covered = int(self.covered_prefix[right] - self.covered_prefix[left])
                window_frames = [
                    float(self.frame_prefixes[frame][right] - self.frame_prefixes[frame][left])
                    for frame in frame_order
                ]
                window_total = sum(window_frames)
                window_frame0 = window_frames[0] / window_total if window_total else 0.0
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


def _empty_segment_features() -> SegmentFeatures:
    """Return an all-zero segment feature record."""
    return SegmentFeatures(
        0.0,
        0.0,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        False,
        0.0,
    )
