#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Calculate and classify family Ribo-seq features.
# Input: Family geometry and sparse P-site signal.
# Output: Coverage, periodicity, window, and extent features.

"""Calculate and classify family Ribo-seq features."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Mapping, Sequence

import numpy as np

from utils.smorf.density import (
    ChromDensity,
)

from .config import (
    LEVEL_RANK,
    STOP_CODONS,
)
from .input import _oriented_blocks
from .models import (
    DensityTrack,
    EngineConfig,
    FamilyGeometry,
    RepresentativeGeometry,
    ScaffoldProfile,
    SegmentFeatures,
    Thresholds,
    _FamilyAccumulator,
)


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
        support = features.supported_windows >= 1 or features.coverage_ratio >= max(
            0.20, thresholds.min_coverage_ratio * 1.5
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
    return features.distributed_windows >= 1 or features.signal_span >= max(
        0.20, config.min_signal_span * 0.75
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
    return relative_density <= config.silent_extension_density_ratio and (
        weak_abundance or weak_phase
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
        and canonical_body.rpf_per_codon >= thresholds.min_rpf_per_codon
        and canonical_body.frame0_ratio >= thresholds.moderate_periodicity
    ):
        return False
    density_floor = max(
        thresholds.min_rpf_per_codon,
        canonical_body.rpf_per_codon * config.noncanonical_override_density_ratio,
    )
    if not (
        extension.rpf_sum
        >= max(
            thresholds.min_rpf_sum,
            thresholds.min_window_rpf,
        )
        and extension.rpf_per_codon >= density_floor
        and extension.covered_codon >= thresholds.min_window_covered
        and extension.coverage_ratio
        >= max(
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
                extension.frame0_density + extension.frame1_density + extension.frame2_density,
                1e-12,
            )
            >= config.noncanonical_override_min_frame_margin
        )
        and not extension.localized_only
    ):
        return False
    if extension_codons < config.long_min_codons:
        return True
    return extension.distributed_windows >= 1 or extension.signal_span >= max(
        0.20, config.min_signal_span * 0.75
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
    boundaries = [(index * codons) // effective_bins for index in range(effective_bins + 1)]
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
        and features.covered_codon >= max(1, math.ceil(thresholds.min_window_covered * 0.5))
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
        "",
        "complete",
        "cmpl",
        "full",
        "true",
        "yes",
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
            or (index < len(pooled) and _pooled_segment_is_supported(pooled[index], config))
        )
    ]


def _candidate_boundary_support(
    accumulator: _FamilyAccumulator,
    candidate_index: int,
    config: EngineConfig,
) -> tuple[bool, bool]:
    """Return pooled start- and stop-window support for one candidate."""
    boundary = accumulator.candidate_boundary_pooled_stats[candidate_index]
    start_supported = bool(boundary and _pooled_segment_is_supported(boundary[0], config))
    end_supported = bool(len(boundary) > 1 and _pooled_segment_is_supported(boundary[1], config))
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

    span_fraction = (supported_bins[-1] - supported_bins[0] + 1) / len(bin_counts)
    if supported_bins[-1] == len(bin_counts) - 1 and end_supported:
        return "Supported"
    downstream_start = max(1, math.floor(len(bin_counts) * 0.60))
    if span_fraction >= max(0.60, config.min_signal_span) and any(
        index >= downstream_start for index in supported_bins
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
    indices = range(len(pooled)) if supported_bins is None else supported_bins
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


def _frame_margin(frame_density: Sequence[float]) -> float:
    """Calculate the candidate frame-0 advantage over the alternative frames."""
    total = float(sum(frame_density))
    if total <= 0:
        return 0.0
    ratios = [float(value) / total for value in frame_density]
    return ratios[0] - max(ratios[1], ratios[2])
