#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Resolve translated extents and candidate start sites.
# Input: Family accumulators and evidence features.
# Output: Resolved primary ORFs and extent audit results.

"""Resolve translated extents and candidate start sites."""

from __future__ import annotations

import math
from typing import Sequence

from .features import (
    _candidate_boundary_support,
    _candidate_cohort_extent_status,
    _candidate_cohort_sample_count,
    _candidate_is_complete,
    _candidate_segmented_frame_density,
    _candidate_supported_bin_indices,
    _canonical_anchor_index,
    _frame_margin,
)
from .models import (
    EngineConfig,
    ExtentResolution,
    FamilyGeometry,
    _FamilyAccumulator,
)


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
                [0]
                * min(
                    family.representatives[index].coding_nt_length // 3,
                    config.extent_bins,
                )
                if family.representatives[index].coding_nt_length // 3 >= config.long_min_codons
                else []
            )
            for index in range(candidate_count)
        ],
        candidate_bin_support_samples=[[] for _ in range(candidate_count)],
        candidate_cohort_support_samples=[[] for _ in range(candidate_count)],
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
                if family.representatives[index].coding_nt_length // 3 >= config.long_min_codons
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
                if family.representatives[index].coding_nt_length // 3 >= config.long_min_codons
                else []
            )
            for index in range(candidate_count)
        ],
        candidate_best_full_features=[None] * candidate_count,
        candidate_best_extension_features=[None] * candidate_count,
        candidate_pooled_frame_density=[[0.0, 0.0, 0.0] for _ in range(candidate_count)],
        candidate_leading_support_counts=[0] * candidate_count,
        candidate_leading_negative_counts=[0] * candidate_count,
        candidate_noncanonical_extension_support_counts=[0] * candidate_count,
        candidate_best_noncanonical_density_ratios=[0.0] * candidate_count,
    )


def _output_columns() -> tuple[str, ...]:
    """Return focused family output columns."""
    return (
        "family_id",
        "gene_id",
        "chrom",
        "strand",
        "category",
        "family_type",
        "family_size",
        "structural_primary",
        "common_body_orf",
        "provisional_primary",
        "quant_primary",
        "evidence_primary",
        "translation_unit_id",
        "translation_unit_index",
        "translation_unit_count",
        "translation_unit_status",
        "evidence_status",
        "family_translation_status",
        "translated_extent_status",
        "extent_supported_sample_count",
        "extent_supporting_samples",
        "extent_supported_bins",
        "extension_rpf_sum",
        "extension_frame0_ratio",
        "family_translation_evidence",
        "reliability_reason",
        "start_site_status",
        "start_site_reason",
        "start_interval_orf_ids",
        "selection_policy",
        "selection_reason",
        "canonical_anchor_orf",
        "longest_candidate_orf",
        "prior_override_status",
        "leading_support_sample_count",
        "leading_negative_sample_count",
        "noncanonical_extension_support_sample_count",
        "noncanonical_extension_density_ratio",
        "pooled_extent_rpf_sum",
        "pooled_extent_frame0_ratio",
        "pooled_start_rpf_sum",
        "pooled_end_rpf_sum",
        "pooled_extent_sample_count",
        "nested_competition_status",
        "dominant_parent_orf",
        "competition_overlap_fraction",
        "competition_frame_relation",
        "competition_reason",
        "frame_margin",
        "reliable_group",
        "supported_sample_count",
        "max_group_sample_support",
        "high_confidence_sample_count",
        "supporting_group_count",
        "supporting_groups",
        "supporting_samples",
        "best_sample",
        "best_sample_evidence",
        "best_rpf_sum",
        "best_rpf_per_codon",
        "best_covered_codon",
        "best_coverage_ratio",
        "best_frame0_ratio",
        "best_supported_windows",
        "best_distributed_windows",
        "best_signal_span",
        "best_top_window_fraction",
        "best_localized_only",
        "representative_count",
        "representative_orf_ids",
        "representative_start_codons",
        "invalid_family_reason",
        "complete_candidate_count",
        "complete_candidate_orf_ids",
        "candidate_audit_status",
        "candidate_failure_reasons",
        "pooled_frame0_density",
        "pooled_frame1_density",
        "pooled_frame2_density",
        "pooled_dominant_frame",
        "suggested_psite_shift_nt",
        "phase_audit_status",
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
    if accumulator.candidate_full_support_counts[candidate_index] >= config.reliable_sample:
        return True
    if (
        codon_count < config.long_min_codons
        and _frame_margin(accumulator.candidate_pooled_frame_density[candidate_index])
        < config.min_frame_margin
    ):
        return False
    extension_codons = max(
        0,
        (family.common_body_start - representative.start_offset) // 3,
    )
    extension_supported = (
        extension_codons <= config.start_resolution_codons
        or accumulator.candidate_extension_support_counts[candidate_index] >= config.reliable_sample
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
            f"candidate_sample_replication_{sample_count}_below_{config.reliable_sample}"
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
            reasons.append(f"supported_extent_bins_{len(supported_bins)}_below_{minimum_bins}")
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
        if supported_bins and supported_bins[-1] == len(bin_counts) - 1 and not end_supported:
            reasons.append("missing_candidate_terminal_phase_support")
        segmented_frames = _candidate_segmented_frame_density(
            accumulator,
            candidate_index,
            supported_bins,
        )
        if supported_bins and _frame_margin(segmented_frames) < config.min_frame_margin:
            reasons.append("candidate_segmented_frame_margin_below_threshold")
    elif sum(pooled_frames) > 0 and _frame_margin(pooled_frames) < config.min_frame_margin:
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


def _resolve_candidate_extent(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    config: EngineConfig,
) -> ExtentResolution:
    """Select an ORF using canonical/length priors plus replicated evidence."""
    selection_policy = "canonical_ATG_then_longest_unless_replicated_RPF_override"
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
            accumulator.candidate_leading_negative_counts[index] >= config.reliable_sample
            and accumulator.candidate_leading_support_counts[index] < config.reliable_sample
        )

    selected_index = prior_index
    if canonical_index is not None:
        strong_upstream_noncanonical = [
            index
            for index, representative in enumerate(family.representatives)
            if representative.start_offset < family.representatives[canonical_index].start_offset
            and representative.start_codon.upper() != "ATG"
            and accumulator.candidate_noncanonical_extension_support_counts[index]
            >= config.reliable_sample
            and (
                accumulator.candidate_leading_support_counts[index] >= config.reliable_sample
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
                "upstream_noncanonical_start_has_strong_exclusive_replicated_translation"
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
                "longer_prior_rejected_by_replicated_N_terminal_silence_and_internal_extent_support"
            )
            override_status = "ShortenedByReplicatedNterminalSilence"
        else:
            diagnostic_frames = accumulator.candidate_pooled_frame_density[prior_index]
            diagnostic_total = float(sum(diagnostic_frames))
            diagnostic_boundary = accumulator.candidate_boundary_pooled_stats[prior_index]
            diagnostic_extension = accumulator.candidate_best_extension_features[prior_index]
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
                    0.0 if diagnostic_extension is None else diagnostic_extension.rpf_sum
                ),
                extension_frame0_ratio=(
                    0.0 if diagnostic_extension is None else diagnostic_extension.frame0_ratio
                ),
                start_site_status="Unresolved",
                start_site_reason=(
                    "longer_prior_contradicted_without_supported_internal_replacement"
                ),
                start_interval_orf_ids=tuple(
                    family.representatives[index].orf_id for index in supported_indices
                ),
                frame_margin=_frame_margin(diagnostic_frames),
                selection_policy=selection_policy,
                selection_reason=("replicated_N_terminal_silence_without_replacement"),
                canonical_anchor_orf=(
                    ""
                    if canonical_index is None
                    else family.representatives[canonical_index].orf_id
                ),
                longest_candidate_orf=(family.representatives[longest_index].orf_id),
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
                    diagnostic_frames[0] / diagnostic_total if diagnostic_total > 0 else 0.0
                ),
                pooled_start_rpf_sum=(
                    float(diagnostic_boundary[0][0]) if diagnostic_boundary else 0.0
                ),
                pooled_end_rpf_sum=(
                    float(diagnostic_boundary[1][0]) if len(diagnostic_boundary) > 1 else 0.0
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
            + [family.representatives[index].orf_id for index in supported_indices]
        )
    )
    nearby = [
        item.orf_id
        for item in family.representatives
        if abs(item.start_offset - selected.start_offset) <= config.start_resolution_codons * 3
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
    extension_total = 0.0 if extension_features is None else extension_features.rpf_sum
    extension_frame0 = 0.0 if extension_features is None else extension_features.frame0_ratio
    pooled_total = float(sum(pooled_frames))
    boundary = accumulator.candidate_boundary_pooled_stats[selected_index]
    return ExtentResolution(
        quant_primary=selected.orf_id,
        evidence_primary=evidence_primary,
        translated_extent_status=(
            "Supported"
            if selected_is_directly_supported
            else (cohort_extent_status if cohort_extent_status else "Compatible")
        ),
        extent_supported_sample_count=max(
            accumulator.candidate_full_support_counts[selected_index],
            accumulator.candidate_extension_support_counts[selected_index],
            accumulator.candidate_leading_support_counts[selected_index],
            len(set(accumulator.candidate_bin_support_samples[selected_index])),
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
            "" if canonical_index is None else family.representatives[canonical_index].orf_id
        ),
        longest_candidate_orf=family.representatives[longest_index].orf_id,
        prior_override_status=override_status,
        leading_support_sample_count=(accumulator.candidate_leading_support_counts[selected_index]),
        leading_negative_sample_count=(
            accumulator.candidate_leading_negative_counts[selected_index]
        ),
        noncanonical_extension_support_sample_count=(
            accumulator.candidate_noncanonical_extension_support_counts[selected_index]
        ),
        noncanonical_extension_density_ratio=(
            accumulator.candidate_best_noncanonical_density_ratios[selected_index]
        ),
        pooled_extent_rpf_sum=pooled_total,
        pooled_extent_frame0_ratio=(pooled_frames[0] / pooled_total if pooled_total > 0 else 0.0),
        pooled_start_rpf_sum=(float(boundary[0][0]) if boundary else 0.0),
        pooled_end_rpf_sum=(float(boundary[1][0]) if len(boundary) > 1 else 0.0),
        pooled_extent_sample_count=_candidate_cohort_sample_count(
            accumulator,
            selected_index,
        ),
    )
