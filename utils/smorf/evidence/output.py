#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Build and write family-evidence outputs.
# Input: Evaluated chromosome families.
# Output: Evidence rows and reliable smORF genePred.

"""Build and write family-evidence outputs."""

from __future__ import annotations

import gzip
import json
import math
import os
import sqlite3
from collections import defaultdict
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence, TextIO

from utils.ribo.ArgsParser import progress_print
from utils.smorf.density import (
    ChromDensity,
)

from .config import (
    ANNOTATED_CATEGORIES,
    LEVEL_RANK,
)
from .extent import (
    _candidate_extent_supported,
    _candidate_failure_reasons,
    _new_accumulator,
    _output_columns,
    _phase_audit,
    _prior_candidate_index,
    _resolve_candidate_extent,
)
from .features import (
    _accumulate_pooled_segment,
    _bin_is_supported,
    _candidate_bin_ranges,
    _candidate_cohort_extent_status,
    _candidate_cohort_sample_count,
    _candidate_is_complete,
    _canonical_anchor_index,
    _extract_sparse_signals,
    _frame_margin,
    _leading_region_end,
    _leading_region_is_negative,
    _leading_region_is_supported,
    _noncanonical_extension_is_strong,
    _track_for_family,
    classify_candidate_extent,
    classify_sample,
)
from .input import (
    _iter_chromosome_family_batches,
    _load_cached_density,
    _safe_name,
    _smart_open,
)
from .models import (
    ChromosomeResult,
    ChromosomeTask,
    DensityTrack,
    EngineConfig,
    EvidenceEngineError,
    FamilyGeometry,
    InvalidFamily,
    ScaffoldProfile,
    _FamilyAccumulator,
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
    supporting_group_indices = sorted(accumulator.group_support_counts)
    supporting_groups = ",".join(str(group_names[index]) for index in supporting_group_indices)
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
    common_body_reliable = accumulator.support_count >= config.reliable_sample
    candidate_extent_reliable = bool(candidate_extent_indices)
    cohort_extent_reliable = bool(cohort_extent_candidates)

    if common_body_reliable or candidate_extent_reliable or cohort_extent_reliable:
        status = "Reliable"
        final_evidence = (
            "HighConfidence"
            if common_body_reliable and accumulator.high_count >= config.reliable_sample
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
                else [(index, "Supported") for index in candidate_extent_indices]
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
            reliable_group = str(group_names[supporting_group_indices[0]])
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
        item.orf_id for item in family.representatives if _candidate_is_complete(item)
    )
    failure_reasons = _candidate_failure_reasons(
        family,
        accumulator,
        selected_index,
        config,
    )
    selected_frames = accumulator.candidate_pooled_frame_density[selected_index]
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
        "extent_supported_sample_count": (resolution.extent_supported_sample_count),
        "extent_supporting_samples": ",".join(resolution.extent_supporting_samples),
        "extent_supported_bins": ",".join(str(value) for value in resolution.extent_supported_bins),
        "extension_rpf_sum": resolution.extension_rpf_sum,
        "extension_frame0_ratio": resolution.extension_frame0_ratio,
        "family_translation_evidence": final_evidence,
        "reliability_reason": reason,
        "start_site_status": resolution.start_site_status,
        "start_site_reason": resolution.start_site_reason,
        "start_interval_orf_ids": ",".join(resolution.start_interval_orf_ids),
        "selection_policy": resolution.selection_policy,
        "selection_reason": resolution.selection_reason,
        "canonical_anchor_orf": resolution.canonical_anchor_orf,
        "longest_candidate_orf": resolution.longest_candidate_orf,
        "prior_override_status": resolution.prior_override_status,
        "leading_support_sample_count": (resolution.leading_support_sample_count),
        "leading_negative_sample_count": (resolution.leading_negative_sample_count),
        "noncanonical_extension_support_sample_count": (
            resolution.noncanonical_extension_support_sample_count
        ),
        "noncanonical_extension_density_ratio": (resolution.noncanonical_extension_density_ratio),
        "pooled_extent_rpf_sum": resolution.pooled_extent_rpf_sum,
        "pooled_extent_frame0_ratio": resolution.pooled_extent_frame0_ratio,
        "pooled_start_rpf_sum": resolution.pooled_start_rpf_sum,
        "pooled_end_rpf_sum": resolution.pooled_end_rpf_sum,
        "pooled_extent_sample_count": resolution.pooled_extent_sample_count,
        "nested_competition_status": (
            "Pending"
            if status == "Reliable"
            and resolution.translated_extent_status in {"Supported", "Compatible"}
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
        "supporting_samples": ",".join(accumulator.supporting_samples),
        "best_sample": accumulator.best_sample,
        "best_sample_evidence": accumulator.best_level,
        "best_rpf_sum": (0.0 if features is None else features.rpf_sum),
        "best_rpf_per_codon": (0.0 if features is None else features.rpf_per_codon),
        "best_covered_codon": (0 if features is None else features.covered_codon),
        "best_coverage_ratio": (0.0 if features is None else features.coverage_ratio),
        "best_frame0_ratio": (0.0 if features is None else features.frame0_ratio),
        "best_supported_windows": (0 if features is None else features.supported_windows),
        "best_distributed_windows": (0 if features is None else features.distributed_windows),
        "best_signal_span": (0.0 if features is None else features.signal_span),
        "best_top_window_fraction": (0.0 if features is None else features.top_window_fraction),
        "best_localized_only": (False if features is None else features.localized_only),
        "representative_count": len(family.representatives),
        "representative_orf_ids": ",".join(item.orf_id for item in family.representatives),
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
        "selection_policy": ("canonical_ATG_then_longest_unless_replicated_RPF_override"),
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
            "\t".join(_format_value(row.get(column, "")) for column in columns) + "\n"
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
            accumulators = [_new_accumulator(family, task.config) for family in valid_families]
            canonical_indices = [_canonical_anchor_index(family) for family in valid_families]
            diagnostic_index_sets = []
            for family, canonical_index in zip(
                valid_families,
                canonical_indices,
            ):
                if canonical_index is None:
                    diagnostic_index_sets.append({0})
                    continue
                canonical_start = family.representatives[canonical_index].start_offset
                diagnostic_index_sets.append(
                    {canonical_index}
                    | {
                        index
                        for index, representative in enumerate(family.representatives)
                        if representative.start_offset < canonical_start
                        and representative.start_codon.upper() != "ATG"
                    }
                )

            for sample in samples:
                threshold = thresholds_by_sample[sample]
                for strand in ("+", "-"):
                    indices = [
                        index
                        for index, family in enumerate(valid_families)
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
                        codon_count = (family.common_body_end - family.common_body_start) // 3
                        level, _ = classify_sample(
                            features,
                            codon_count,
                            threshold,
                            task.config,
                        )
                        if LEVEL_RANK[level] >= LEVEL_RANK["LowConfidence"]:
                            accumulator.any_evidence = True
                        if LEVEL_RANK[level] > LEVEL_RANK[accumulator.best_level] or (
                            LEVEL_RANK[level] == LEVEL_RANK[accumulator.best_level]
                            and features.score > accumulator.best_score
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
                                accumulator.group_support_counts.get(group_id, 0) + 1
                            )
                            if level == "HighConfidence":
                                accumulator.high_count += 1
                                accumulator.group_high_counts[group_id] = (
                                    accumulator.group_high_counts.get(group_id, 0) + 1
                                )

                        canonical_index = canonical_indices[family_index]
                        diagnostic_indices = diagnostic_index_sets[family_index]
                        canonical_features = features
                        canonical_start: int | None = None
                        if canonical_index is not None:
                            canonical_start = family.representatives[canonical_index].start_offset
                            if canonical_start != family.common_body_start:
                                canonical_features = profile.region_features(
                                    canonical_start,
                                    family.common_body_end,
                                    threshold,
                                    task.config,
                                )

                        for candidate_index, representative in enumerate(family.representatives):
                            candidate_start = representative.start_offset
                            candidate_end = family.common_body_end
                            extension_end = family.common_body_start
                            full = (
                                features
                                if candidate_start == family.common_body_start
                                else profile.region_features(
                                    candidate_start,
                                    candidate_end,
                                    threshold,
                                    task.config,
                                )
                            )
                            if full.rpf_sum > 0:
                                accumulator.any_signal = True
                            extension = profile.region_features(
                                candidate_start,
                                extension_end,
                                threshold,
                                task.config,
                            )
                            full_codons = max(0, (candidate_end - candidate_start) // 3)
                            extension_codons = max(0, (extension_end - candidate_start) // 3)
                            full_supported, extension_supported = classify_candidate_extent(
                                full,
                                extension,
                                full_codons,
                                extension_codons,
                                threshold,
                                task.config,
                            )
                            if full_supported:
                                accumulator.candidate_full_support_counts[candidate_index] += 1
                            if extension_supported:
                                accumulator.candidate_extension_support_counts[candidate_index] += 1
                            if full_supported or (extension_supported and extension.rpf_sum > 0):
                                accumulator.candidate_support_samples[candidate_index].append(
                                    sample
                                )

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
                                if LEVEL_RANK[level] >= LEVEL_RANK[
                                    "MediumConfidence"
                                ] and _leading_region_is_negative(
                                    leading,
                                    features,
                                    leading_codons,
                                    threshold,
                                    task.config,
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
                                noncanonical_extension = profile.region_features(
                                    candidate_start,
                                    canonical_start,
                                    threshold,
                                    task.config,
                                )
                                noncanonical_codons = (canonical_start - candidate_start) // 3
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

                            previous_full = accumulator.candidate_best_full_features[
                                candidate_index
                            ]
                            if previous_full is None or full.score > previous_full.score:
                                accumulator.candidate_best_full_features[candidate_index] = full
                            previous_extension = accumulator.candidate_best_extension_features[
                                candidate_index
                            ]
                            if (
                                previous_extension is None
                                or extension.score > previous_extension.score
                            ):
                                accumulator.candidate_best_extension_features[candidate_index] = (
                                    extension
                                )
                            pooled = accumulator.candidate_pooled_frame_density[candidate_index]
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
                            for bin_index, (bin_start, bin_end) in enumerate(bin_ranges):
                                bin_features = profile.region_features(
                                    bin_start,
                                    bin_end,
                                    threshold,
                                    task.config,
                                )
                                _accumulate_pooled_segment(
                                    accumulator.candidate_bin_pooled_stats[candidate_index][
                                        bin_index
                                    ],
                                    bin_features,
                                )
                                if _bin_is_supported(
                                    bin_features,
                                    threshold,
                                ):
                                    accumulator.candidate_bin_support_counts[candidate_index][
                                        bin_index
                                    ] += 1
                                    sample_has_supported_bin = True
                            if sample_has_supported_bin:
                                accumulator.candidate_bin_support_samples[candidate_index].append(
                                    sample
                                )
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
                                    accumulator.candidate_boundary_pooled_stats[candidate_index][0],
                                    start_boundary,
                                )
                                _accumulate_pooled_segment(
                                    accumulator.candidate_boundary_pooled_stats[candidate_index][1],
                                    end_boundary,
                                )

            rows: list[Mapping[str, Any]] = []
            for family, accumulator in zip(valid_families, accumulators):
                row = _family_output_row(family, accumulator, task.group_names, task.config)
                rows.append(row)
                counts["total_families"] += 1
                if row["evidence_status"] == "Reliable":
                    counts["reliable_families"] += 1
                    if (
                        row["translated_extent_status"] in {"Supported", "Compatible"}
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
        query = f"SELECT orf_id, matched FROM reliable_export WHERE orf_id IN ({placeholders})"
        for row in connection.execute(query, chunk):
            matches[str(row["orf_id"])] = int(row["matched"])
    return matches


def _export_reliable_genepred(
    source_path: str | Path,
    temporary_path: Path,
    connection: sqlite3.Connection,
) -> int:
    """Filter scanner genePred records using the validated Reliable ID index."""
    expected = int(connection.execute("SELECT COUNT(*) FROM reliable_export").fetchone()[0])
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
                    duplicate_examples.append(f"{orf_id}@line{line_number}")
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

    with (
        _smart_open(source_path, "rt") as source,
        temporary_path.open(
            "w",
            encoding="utf-8",
            buffering=8 * 1024 * 1024,
        ) as target,
    ):
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
        connection.execute("SELECT COUNT(*) FROM reliable_export WHERE matched=0").fetchone()[0]
    )
    if missing_count:
        examples = [
            str(row[0])
            for row in connection.execute(
                "SELECT orf_id FROM reliable_export WHERE matched=0 ORDER BY orf_id LIMIT 10"
            )
        ]
        raise EvidenceEngineError(
            f"Scanner genePred is missing {missing_count:,} Reliable ORF(s). "
            "Examples: " + ", ".join(examples)
        )
    return expected
