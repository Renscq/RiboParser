#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Calibrate evidence thresholds with annotated ORFs.
# Input: Annotated controls and density tracks.
# Output: Length models and sample-specific thresholds.

"""Calibrate evidence thresholds with annotated ORFs."""

from __future__ import annotations

import json
import math
import sqlite3
from collections import defaultdict
from dataclasses import asdict, replace
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from utils.ribo.ArgsParser import progress_print

from .config import (
    ADVANCED_DEFAULTS,
    LONG_MODEL_CEILING,
    LONG_MODEL_FLOOR,
    MANUAL_THRESHOLD_DEFAULTS,
    SHORT_MODEL_CEILING,
    SHORT_MODEL_FLOOR,
    TARGET_LENGTH_CATEGORIES,
)
from .features import (
    _extract_sparse_signals,
    _manual_thresholds,
    segment_features,
)
from .input import (
    _build_family_geometry,
    _coding_nt_length,
    _load_cached_density,
)
from .models import (
    DensityTrack,
    EngineConfig,
    FamilyGeometry,
    Thresholds,
    _StageLogger,
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
        int(round(float(np.median(short_values)))) if short_values else config.short_max_codons
    )
    estimated_long = (
        int(round(float(np.median(long_values)))) if long_values else config.long_min_codons
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
    selected = [item for item in ordered if float(item["rpf_per_codon"]) >= cutoff]
    if len(selected) < config.positive_min_controls:
        selected = ordered[-config.positive_min_controls :]
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
        supported_windows = int(ADVANCED_DEFAULTS["min_supported_windows"])

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
                    float(ADVANCED_DEFAULTS["localized_top_window_fraction"]),
                    float(np.quantile(top_fractions, 1.0 - quantile)),
                ),
                float(ADVANCED_DEFAULTS["localized_top_window_fraction"]),
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
        thresholds = [_manual_thresholds(sample, config) for sample in samples]
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
        logger.write("Manual evidence mode: use the eight user-defined thresholds.")
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
                subset = [family for family in families if family.strand == strand]
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
                            "top_window_fraction": (features.top_window_fraction),
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
        progress_print(f"canonical controls: {sample_number}/{len(samples)}")

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

        moderate = float(np.clip(np.quantile(frame0, q), 0.40, 0.85))
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
