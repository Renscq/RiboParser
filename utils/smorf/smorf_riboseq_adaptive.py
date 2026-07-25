#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.006
# Function: Calculate length-adaptive and sliding-window Ribo-seq evidence.
# Input: Transcript-oriented nucleotide and codon P-site profiles.
# Output: Adaptive ORF evidence metrics and evidence classifications.

"""Length-adaptive Ribo-seq evidence metrics for smORF families."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from statistics import median
from typing import Final, Iterable, Mapping, Sequence

import numpy as np

EVIDENCE_LEVELS: Final[dict[str, int]] = {
    "NoEvidence": 0,
    "LowConfidence": 1,
    "MediumConfidence": 2,
    "HighConfidence": 3,
}


@dataclass(frozen=True, slots=True)
class AdaptiveEvidenceConfig:
    """Store length-adaptive evidence parameters."""

    min_rpf_sum: float = 5.0
    min_covered_codon: int = 3
    min_codon_coverage: float = 0.15
    moderate_periodicity: float = 0.55
    strong_periodicity: float = 0.70
    short_max_codons: int = 30
    long_min_codons: int = 100
    window_codons: int = 20
    window_step_codons: int = 5
    min_window_rpf: float = 3.0
    min_window_covered_codons: int = 3
    min_supported_windows: int = 2
    min_window_gap_codons: int = 15
    min_coverage_span: float = 0.35
    localized_span_max: float = 0.20
    localized_top_window_fraction: float = 0.70
    start_window_codons: int = 5
    body_margin_codons: int = 5
    pseudocount: float = 0.10
    positive_quantile: float = 0.20
    positive_min_controls: int = 30
    positive_min_rpf_sum: float = 10.0

    def validate(self) -> None:
        """Validate parameter ranges and ordering."""
        if self.min_rpf_sum < 0 or self.min_window_rpf < 0:
            raise ValueError("RPF thresholds must be >= 0.")
        if self.min_covered_codon < 1:
            raise ValueError("min_covered_codon must be >= 1.")
        if self.min_window_covered_codons < 1:
            raise ValueError("min_window_covered_codons must be >= 1.")
        if not 0 <= self.min_codon_coverage <= 1:
            raise ValueError("min_codon_coverage must be in [0, 1].")
        if not (
            0 <= self.moderate_periodicity
            <= self.strong_periodicity
            <= 1
        ):
            raise ValueError(
                "Require moderate_periodicity <= strong_periodicity in [0, 1]."
            )
        if self.short_max_codons < 3:
            raise ValueError("short_max_codons must be >= 3.")
        if self.long_min_codons <= self.short_max_codons:
            raise ValueError("long_min_codons must exceed short_max_codons.")
        if self.window_codons < 6:
            raise ValueError("window_codons must be >= 6.")
        if self.window_step_codons < 1:
            raise ValueError("window_step_codons must be >= 1.")
        if self.min_supported_windows < 1:
            raise ValueError("min_supported_windows must be >= 1.")
        if self.min_window_gap_codons < 0:
            raise ValueError("min_window_gap_codons must be >= 0.")
        for name in (
            "min_coverage_span",
            "localized_span_max",
            "localized_top_window_fraction",
            "positive_quantile",
        ):
            value = float(getattr(self, name))
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be in [0, 1].")
        if self.start_window_codons < 1:
            raise ValueError("start_window_codons must be >= 1.")
        if self.body_margin_codons < 0:
            raise ValueError("body_margin_codons must be >= 0.")
        if self.pseudocount < 0:
            raise ValueError("pseudocount must be >= 0.")
        if self.positive_min_controls < 1:
            raise ValueError("positive_min_controls must be >= 1.")

    def with_calibration(
        self,
        controls: Sequence[Mapping[str, float | int | str]],
    ) -> tuple["AdaptiveEvidenceConfig", dict[str, float | int | str]]:
        """Calibrate permissive lower bounds from annotated mORF controls.

        Calibration can relax a preset for a low-quality library, but never
        makes the preset more stringent. Very poor controls are prevented from
        lowering frame periodicity below 0.45.
        """
        eligible = [
            control
            for control in controls
            if float(control.get("rpf_sum", 0.0))
            >= self.positive_min_rpf_sum
            and bool(control.get("periodicity_evaluable", False))
        ]
        summary: dict[str, float | int | str] = {
            "positive_control_count": len(controls),
            "positive_control_eligible": len(eligible),
            "calibration_status": "insufficient_controls",
            "base_moderate_periodicity": self.moderate_periodicity,
            "base_strong_periodicity": self.strong_periodicity,
            "calibrated_moderate_periodicity": self.moderate_periodicity,
            "calibrated_strong_periodicity": self.strong_periodicity,
            "calibrated_min_window_rpf": self.min_window_rpf,
        }
        if len(eligible) < self.positive_min_controls:
            return self, summary

        frame_ratios = np.asarray(
            [float(item["frame0_ratio"]) for item in eligible],
            dtype=np.float64,
        )
        rpf_per_codon = np.asarray(
            [float(item.get("rpf_per_codon", 0.0)) for item in eligible],
            dtype=np.float64,
        )
        q = self.positive_quantile
        lower_frame = float(np.quantile(frame_ratios, q))
        median_frame = float(np.quantile(frame_ratios, 0.50))
        lower_rpf = float(np.quantile(rpf_per_codon, q))

        moderate = max(0.45, min(self.moderate_periodicity, lower_frame))
        strong = max(
            moderate + 0.08,
            min(self.strong_periodicity, median_frame),
        )
        strong = min(strong, 0.90)
        window_rpf = max(
            2.0,
            min(
                self.min_window_rpf,
                lower_rpf * min(self.window_codons, 20),
            ),
        )
        calibrated = replace(
            self,
            moderate_periodicity=moderate,
            strong_periodicity=strong,
            min_window_rpf=window_rpf,
        )
        calibrated.validate()
        summary.update(
            {
                "calibration_status": "calibrated",
                "positive_frame0_q20": lower_frame,
                "positive_frame0_median": median_frame,
                "positive_rpf_per_codon_q20": lower_rpf,
                "calibrated_moderate_periodicity": moderate,
                "calibrated_strong_periodicity": strong,
                "calibrated_min_window_rpf": window_rpf,
            }
        )
        return calibrated, summary


@dataclass(frozen=True, slots=True)
class WindowEvidence:
    """Store one sliding-window evidence record."""

    start_codon: int
    end_codon: int
    rpf_sum: float
    covered_codon: int
    frame0_ratio: float
    score: float
    supported: bool


@dataclass(frozen=True, slots=True)
class AdaptiveEvidenceResult:
    """Store adaptive evidence metrics for one ORF/profile."""

    metrics: dict[str, float | int | str | bool]
    windows: tuple[WindowEvidence, ...]


def _safe_divide(
    numerator: float,
    denominator: float,
    default: float = 0.0,
) -> float:
    """Return a finite division result."""
    if not math.isfinite(float(denominator)) or float(denominator) == 0.0:
        return float(default)
    return float(numerator) / float(denominator)


def _frame_metrics(nt_values: np.ndarray) -> tuple[float, float, float, float]:
    """Return frame densities and frame-0 ratio."""
    usable = len(nt_values) - (len(nt_values) % 3)
    if usable <= 0:
        return 0.0, 0.0, 0.0, 0.0
    coding = np.asarray(nt_values[:usable], dtype=np.float64)
    frame0 = float(coding[0::3].sum())
    frame1 = float(coding[1::3].sum())
    frame2 = float(coding[2::3].sum())
    total = frame0 + frame1 + frame2
    return frame0, frame1, frame2, _safe_divide(frame0, total)


def _window_starts(codon_count: int, config: AdaptiveEvidenceConfig) -> list[int]:
    """Return deterministic sliding-window start coordinates."""
    if codon_count <= config.window_codons:
        return [0] if codon_count > 0 else []
    starts = list(
        range(
            0,
            codon_count - config.window_codons + 1,
            config.window_step_codons,
        )
    )
    final_start = codon_count - config.window_codons
    if starts[-1] != final_start:
        starts.append(final_start)
    return starts


def calculate_window_evidence(
    nt_profile: np.ndarray,
    codon_profile: np.ndarray,
    config: AdaptiveEvidenceConfig,
) -> tuple[WindowEvidence, ...]:
    """Calculate sliding-window abundance, coverage, and periodicity."""
    codon_values = np.asarray(codon_profile, dtype=np.float64)
    windows: list[WindowEvidence] = []
    for start in _window_starts(len(codon_values), config):
        end = min(len(codon_values), start + config.window_codons)
        codon_slice = codon_values[start:end]
        nt_slice = np.asarray(nt_profile[start * 3 : end * 3], dtype=np.float64)
        rpf_sum = float(nt_slice.sum())
        covered = int(np.count_nonzero(codon_slice > 0))
        _frame0, _frame1, _frame2, frame0_ratio = _frame_metrics(nt_slice)
        abundance_score = min(1.0, _safe_divide(rpf_sum, config.min_window_rpf))
        coverage_score = min(
            1.0,
            _safe_divide(covered, config.min_window_covered_codons),
        )
        periodicity_score = max(
            0.0,
            min(
                1.0,
                _safe_divide(frame0_ratio - (1.0 / 3.0), 2.0 / 3.0),
            ),
        )
        score = (
            0.35 * abundance_score
            + 0.25 * coverage_score
            + 0.40 * periodicity_score
        )
        supported = (
            rpf_sum >= config.min_window_rpf
            and covered >= min(
                config.min_window_covered_codons,
                max(1, len(codon_slice)),
            )
            and frame0_ratio >= config.moderate_periodicity
        )
        windows.append(
            WindowEvidence(
                start_codon=start,
                end_codon=end,
                rpf_sum=rpf_sum,
                covered_codon=covered,
                frame0_ratio=frame0_ratio,
                score=score,
                supported=supported,
            )
        )
    return tuple(windows)


def _supported_span(
    supported_windows: Sequence[WindowEvidence],
    codon_count: int,
) -> tuple[float, int, int, int]:
    """Return supported span ratio and spatially separated window count."""
    if not supported_windows or codon_count <= 0:
        return 0.0, -1, -1, 0
    first = min(window.start_codon for window in supported_windows)
    last = max(window.end_codon for window in supported_windows)
    separated = 1
    previous = min(supported_windows, key=lambda item: item.start_codon)
    for window in sorted(supported_windows, key=lambda item: item.start_codon)[1:]:
        if window.start_codon - previous.start_codon >= 1:
            separated += 1
            previous = window
    return min(1.0, (last - first) / codon_count), first, last, separated


def _max_separated_windows(
    windows: Sequence[WindowEvidence],
    gap_codons: int,
) -> int:
    """Count greedily selected supported windows separated by a minimum gap."""
    selected = 0
    last_start: int | None = None
    for window in sorted(windows, key=lambda item: item.start_codon):
        if not window.supported:
            continue
        if last_start is None or window.start_codon - last_start >= gap_codons:
            selected += 1
            last_start = window.start_codon
    return selected


def _length_scaled_min_rpf(codon_count: int, config: AdaptiveEvidenceConfig) -> float:
    """Return a sublinear abundance threshold for long ORFs."""
    scale = math.sqrt(max(1.0, codon_count / max(1, config.short_max_codons)))
    return config.min_rpf_sum * min(4.0, scale)


def evaluate_adaptive_evidence(
    nt_profile: np.ndarray,
    codon_profile: np.ndarray,
    *,
    category: str,
    config: AdaptiveEvidenceConfig,
    release_label: str = "NA",
    release_evaluable: bool = False,
) -> AdaptiveEvidenceResult:
    """Classify one ORF with length-adaptive and distributed evidence rules."""
    config.validate()
    nt_values = np.asarray(nt_profile, dtype=np.float64)
    codon_values = np.asarray(codon_profile, dtype=np.float64)
    codon_count = int(len(codon_values))
    rpf_sum = float(nt_values.sum())
    covered_codon = int(np.count_nonzero(codon_values > 0))
    covered_ratio = _safe_divide(covered_codon, codon_count)
    frame0, frame1, frame2, frame0_ratio = _frame_metrics(nt_values)
    periodicity_evaluable = (
        rpf_sum >= config.min_rpf_sum
        and covered_codon >= config.min_covered_codon
    )
    windows = calculate_window_evidence(nt_values, codon_values, config)
    supported = tuple(window for window in windows if window.supported)
    span_ratio, first_window, last_window, _raw_separated = _supported_span(
        supported,
        codon_count,
    )
    separated_count = _max_separated_windows(
        supported,
        config.min_window_gap_codons,
    )
    top_window_rpf = max((window.rpf_sum for window in windows), default=0.0)
    top_window_fraction = _safe_divide(top_window_rpf, rpf_sum)
    covered_positions = np.flatnonzero(codon_values > 0)
    if covered_positions.size:
        signal_span_ratio = min(
            1.0,
            (int(covered_positions[-1]) - int(covered_positions[0]) + 1)
            / max(1, codon_count),
        )
        signal_bins = {
            min(2, int(position * 3 / max(1, codon_count)))
            for position in covered_positions
        }
        signal_bin_count = len(signal_bins)
    else:
        signal_span_ratio = 0.0
        signal_bin_count = 0
    top_scores = sorted((window.score for window in windows), reverse=True)[:3]
    top3_median = float(median(top_scores)) if top_scores else 0.0

    start_end = min(codon_count, config.start_window_codons)
    start_nt = nt_values[: start_end * 3]
    start_rpf = float(start_nt.sum())
    _s0, _s1, _s2, start_frame0_ratio = _frame_metrics(start_nt)
    start_support = (
        start_rpf >= max(2.0, config.min_window_rpf * 0.5)
        and start_frame0_ratio >= config.moderate_periodicity
    )

    body_start = min(codon_count, config.body_margin_codons)
    body_end = max(body_start, codon_count - config.body_margin_codons)
    body_supported = any(
        window.supported
        and window.end_codon > body_start
        and window.start_codon < body_end
        for window in windows
    )
    release_support = release_evaluable and release_label in {"Moderate", "Strong"}

    global_core = (
        rpf_sum >= config.min_rpf_sum
        and covered_codon >= config.min_covered_codon
        and covered_ratio >= config.min_codon_coverage
        and frame0_ratio >= config.moderate_periodicity
    )
    global_strong = (
        global_core
        and frame0_ratio >= config.strong_periodicity
        and rpf_sum >= config.min_rpf_sum * 2.0
    )
    distributed_support = (
        separated_count >= config.min_supported_windows
        and span_ratio >= config.min_coverage_span
    )
    boundary_support = start_support and body_supported and release_support
    distributed_sparse_support = (
        codon_count >= config.long_min_codons
        and rpf_sum >= config.min_rpf_sum
        and covered_codon >= config.min_covered_codon
        and frame0_ratio >= config.strong_periodicity
        and signal_span_ratio >= config.min_coverage_span
        and signal_bin_count >= 2
    )
    high_depth_support = (
        rpf_sum >= _length_scaled_min_rpf(codon_count, config)
        and frame0_ratio >= config.moderate_periodicity
        and max(span_ratio, signal_span_ratio) >= config.min_coverage_span
    )
    localized_only = (
        codon_count >= config.long_min_codons
        and (bool(supported) or covered_codon > 0)
        and (
            signal_span_ratio < config.localized_span_max
            or top_window_fraction >= config.localized_top_window_fraction
        )
        and not distributed_sparse_support
    )

    category_key = str(category).strip().lower()
    is_lnc = category_key == "lncorf"
    length_class: str
    if codon_count <= config.short_max_codons:
        length_class = "short"
        if global_strong:
            evidence = "HighConfidence"
            reason = "short_global_strong"
        elif global_core:
            evidence = "MediumConfidence"
            reason = "short_global_core"
        elif supported and frame0_ratio >= config.moderate_periodicity:
            evidence = "LowConfidence"
            reason = "short_local_support"
        else:
            evidence = "NoEvidence"
            reason = "short_insufficient"
    elif codon_count < config.long_min_codons:
        length_class = "medium"
        if global_strong and (distributed_support or len(supported) >= 2):
            evidence = "HighConfidence"
            reason = "medium_global_distributed"
        elif global_core or distributed_support or boundary_support:
            evidence = "MediumConfidence"
            reason = "medium_adaptive_support"
        elif supported:
            evidence = "LowConfidence"
            reason = "medium_local_support"
        else:
            evidence = "NoEvidence"
            reason = "medium_insufficient"
    else:
        length_class = "long"
        long_whole_support = (
            distributed_support
            or distributed_sparse_support
            or boundary_support
            or high_depth_support
        )
        if localized_only:
            evidence = "NoEvidence" if is_lnc else "LowConfidence"
            reason = (
                "lncORF_localized_only_not_whole_orf"
                if is_lnc
                else "long_localized_only"
            )
        elif is_lnc and not (
            distributed_support
            or distributed_sparse_support
            or boundary_support
        ):
            evidence = "LowConfidence" if supported else "NoEvidence"
            reason = (
                "lncORF_requires_distributed_or_boundary_support"
                if supported
                else "long_insufficient"
            )
        elif global_strong and long_whole_support:
            evidence = "HighConfidence"
            reason = "long_global_distributed"
        elif long_whole_support:
            evidence = "MediumConfidence"
            reason = "long_adaptive_support"
        elif supported:
            evidence = "LowConfidence"
            reason = "long_local_support"
        else:
            evidence = "NoEvidence"
            reason = "long_insufficient"

    if rpf_sum <= 0:
        evaluability = "NoSignal"
    elif not periodicity_evaluable and not supported:
        evaluability = "InsufficientDepth"
    elif localized_only:
        evaluability = "LocalizedOnly"
    elif evidence == "NoEvidence":
        evaluability = "EvaluableNotSupported"
    else:
        evaluability = "Translated"

    metrics: dict[str, float | int | str | bool] = {
        "coding_codon_count": codon_count,
        "rpf_sum": rpf_sum,
        "rpf_per_codon": _safe_divide(rpf_sum, codon_count),
        "covered_codon": covered_codon,
        "covered_codon_ratio": covered_ratio,
        "frame0_density": frame0,
        "frame1_density": frame1,
        "frame2_density": frame2,
        "frame0_ratio": frame0_ratio,
        "periodicity_evaluable": periodicity_evaluable,
        "length_class": length_class,
        "window_count": len(windows),
        "supported_window_count": len(supported),
        "separated_supported_window_count": separated_count,
        "supported_window_fraction": _safe_divide(len(supported), len(windows)),
        "best_window_rpf": top_window_rpf,
        "best_window_frame0_ratio": max(
            (window.frame0_ratio for window in windows),
            default=0.0,
        ),
        "top3_window_score_median": top3_median,
        "coverage_span_ratio": max(span_ratio, signal_span_ratio),
        "window_coverage_span_ratio": span_ratio,
        "signal_span_ratio": signal_span_ratio,
        "signal_bin_count": signal_bin_count,
        "first_supported_window": first_window,
        "last_supported_window": last_window,
        "top_window_rpf_fraction": top_window_fraction,
        "start_rpf": start_rpf,
        "start_frame0_ratio": start_frame0_ratio,
        "start_support": start_support,
        "body_support": body_supported,
        "release_support": release_support,
        "global_core_support": global_core,
        "distributed_support": distributed_support,
        "distributed_sparse_support": distributed_sparse_support,
        "boundary_support": boundary_support,
        "high_depth_support": high_depth_support,
        "localized_only": localized_only,
        "evidence_evaluability": evaluability,
        "evidence_reason": reason,
        "translation_evidence": evidence,
    }
    return AdaptiveEvidenceResult(metrics=metrics, windows=windows)


def evidence_max(labels: Iterable[str]) -> str:
    """Return the strongest evidence label in an iterable."""
    return max(
        labels,
        key=lambda label: EVIDENCE_LEVELS.get(str(label), -1),
        default="NoEvidence",
    )
