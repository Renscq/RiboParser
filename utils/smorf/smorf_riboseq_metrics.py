#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Calculate and classify smORF Ribo-seq evidence metrics.
# Input: Transcript-oriented nucleotide, codon, and downstream P-site profiles.
# Output: Quantification, periodicity, pause, release, shape, and evidence labels.

"""Metric functions for smORF Ribo-seq evidence analysis.

Metrics are divided into three groups:

1. Core evidence: coding-region abundance, codon coverage, and frame periodicity.
2. Supporting evidence: coverage shape and transcript-aware release.
3. Descriptive signals: start and stop pausing.

Metrics that are unreliable for short or low-count ORFs are marked ``NA``
instead of being forced to ``Weak``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

from .smorf_riboseq_constants import EvidenceThresholds


def safe_divide(
    numerator: float,
    denominator: float,
    default: float = 0.0,
) -> float:
    """Safely divide finite values.

    Args:
        numerator: Numerator.
        denominator: Denominator.
        default: Value returned for invalid or zero denominators.

    Returns:
        Division result.
    """
    if (
        not math.isfinite(float(denominator))
        or float(denominator) == 0.0
    ):
        return float(default)
    return float(numerator) / float(denominator)


def gini_index(values: np.ndarray) -> float:
    """Calculate the Gini index for non-negative values.

    Args:
        values: Input vector.

    Returns:
        Gini index in [0, 1].
    """
    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return 0.0
    if np.any(array < 0):
        raise ValueError("Gini input must be non-negative.")
    total = float(array.sum())
    if total <= 0:
        return 0.0

    array = np.sort(array)
    count = array.size
    index = np.arange(1, count + 1, dtype=np.float64)
    value = (
        2.0 * np.sum(index * array) / (count * total)
        - (count + 1.0) / count
    )
    return float(max(0.0, min(1.0, value)))


def window_mean(
    values: np.ndarray,
    start: int,
    end: int,
) -> float:
    """Calculate a clipped interval mean.

    Args:
        values: Input vector.
        start: Inclusive start.
        end: Exclusive end.

    Returns:
        Mean value or zero for an empty interval.
    """
    length = len(values)
    clipped_start = max(0, int(start))
    clipped_end = min(length, int(end))
    if clipped_end <= clipped_start:
        return 0.0
    return float(np.mean(values[clipped_start:clipped_end]))


def quantify_profile(
    nt_profile: np.ndarray,
    codon_profile: np.ndarray,
    coding_nt_length: int | None = None,
) -> dict[str, float | int]:
    """Calculate coding-region abundance and coverage.

    ``rpf_sum`` is calculated from the coding region and excludes a terminal
    stop codon when ``coding_nt_length`` is supplied.

    Args:
        nt_profile: ORF nucleotide P-site profile.
        codon_profile: Coding codon profile.
        coding_nt_length: Optional coding length excluding stop codon.

    Returns:
        Quantification metrics.
    """
    usable_length = (
        len(nt_profile)
        if coding_nt_length is None
        else max(0, min(int(coding_nt_length), len(nt_profile)))
    )
    coding_nt = np.asarray(
        nt_profile[:usable_length],
        dtype=np.float64,
    )
    codon_values = np.asarray(codon_profile, dtype=np.float64)

    rpf_sum = float(coding_nt.sum())
    nt_count = int(coding_nt.size)
    codon_count = int(codon_values.size)
    covered_nt = int(np.count_nonzero(coding_nt > 0))
    covered_codon = int(np.count_nonzero(codon_values > 0))

    return {
        "rpf_sum": rpf_sum,
        "rpf_mean": safe_divide(rpf_sum, nt_count),
        "rpf_per_codon": safe_divide(rpf_sum, codon_count),
        "covered_nt": covered_nt,
        "coverage_ratio": safe_divide(covered_nt, nt_count),
        "covered_codon": covered_codon,
        "covered_codon_ratio": safe_divide(
            covered_codon,
            codon_count,
        ),
        "max_density": (
            float(coding_nt.max()) if nt_count else 0.0
        ),
    }


def calculate_periodicity(
    nt_profile: np.ndarray,
    coding_nt_length: int,
    thresholds: EvidenceThresholds,
) -> dict[str, float | str | bool]:
    """Calculate frame-specific P-site periodicity.

    Periodicity is evaluable only after the coding region satisfies the same
    minimum abundance and covered-codon requirements used by the core evidence
    gate.

    Args:
        nt_profile: Transcript-oriented nucleotide profile.
        coding_nt_length: Coding length excluding stop codon.
        thresholds: Evidence thresholds.

    Returns:
        Frame densities, ratios, normalized score, label, and evaluability.
    """
    usable_length = max(
        0,
        min(int(coding_nt_length), len(nt_profile)),
    )
    usable_length -= usable_length % 3

    empty_result = {
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
    }
    if usable_length <= 0:
        return empty_result

    coding = np.asarray(
        nt_profile[:usable_length],
        dtype=np.float64,
    )
    frame0 = float(coding[0::3].sum())
    frame1 = float(coding[1::3].sum())
    frame2 = float(coding[2::3].sum())
    total = frame0 + frame1 + frame2

    frame0_ratio = safe_divide(frame0, total)
    frame1_ratio = safe_divide(frame1, total)
    frame2_ratio = safe_divide(frame2, total)
    alt_mean = (frame1 + frame2) / 2.0
    frame0_vs_alt = safe_divide(
        frame0 + thresholds.pseudocount,
        alt_mean + thresholds.pseudocount,
    )
    periodicity_score = max(
        0.0,
        min(
            1.0,
            safe_divide(
                frame0_ratio - (1.0 / 3.0),
                2.0 / 3.0,
            ),
        ),
    )

    codon_profile = coding.reshape(-1, 3).sum(axis=1)
    covered_codon = int(np.count_nonzero(codon_profile > 0))
    evaluable = (
        total >= thresholds.min_rpf_sum
        and covered_codon >= thresholds.min_covered_codon
    )

    if not evaluable:
        label = "NA"
    elif frame0_ratio >= thresholds.strong_periodicity:
        label = "Strong"
    elif frame0_ratio >= thresholds.moderate_periodicity:
        label = "Moderate"
    else:
        label = "Weak"

    return {
        "frame0_density": frame0,
        "frame1_density": frame1,
        "frame2_density": frame2,
        "frame0_ratio": frame0_ratio,
        "frame1_ratio": frame1_ratio,
        "frame2_ratio": frame2_ratio,
        "frame0_vs_alt_ratio": frame0_vs_alt,
        "periodicity_score": periodicity_score,
        "periodicity_label": label,
        "periodicity_evaluable": evaluable,
    }


def calculate_pausing(
    codon_profile: np.ndarray,
    thresholds: EvidenceThresholds,
    pseudocount: float | None = None,
) -> dict[str, float | str | bool]:
    """Calculate descriptive start and pre-stop pausing ratios.

    Pausing is not used as a mandatory translation criterion because it is
    sensitive to library protocol and ORF length.

    Args:
        codon_profile: Coding codon profile.
        thresholds: Evidence thresholds.
        pseudocount: Optional ratio pseudocount.

    Returns:
        Pause means, ratios, labels, and evaluability.
    """
    values = np.asarray(codon_profile, dtype=np.float64)
    count = len(values)
    pseudo = (
        thresholds.pseudocount
        if pseudocount is None
        else float(pseudocount)
    )

    if count < thresholds.min_boundary_codons:
        return {
            "start_codon_mean": 0.0,
            "body_codon_mean": 0.0,
            "pre_stop_codon_mean": 0.0,
            "start_pause_ratio": 0.0,
            "stop_pause_ratio": 0.0,
            "start_pause_label": "NA",
            "stop_pause_label": "NA",
            "pausing_label": "NA",
            "pausing_evaluable": False,
        }

    start_mean = window_mean(values, 0, 3)
    body_mean = window_mean(values, 3, count - 3)
    pre_stop_mean = window_mean(values, count - 3, count)

    start_ratio = safe_divide(
        start_mean + pseudo,
        body_mean + pseudo,
    )
    stop_ratio = safe_divide(
        pre_stop_mean + pseudo,
        body_mean + pseudo,
    )

    def classify(
        ratio: float,
        moderate: float,
        strong: float,
    ) -> str:
        """Return a pause label for one ratio."""
        if ratio >= strong:
            return "Strong"
        if ratio >= moderate:
            return "Moderate"
        return "Weak"

    start_label = classify(
        start_ratio,
        thresholds.moderate_start_pause,
        thresholds.strong_start_pause,
    )
    stop_label = classify(
        stop_ratio,
        thresholds.moderate_stop_pause,
        thresholds.strong_stop_pause,
    )
    combined_label = (
        "Strong"
        if "Strong" in {start_label, stop_label}
        else (
            "Moderate"
            if "Moderate" in {start_label, stop_label}
            else "Weak"
        )
    )

    return {
        "start_codon_mean": start_mean,
        "body_codon_mean": body_mean,
        "pre_stop_codon_mean": pre_stop_mean,
        "start_pause_ratio": start_ratio,
        "stop_pause_ratio": stop_ratio,
        "start_pause_label": start_label,
        "stop_pause_label": stop_label,
        "pausing_label": combined_label,
        "pausing_evaluable": True,
    }


def calculate_release(
    codon_profile: np.ndarray,
    downstream_nt_profile: np.ndarray,
    post_stop_codons: int,
    thresholds: EvidenceThresholds,
    pseudocount: float | None = None,
    context: str = "transcript",
) -> dict[str, float | str | bool]:
    """Calculate transcript-aware release after the stop codon.

    Release is evaluable only when a transcript-spliced downstream profile is
    available. Genomic downstream sequence is intentionally not used because
    internal ORFs may be followed by introns or downstream exons.

    Args:
        codon_profile: Coding codon profile.
        downstream_nt_profile: Transcript-downstream nucleotide profile.
        post_stop_codons: Requested downstream codons.
        thresholds: Evidence thresholds.
        pseudocount: Optional ratio pseudocount.
        context: Downstream context label.

    Returns:
        Release metrics and evaluability.
    """
    pseudo = (
        thresholds.pseudocount
        if pseudocount is None
        else float(pseudocount)
    )
    coding_values = np.asarray(codon_profile, dtype=np.float64)
    downstream_values = np.asarray(
        downstream_nt_profile,
        dtype=np.float64,
    )

    requested_nt = max(0, int(post_stop_codons) * 3)
    usable_nt = min(len(downstream_values), requested_nt)
    usable_nt -= usable_nt % 3
    evaluable = (
        context == "transcript"
        and len(coding_values) >= 3
        and usable_nt >= 3
    )

    if not evaluable:
        return {
            "post_stop_mean": 0.0,
            "release_ratio": 0.0,
            "release_drop_score": 0.0,
            "release_label": "NA",
            "release_evaluable": False,
            "release_context": context,
            "post_stop_nt": int(usable_nt),
        }

    pre_stop_mean = float(np.mean(coding_values[-3:]))
    downstream_codon = (
        downstream_values[:usable_nt]
        .reshape(-1, 3)
        .sum(axis=1)
    )
    post_stop_mean = float(np.mean(downstream_codon))
    release_ratio = safe_divide(
        pre_stop_mean + pseudo,
        post_stop_mean + pseudo,
    )
    drop_score = max(
        0.0,
        min(
            1.0,
            1.0 - safe_divide(
                post_stop_mean,
                pre_stop_mean,
                default=1.0,
            ),
        ),
    )

    if pre_stop_mean <= 0:
        label = "Weak"
    elif release_ratio >= thresholds.strong_release:
        label = "Strong"
    elif release_ratio >= thresholds.moderate_release:
        label = "Moderate"
    else:
        label = "Weak"

    return {
        "post_stop_mean": post_stop_mean,
        "release_ratio": release_ratio,
        "release_drop_score": drop_score,
        "release_label": label,
        "release_evaluable": True,
        "release_context": context,
        "post_stop_nt": int(usable_nt),
    }


def classify_coverage_shape(
    codon_profile: np.ndarray,
    thresholds: EvidenceThresholds,
) -> dict[str, float | str | bool]:
    """Classify codon-level coverage shape.

    Shape is not evaluated for very short ORFs or zero-density profiles because
    Gini and top-fraction statistics are unstable in those settings.

    Args:
        codon_profile: Coding codon profile.
        thresholds: Evidence thresholds.

    Returns:
        Shape metrics, label, and evaluability.
    """
    values = np.asarray(codon_profile, dtype=np.float64)
    count = len(values)
    total = float(values.sum())

    if count < thresholds.min_shape_codons or total <= 0:
        return {
            "codon_gini": 0.0,
            "max_to_mean_ratio": 0.0,
            "top10_fraction": 0.0,
            "coverage_shape": "NA",
            "shape_evaluable": False,
        }

    coverage_ratio = safe_divide(
        int(np.count_nonzero(values > 0)),
        count,
    )
    mean_value = float(values.mean())
    max_value = float(values.max())
    max_to_mean = safe_divide(max_value, mean_value)

    top_count = max(1, int(math.ceil(count * 0.10)))
    sorted_values = np.sort(values)
    top_fraction = safe_divide(
        float(sorted_values[-top_count:].sum()),
        total,
    )
    gini = gini_index(values)

    if (
        coverage_ratio >= thresholds.uniform_coverage_ratio
        and gini <= thresholds.uniform_gini
        and max_to_mean <= thresholds.uniform_max_to_mean
    ):
        label = "Uniform"
    elif (
        max_to_mean >= thresholds.skewed_max_to_mean
        or top_fraction >= thresholds.skewed_top_fraction
    ):
        label = "Skewed"
    else:
        label = "Intermediate"

    return {
        "codon_gini": gini,
        "max_to_mean_ratio": max_to_mean,
        "top10_fraction": top_fraction,
        "coverage_shape": label,
        "shape_evaluable": True,
    }


def classify_translation_evidence(
    quant: Mapping[str, float | int],
    periodicity: Mapping[str, float | str | bool],
    release: Mapping[str, float | str | bool],
    coverage_shape: Mapping[str, float | str | bool],
    thresholds: EvidenceThresholds,
) -> tuple[str, str, str]:
    """Assign a transparent sample-level translation-evidence label.

    The core gate requires abundance, covered codons, and codon coverage.
    Periodicity is the primary translation signal. Uniform coverage or
    transcript-aware release can promote a strong-periodicity ORF to high
    confidence. Release is optional when transcript context is unavailable.

    Args:
        quant: Quantification metrics.
        periodicity: Periodicity metrics.
        release: Release metrics.
        coverage_shape: Coverage-shape metrics.
        thresholds: Evidence thresholds.

    Returns:
        Evidence label, semicolon-delimited evidence components, and
        semicolon-delimited failure reasons.
    """
    rpf_sum = float(quant["rpf_sum"])
    covered_codon = int(quant["covered_codon"])
    codon_coverage = float(quant["covered_codon_ratio"])

    if rpf_sum <= 0:
        return "NoEvidence", "none", "zero_rpf"

    failures = []
    if rpf_sum < thresholds.min_rpf_sum:
        failures.append("low_rpf_sum")
    if covered_codon < thresholds.min_covered_codon:
        failures.append("few_covered_codons")
    if codon_coverage < thresholds.min_codon_coverage:
        failures.append("low_codon_coverage")
    if failures:
        return "LowConfidence", "abundance_only", ";".join(failures)

    periodicity_label = str(periodicity["periodicity_label"])
    shape_label = str(coverage_shape["coverage_shape"])
    release_label = str(release["release_label"])

    components = [
        "abundance",
        "codon_coverage",
    ]
    if periodicity_label in {"Moderate", "Strong"}:
        components.append(
            f"periodicity_{periodicity_label.lower()}"
        )
    if shape_label == "Uniform":
        components.append("uniform_coverage")
    elif shape_label == "Skewed":
        components.append("skewed_coverage")
    if release_label in {"Moderate", "Strong"}:
        components.append(
            f"release_{release_label.lower()}"
        )

    if (
        periodicity_label == "Strong"
        and shape_label != "Skewed"
        and (
            shape_label == "Uniform"
            or release_label in {"Moderate", "Strong"}
        )
    ):
        label = "HighConfidence"
        reason = "PASS"
    elif (
        periodicity_label in {"Moderate", "Strong"}
        and shape_label != "Skewed"
    ):
        label = "MediumConfidence"
        reason = "PASS"
    else:
        label = "LowConfidence"
        low_reasons = []
        if periodicity_label in {"Weak", "NA"}:
            low_reasons.append("weak_or_unevaluable_periodicity")
        if shape_label == "Skewed":
            low_reasons.append("skewed_coverage")
        reason = ";".join(low_reasons) or "insufficient_support"

    return label, ";".join(components), reason
