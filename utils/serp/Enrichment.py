#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev001
# Function: Calculate matched-replicate local SeRP enrichment profiles.
# Input: Per-codon normalized control and IP density profiles.
# Output: Replicate-specific ratios and median consensus enrichment profiles.

"""Shared enrichment calculations for selective ribosome profiling analyses."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd


def calculate_matched_local_enrichment(
    raw_gene_rpm: pd.DataFrame,
    control_samples: Sequence[str],
    ip_samples: Sequence[str],
    window: int = 5,
    pseudocount: float = 0.1,
    threshold: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate matched-replicate local IP/control enrichment.

    Parameters
    ----------
    raw_gene_rpm : pandas.DataFrame
        Per-position normalized densities. The index defines positional order.
    control_samples : Sequence[str]
        Ordered control sample names.
    ip_samples : Sequence[str]
        Ordered IP sample names. Each IP sample is paired with the control
        sample at the same sequence position.
    window : int, default=5
        Centered rolling-sum window in codons. Use 1 for unsmoothed codon-level
        enrichment.
    pseudocount : float, default=0.1
        Positive pseudocount added to normalized density before ratio
        calculation.
    threshold : float, optional
        Enrichment threshold used only to calculate pointwise replicate support.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Replicate-specific enrichment profiles and a consensus profile. The
        consensus ``enrich`` column is the median across matched replicate
        ratios. A ``support`` column is included when ``threshold`` is supplied.

    Notes
    -----
    The ratio is calculated from rolling sums rather than from a rolling mean of
    per-codon ratios. This is less sensitive to isolated zero or near-zero
    control positions while retaining a direct IP/control interpretation.
    """
    control_samples = [str(sample) for sample in control_samples]
    ip_samples = [str(sample) for sample in ip_samples]

    if len(control_samples) != len(ip_samples):
        raise ValueError(
            "Matched enrichment requires equal numbers of control and IP samples."
        )
    if not control_samples:
        raise ValueError("Matched enrichment requires at least one replicate pair.")
    if int(window) < 1:
        raise ValueError("Enrichment window must be >= 1.")
    if float(pseudocount) <= 0:
        raise ValueError("Enrichment pseudocount must be > 0.")

    required_columns = control_samples + ip_samples
    missing_columns = [
        sample for sample in required_columns if sample not in raw_gene_rpm.columns
    ]
    if missing_columns:
        raise ValueError(
            "Normalized density table is missing sample(s): {0}".format(
                ", ".join(missing_columns)
            )
        )

    ordered_rpm = raw_gene_rpm.sort_index()
    replicate_ratio = pd.DataFrame(index=ordered_rpm.index, dtype=float)
    window = int(window)
    pseudocount = float(pseudocount)

    for control_sample, ip_sample in zip(control_samples, ip_samples):
        control_signal = ordered_rpm[control_sample].astype(float) + pseudocount
        ip_signal = ordered_rpm[ip_sample].astype(float) + pseudocount

        if window > 1:
            control_signal = control_signal.rolling(
                window=window,
                center=True,
                min_periods=1,
            ).sum()
            ip_signal = ip_signal.rolling(
                window=window,
                center=True,
                min_periods=1,
            ).sum()

        ratio_name = "{0}_{1}".format(ip_sample, control_sample)
        replicate_ratio[ratio_name] = ip_signal.div(control_signal)

    consensus_profile = pd.DataFrame(index=replicate_ratio.index)
    consensus_profile["enrich"] = replicate_ratio.median(axis=1)
    if threshold is not None:
        consensus_profile["support"] = (
            replicate_ratio >= float(threshold)
        ).mean(axis=1)

    return replicate_ratio, consensus_profile
