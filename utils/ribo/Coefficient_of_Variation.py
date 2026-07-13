#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.003
# Function: Calculate transcript-level CDS coefficient of variation from RPF density.
# Input: RPF density file in JSONL or TXT format, optional transcript list, and optional sample-group table.
# Output: Gene-level CoV table, outlier table, group comparison table, fitted mean-CoV curves, and summary JSON.

"""Core functions for transcript-level RPF coefficient-of-variation analysis."""

from __future__ import annotations

import json
import math
import os
import re
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from typing import Any

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import RPFs


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
RPM_SCALE = 1_000_000.0
FIT_MIN_POINTS = 20
PLOT_POINT_LIMIT = 50000


def _normal_two_sided_p(z_value: float) -> float:
    """Return a two-sided normal-approximation p-value."""
    if not np.isfinite(z_value):
        return np.nan
    return float(math.erfc(abs(float(z_value)) / np.sqrt(2.0)))


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Return one-based average ranks with tie handling."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        average_rank = 0.5 * ((start + 1) + stop)
        ranks[order[start:stop]] = average_rank
        start = stop
    return ranks


def _paired_t_normal_approx(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Calculate a paired t statistic with a normal-approximation p-value."""
    diff = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    diff = diff[np.isfinite(diff)]
    if len(diff) < 2:
        return np.nan, np.nan
    sd = float(np.std(diff, ddof=1))
    if sd == 0:
        statistic = 0.0 if float(np.mean(diff)) == 0 else np.sign(np.mean(diff)) * np.inf
    else:
        statistic = float(np.mean(diff) / (sd / np.sqrt(len(diff))))
    return statistic, _normal_two_sided_p(statistic)


def _wilcoxon_normal_approx(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Calculate a paired Wilcoxon signed-rank normal approximation."""
    diff = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    diff = diff[np.isfinite(diff) & (diff != 0)]
    n = len(diff)
    if n == 0:
        return 0.0, 1.0
    abs_diff = np.abs(diff)
    ranks = _average_ranks(abs_diff)
    w_plus = float(ranks[diff > 0].sum())
    w_minus = float(ranks[diff < 0].sum())
    statistic = min(w_plus, w_minus)
    mean_w = n * (n + 1) / 4.0
    _, counts = np.unique(abs_diff, return_counts=True)
    tie_term = float(np.sum(counts * (counts + 1) * (2 * counts + 1)))
    variance = (n * (n + 1) * (2 * n + 1) - 0.5 * tie_term) / 24.0
    if variance <= 0:
        return statistic, 1.0
    correction = 0.5 * np.sign(w_plus - mean_w)
    z_value = (w_plus - mean_w - correction) / np.sqrt(variance)
    return statistic, _normal_two_sided_p(z_value)


def _ks_2sample_asymptotic(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Calculate a two-sample KS statistic and asymptotic p-value."""
    x = np.sort(np.asarray(x, dtype=float))
    y = np.sort(np.asarray(y, dtype=float))
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    values = np.sort(np.unique(np.concatenate([x, y])))
    cdf_x = np.searchsorted(x, values, side="right") / len(x)
    cdf_y = np.searchsorted(y, values, side="right") / len(y)
    statistic = float(np.max(np.abs(cdf_x - cdf_y)))
    effective_n = len(x) * len(y) / (len(x) + len(y))
    if statistic <= 0 or effective_n <= 0:
        return statistic, 1.0
    lam = (np.sqrt(effective_n) + 0.12 + 0.11 / np.sqrt(effective_n)) * statistic
    terms = [(-1) ** (k - 1) * np.exp(-2.0 * (k * lam) ** 2) for k in range(1, 101)]
    p_value = float(np.clip(2.0 * np.sum(terms), 0.0, 1.0))
    return statistic, p_value


def _fit_nonnegative_mean_cov(mean: np.ndarray, cov: np.ndarray) -> tuple[float, float]:
    """Fit CV^2 = alpha + beta / mean with alpha and beta constrained >= 0."""
    x = np.asarray(mean, dtype=float)
    y = np.square(np.asarray(cov, dtype=float))
    design = np.column_stack([np.ones_like(x), 1.0 / x])
    candidates = []
    full, *_ = np.linalg.lstsq(design, y, rcond=None)
    if np.all(full >= 0):
        candidates.append(full)
    alpha_only = np.array([max(float(np.mean(y)), 0.0), 0.0])
    candidates.append(alpha_only)
    inv_mean = design[:, 1]
    denominator = float(np.dot(inv_mean, inv_mean))
    beta = max(float(np.dot(inv_mean, y) / denominator), 0.0) if denominator > 0 else 0.0
    candidates.append(np.array([0.0, beta]))
    candidates.append(np.array([0.0, 0.0]))
    best = min(candidates, key=lambda p: float(np.sum((y - design @ p) ** 2)))
    return float(best[0]), float(best[1])


class CoV(object):
    """Calculate positional RPF coefficient of variation within CDS regions.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_CoV``.
    """

    def __init__(self, args):
        self.rpf = args.rpf
        self.output = args.output
        self.transcript = args.list
        self.group_file = args.group

        self.site = args.site
        self.frame = args.frame
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts
        self.norm = args.normal
        self.ddof = args.ddof
        self.min_codons = args.min_codons
        self.thread = max(1, int(args.thread))

        self.remove_outlier = args.remove_outlier
        self.outlier_iqr = args.outlier_iqr
        self.outlier_window = args.outlier_window
        self.outlier_local_fold = args.outlier_local_fold

        self.fit_model = args.fit_model
        self.fit_quantile = args.fit_quantile
        self.plot_transform = args.plot_transform

        self.rpf_data = None
        self.raw_rpf = None
        self.merged_rpf = None
        self.sample_name: list[str] = []
        self.sample_num = 0
        self.total_rpf_num = None
        self.file_format = None

        self.sample_high_genes: dict[str, pd.Index] = {}
        self.cov_long = pd.DataFrame()
        self.cov_wide = pd.DataFrame()
        self.outliers = pd.DataFrame()
        self.group = None
        self.group_summary = pd.DataFrame()
        self.group_test = pd.DataFrame()
        self.fit_table = pd.DataFrame()
        self.output_files = OrderedDict()

    # ------------------------------------------------------------------
    # Input
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import JSONL or TXT RPF density through ``RPFData``."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf,
            sample_name=None,
            gene=self.transcript,
            tis=self.tis,
            tts=self.tts,
        )
        self.raw_rpf = self.rpf_data.raw_rpf
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf_num = self.rpf_data.total_rpf_num
        self.file_format = self.rpf_data.file_format

        self.merged_rpf = self.rpf_data.get_frame(frame=self.frame)
        self.merged_rpf = RPFs.shift_site(
            self.merged_rpf,
            self.sample_name,
            RPFs.set_codon_shift(self.site),
        )
        self.merged_rpf = self.merged_rpf.loc[
            self.merged_rpf["region"].eq("cds"), BASE_COLUMNS + self.sample_name
        ].copy()
        self.merged_rpf.reset_index(drop=True, inplace=True)

        if self.merged_rpf.empty:
            raise ValueError("No CDS density remains after filtering.")

        # Select high-expression transcripts from raw counts so that -m keeps
        # the same meaning regardless of optional RPM reporting.
        self._set_sample_high_genes()

        if self.norm:
            for sample in self.sample_name:
                total = float(self.total_rpf_num.get(sample, 0.0))
                if total > 0:
                    self.merged_rpf[sample] = (
                        self.merged_rpf[sample].astype(float) * RPM_SCALE / total
                    )
        print(
            "Imported {samples} sample(s), format={fmt}, CDS rows={rows:,}.".format(
                samples=self.sample_num,
                fmt=self.file_format,
                rows=len(self.merged_rpf),
            ),
            flush=True,
        )

    def _set_sample_high_genes(self) -> None:
        """Select high-expression transcripts independently for each sample."""
        for sample in self.sample_name:
            counts = self.merged_rpf.groupby("name", sort=False)[sample].sum()
            keep = counts.loc[counts >= float(self.rpf_num)].index
            self.sample_high_genes[sample] = keep
            print(
                f"Sample {sample}: retained transcripts={len(keep):,}.",
                flush=True,
            )

    def read_group(self) -> None:
        """Read and validate a two-column sample-group table."""
        if not self.group_file:
            self.group = pd.DataFrame(
                {"Name": self.sample_name, "Group": self.sample_name}
            )
            return

        group = pd.read_csv(self.group_file, sep="\t", header=0)
        lower_map = {str(column).lower(): column for column in group.columns}
        if "name" not in lower_map or "group" not in lower_map:
            raise ValueError("Group table must contain Name and Group columns.")
        group = group.rename(
            columns={lower_map["name"]: "Name", lower_map["group"]: "Group"}
        )
        group = group.loc[:, ["Name", "Group"]].dropna().copy()
        group["Name"] = group["Name"].astype(str)
        group["Group"] = group["Group"].astype(str)

        duplicated = group.loc[group["Name"].duplicated(), "Name"].unique().tolist()
        if duplicated:
            raise ValueError("Duplicated sample names in group table: " + ", ".join(duplicated))

        missing = sorted(set(group["Name"]) - set(self.sample_name))
        if missing:
            raise ValueError("Group table contains unknown samples: " + ", ".join(missing))

        missing_assignment = sorted(set(self.sample_name) - set(group["Name"]))
        if missing_assignment:
            raise ValueError(
                "Samples missing from group table: " + ", ".join(missing_assignment)
            )

        self.group = group.set_index("Name").loc[self.sample_name].reset_index()

    # ------------------------------------------------------------------
    # Outlier filtering
    # ------------------------------------------------------------------

    @staticmethod
    def _local_background(values: np.ndarray, window: int) -> np.ndarray:
        """Return symmetric local mean excluding the center position."""
        values = np.asarray(values, dtype=float)
        result = np.full(values.size, np.nan, dtype=float)
        if values.size == 0:
            return result
        window = max(1, int(window))
        valid = np.isfinite(values)
        clean = np.where(valid, values, 0.0)
        prefix_sum = np.concatenate(([0.0], np.cumsum(clean)))
        prefix_n = np.concatenate(([0], np.cumsum(valid.astype(int))))
        idx = np.arange(values.size)
        left = np.maximum(0, idx - window)
        right = np.minimum(values.size, idx + window + 1)
        total = prefix_sum[right] - prefix_sum[left] - clean
        count = prefix_n[right] - prefix_n[left] - valid.astype(int)
        np.divide(total, count, out=result, where=count > 0)
        return result

    @staticmethod
    def _robust_upper_cutoff(values: pd.Series, multiplier: float) -> float:
        """Return a conservative robust upper cutoff on the log1p scale."""
        positive = values.loc[values > 0].astype(float)
        if positive.size < 4:
            return float("inf")
        x = np.log1p(positive.to_numpy(dtype=float))
        q1, q3 = np.percentile(x, [25, 75])
        median = float(np.median(x))
        mad_sigma = 1.4826 * float(np.median(np.abs(x - median)))
        cutoffs = []
        if q3 > q1:
            cutoffs.append(float(q3 + multiplier * (q3 - q1)))
        if mad_sigma > 0:
            cutoffs.append(float(median + multiplier * mad_sigma))
        return float(np.expm1(max(cutoffs))) if cutoffs else float("inf")

    def _detect_sample_outliers(
        self, sample_table: pd.DataFrame, sample: str
    ) -> tuple[pd.Index, pd.DataFrame]:
        """Detect isolated extreme transcript-position pileups."""
        columns = BASE_COLUMNS + [
            "Sample",
            "RawDensity",
            "LocalBackground",
            "LocalFold",
            "OutlierCutoff",
            "OutlierMethod",
        ]
        if not self.remove_outlier or sample_table.empty:
            return pd.Index([]), pd.DataFrame(columns=columns)

        cutoff = self._robust_upper_cutoff(sample_table[sample], self.outlier_iqr)
        if not np.isfinite(cutoff):
            return pd.Index([]), pd.DataFrame(columns=columns)

        candidate = sample_table.loc[sample_table[sample] > cutoff].copy()
        if candidate.empty:
            return pd.Index([]), pd.DataFrame(columns=columns)

        candidate_groups = candidate.groupby("name", sort=False).groups
        local = pd.Series(np.nan, index=candidate.index, dtype=float)
        restricted = sample_table.loc[sample_table["name"].isin(candidate_groups.keys())]
        for gene, gene_df in restricted.groupby("name", sort=False):
            background = self._local_background(
                gene_df[sample].to_numpy(dtype=float), self.outlier_window
            )
            lookup = pd.Series(background, index=gene_df.index)
            indices = candidate_groups.get(gene, [])
            if len(indices):
                local.loc[indices] = lookup.loc[indices].to_numpy()

        candidate["LocalBackground"] = local.reindex(candidate.index).to_numpy()
        candidate["LocalFold"] = (
            (candidate[sample].astype(float) + 1.0)
            / (candidate["LocalBackground"].fillna(0.0) + 1.0)
        )
        candidate = candidate.loc[
            candidate["LocalFold"] >= float(self.outlier_local_fold)
        ].copy()
        if candidate.empty:
            return pd.Index([]), pd.DataFrame(columns=columns)

        candidate["Sample"] = sample
        candidate["RawDensity"] = candidate[sample].astype(float)
        candidate["OutlierCutoff"] = cutoff
        candidate["OutlierMethod"] = "log1p_robust_cutoff_and_symmetric_local_peak"
        return candidate.index, candidate.reindex(columns=columns)

    # ------------------------------------------------------------------
    # CoV calculation
    # ------------------------------------------------------------------

    def _calculate_sample(self, sample: str) -> tuple[str, pd.DataFrame, pd.DataFrame]:
        """Calculate transcript-level CoV for one sample."""
        keep = self.sample_high_genes[sample]
        table = self.merged_rpf.loc[
            self.merged_rpf["name"].isin(keep), BASE_COLUMNS + [sample]
        ].copy()
        table.reset_index(drop=True, inplace=True)

        outlier_index, outlier_df = self._detect_sample_outliers(table, sample)
        if len(outlier_index):
            table.loc[outlier_index, sample] = np.nan

        grouped = table.groupby("name", sort=False)[sample]
        result = grouped.agg(
            CodonCount="size",
            ObservedCodonCount="count",
            PositiveCodonCount=lambda x: int((x > 0).sum()),
            Sum="sum",
            Mean="mean",
            SD=lambda x: x.std(ddof=self.ddof),
        ).reset_index()

        result["CoV"] = np.divide(
            result["SD"].to_numpy(dtype=float),
            result["Mean"].to_numpy(dtype=float),
            out=np.full(len(result), np.nan, dtype=float),
            where=result["Mean"].to_numpy(dtype=float) > 0,
        )
        result = result.loc[
            (result["ObservedCodonCount"] >= int(self.min_codons))
            & np.isfinite(result["CoV"])
        ].copy()
        result.insert(1, "Sample", sample)
        return sample, result, outlier_df

    def calculate_cov(self) -> None:
        """Calculate sample-specific transcript CoV values."""
        workers = min(self.thread, max(1, self.sample_num))
        results: dict[str, pd.DataFrame] = {}
        outliers: dict[str, pd.DataFrame] = {}

        if workers == 1:
            for sample in self.sample_name:
                name, result, outlier = self._calculate_sample(sample)
                results[name] = result
                outliers[name] = outlier
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(self._calculate_sample, sample): sample
                    for sample in self.sample_name
                }
                for future in as_completed(futures):
                    name, result, outlier = future.result()
                    results[name] = result
                    outliers[name] = outlier

        self.cov_long = pd.concat(
            [results[sample] for sample in self.sample_name], ignore_index=True
        )
        outlier_tables = [outliers[s] for s in self.sample_name if not outliers[s].empty]
        self.outliers = (
            pd.concat(outlier_tables, ignore_index=True)
            if outlier_tables
            else pd.DataFrame()
        )

        metrics = [
            "CodonCount",
            "ObservedCodonCount",
            "PositiveCodonCount",
            "Sum",
            "Mean",
            "SD",
            "CoV",
        ]
        wide_parts = []
        for sample in self.sample_name:
            now = results[sample].set_index("name")[metrics].copy()
            now.columns = [f"{sample}_{metric}" for metric in metrics]
            wide_parts.append(now)
        self.cov_wide = pd.concat(wide_parts, axis=1, join="outer")
        self.cov_wide.index.name = "name"

    # ------------------------------------------------------------------
    # Group statistics and curve fitting
    # ------------------------------------------------------------------

    @staticmethod
    def _fit_function(mean: np.ndarray, alpha: float, beta: float) -> np.ndarray:
        """Mean-CoV model: CV = sqrt(alpha + beta / mean)."""
        mean = np.asarray(mean, dtype=float)
        return np.sqrt(np.maximum(alpha + beta / mean, 0.0))

    def _group_gene_values(self, group_name: str) -> pd.DataFrame:
        """Return per-gene group-average mean and CoV."""
        samples = self.group.loc[self.group["Group"] == group_name, "Name"].tolist()
        mean_columns = [f"{sample}_Mean" for sample in samples]
        cov_columns = [f"{sample}_CoV" for sample in samples]
        result = pd.DataFrame(index=self.cov_wide.index)
        result["Mean"] = self.cov_wide[mean_columns].mean(axis=1, skipna=True)
        result["CoV"] = self.cov_wide[cov_columns].mean(axis=1, skipna=True)
        result["ReplicateCount"] = self.cov_wide[cov_columns].notna().sum(axis=1)
        return result.replace([np.inf, -np.inf], np.nan).dropna(subset=["Mean", "CoV"])

    def compare_groups(self) -> None:
        """Compare group-level CoV distributions using paired gene values."""
        if self.group is None:
            self.read_group()

        summaries = []
        group_values: dict[str, pd.DataFrame] = {}
        for group_name in self.group["Group"].drop_duplicates():
            values = self._group_gene_values(group_name)
            group_values[group_name] = values
            summaries.append(
                {
                    "Group": group_name,
                    "SampleCount": int((self.group["Group"] == group_name).sum()),
                    "GeneCount": int(len(values)),
                    "MedianMean": float(values["Mean"].median()),
                    "MedianCoV": float(values["CoV"].median()),
                    "MeanCoV": float(values["CoV"].mean()),
                }
            )
        self.group_summary = pd.DataFrame(summaries)

        records = []
        for group1, group2 in combinations(group_values.keys(), 2):
            merged = group_values[group1][["CoV"]].join(
                group_values[group2][["CoV"]],
                how="inner",
                lsuffix="_1",
                rsuffix="_2",
            ).dropna()
            x = merged["CoV_1"].to_numpy(dtype=float)
            y = merged["CoV_2"].to_numpy(dtype=float)
            if len(merged) < 2:
                continue
            paired_t_stat, paired_t_p = _paired_t_normal_approx(x, y)
            wilcoxon_stat, wilcoxon_p = _wilcoxon_normal_approx(x, y)
            ks_stat, ks_p = _ks_2sample_asymptotic(x, y)
            records.append(
                {
                    "Group1": group1,
                    "Group2": group2,
                    "PairedGeneCount": int(len(merged)),
                    "MeanCoVGroup1": float(np.mean(x)),
                    "MeanCoVGroup2": float(np.mean(y)),
                    "MeanDifference": float(np.mean(x - y)),
                    "PairedTStatistic": float(paired_t_stat),
                    "PairedTPValue": float(paired_t_p),
                    "WilcoxonStatistic": float(wilcoxon_stat),
                    "WilcoxonPValue": float(wilcoxon_p),
                    "KSStatistic": float(ks_stat),
                    "KSPValue": float(ks_p),
                }
            )
        self.group_test = pd.DataFrame(records)

    def fit_mean_cov(self) -> None:
        """Fit the corrected mean-CoV curve for each group."""
        if self.group is None:
            self.read_group()
        records = []
        for group_name in self.group["Group"].drop_duplicates():
            values = self._group_gene_values(group_name)
            values = values.loc[(values["Mean"] > 0) & (values["CoV"] > 0)].copy()
            if values.empty:
                continue

            lower = float(self.fit_quantile)
            upper = 1.0 - lower
            if 0 < lower < 0.5 and len(values) >= FIT_MIN_POINTS:
                mean_low, mean_high = values["Mean"].quantile([lower, upper])
                cov_low, cov_high = values["CoV"].quantile([lower, upper])
                fit_values = values.loc[
                    values["Mean"].between(mean_low, mean_high)
                    & values["CoV"].between(cov_low, cov_high)
                ].copy()
            else:
                fit_values = values

            alpha = beta = r_squared = np.nan
            success = False
            message = "insufficient_points"
            if len(fit_values) >= FIT_MIN_POINTS:
                x = fit_values["Mean"].to_numpy(dtype=float)
                y = fit_values["CoV"].to_numpy(dtype=float)
                try:
                    alpha, beta = _fit_nonnegative_mean_cov(x, y)
                    predicted = self._fit_function(x, alpha, beta)
                    residual_ss = float(np.sum((y - predicted) ** 2))
                    total_ss = float(np.sum((y - np.mean(y)) ** 2))
                    r_squared = 1.0 - residual_ss / total_ss if total_ss > 0 else np.nan
                    success = True
                    message = "success"
                except Exception as exc:  # pragma: no cover - data dependent
                    message = str(exc)

            records.append(
                {
                    "Group": group_name,
                    "TotalGeneCount": int(len(values)),
                    "FitGeneCount": int(len(fit_values)),
                    "Alpha": alpha,
                    "Beta": beta,
                    "R2": r_squared,
                    "Success": success,
                    "Message": message,
                    "Model": "CV = sqrt(alpha + beta / mean); NumPy constrained least squares",
                }
            )
        self.fit_table = pd.DataFrame(records)

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def output_tables(self) -> None:
        """Write CoV, group-comparison, fitting, and outlier tables."""
        cov_file = self.output + "_CoV.txt"
        self.cov_wide.to_csv(cov_file, sep="\t", index=True)
        self.output_files["cov_table"] = cov_file

        cov_long_file = self.output + "_CoV.long.txt"
        self.cov_long.to_csv(cov_long_file, sep="\t", index=False)
        self.output_files["cov_long_table"] = cov_long_file

        if self.group_summary is not None and not self.group_summary.empty:
            path = self.output + "_CoV.group_summary.txt"
            self.group_summary.to_csv(path, sep="\t", index=False)
            self.output_files["group_summary"] = path

        if self.group_test is not None and not self.group_test.empty:
            path = self.output + "_compared_CoV.txt"
            self.group_test.to_csv(path, sep="\t", index=False)
            self.output_files["group_comparison"] = path

        if self.fit_table is not None and not self.fit_table.empty:
            path = self.output + "_CoV.fit_parameters.txt"
            self.fit_table.to_csv(path, sep="\t", index=False)
            self.output_files["fit_parameters"] = path

        if self.remove_outlier:
            path = self.output + "_CoV.outliers.txt"
            self.outliers.to_csv(path, sep="\t", index=False)
            self.output_files["outlier_table"] = path

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_name(value: str) -> str:
        return re.sub(r"[^0-9A-Za-z._-]+", "_", str(value)).strip("_") or "group"

    def draw_fit_plot(self) -> None:
        """Draw all group mean-CoV scatters with corrected fitted curves."""
        if self.group is None:
            self.read_group()
        if self.fit_table is None or self.fit_table.empty:
            self.fit_mean_cov()

        groups = self.group["Group"].drop_duplicates().tolist()
        n_groups = len(groups)
        ncols = min(3, max(1, n_groups))
        nrows = int(np.ceil(n_groups / ncols))
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(5.2 * ncols, 4.3 * nrows),
            squeeze=False,
        )

        for ax, group_name in zip(axes.flat, groups):
            values = self._group_gene_values(group_name)
            values = values.loc[(values["Mean"] > 0) & (values["CoV"] > 0)].copy()
            if len(values) > PLOT_POINT_LIMIT:
                values = values.sample(PLOT_POINT_LIMIT, random_state=0)

            ax.scatter(
                np.log2(values["Mean"]),
                np.log2(values["CoV"]),
                s=8,
                alpha=0.22,
                edgecolors="none",
            )

            fit_row = self.fit_table.loc[self.fit_table["Group"] == group_name]
            if not fit_row.empty and bool(fit_row.iloc[0]["Success"]):
                alpha = float(fit_row.iloc[0]["Alpha"])
                beta = float(fit_row.iloc[0]["Beta"])
                x_grid = np.geomspace(values["Mean"].min(), values["Mean"].max(), 300)
                y_grid = self._fit_function(x_grid, alpha, beta)
                ax.plot(np.log2(x_grid), np.log2(y_grid), linewidth=2.0)
                label = "α={:.3g}, β={:.3g}, R²={:.3f}".format(
                    alpha, beta, float(fit_row.iloc[0]["R2"])
                )
                ax.text(
                    0.03,
                    0.97,
                    label,
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=8,
                    bbox={"boxstyle": "round,pad=0.25", "fc": "white", "alpha": 0.8, "ec": "none"},
                )

            ax.set_title(str(group_name))
            ax.set_xlabel("log2(mean RPF density)")
            ax.set_ylabel("log2(CoV)")
            ax.grid(linewidth=0.4, alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        for ax in axes.flat[n_groups:]:
            ax.set_axis_off()

        fig.suptitle("Mean–CoV relationship", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        pdf = self.output + "_CoV_fitplot.pdf"
        png = self.output + "_CoV_fitplot.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["fitplot_pdf"] = pdf
        self.output_files["fitplot_png"] = png

    def draw_distribution_plot(self) -> None:
        """Draw group-level CoV distributions and ECDF curves."""
        if self.group is None:
            self.read_group()

        rows = []
        for group_name in self.group["Group"].drop_duplicates():
            values = self._group_gene_values(group_name)["CoV"].dropna()
            rows.extend({"Group": group_name, "CoV": value} for value in values)
        plot_df = pd.DataFrame(rows)
        if plot_df.empty:
            return

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"wspace": 0.28})
        groups = plot_df["Group"].drop_duplicates().tolist()
        data = [plot_df.loc[plot_df["Group"] == group, "CoV"].to_numpy() for group in groups]
        axes[0].boxplot(data, labels=groups, showfliers=False)
        axes[0].set_ylabel("Coefficient of variation")
        axes[0].set_title("Group CoV distribution")
        axes[0].tick_params(axis="x", rotation=30)

        for group in groups:
            x = np.sort(plot_df.loc[plot_df["Group"] == group, "CoV"].to_numpy())
            if x.size:
                y = np.arange(1, x.size + 1) / x.size
                axes[1].plot(x, y, linewidth=1.5, label=group)
        axes[1].set_xlabel("Coefficient of variation")
        axes[1].set_ylabel("Empirical cumulative probability")
        axes[1].set_title("CoV empirical distributions")
        axes[1].legend(frameon=False, fontsize=8)

        for ax in axes:
            ax.grid(linewidth=0.4, alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        pdf = self.output + "_CoV_distribution.pdf"
        png = self.output + "_CoV_distribution.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["distribution_pdf"] = pdf
        self.output_files["distribution_png"] = png

    def write_summary(self) -> None:
        """Write a machine-readable analysis summary."""
        summary_path = self.output + "_CoV.summary.json"
        self.output_files["summary_json"] = summary_path
        summary: dict[str, Any] = OrderedDict(
            [
                ("tool", "rpf_CoV"),
                ("version", "0.2.8-dev.002"),
                ("input_rpf", os.path.abspath(self.rpf)),
                ("input_format", self.file_format),
                ("sample_count", self.sample_num),
                ("samples", self.sample_name),
                ("parameters", OrderedDict([
                    ("site", self.site),
                    ("frame", self.frame),
                    ("min", self.rpf_num),
                    ("tis", self.tis),
                    ("tts", self.tts),
                    ("normal", self.norm),
                    ("ddof", self.ddof),
                    ("min_codons", self.min_codons),
                    ("thread", self.thread),
                    ("remove_outlier", self.remove_outlier),
                    ("outlier_iqr", self.outlier_iqr),
                    ("outlier_window", self.outlier_window),
                    ("outlier_local_fold", self.outlier_local_fold),
                    ("fit_model", self.fit_model),
                    ("fit_quantile", self.fit_quantile),
                ])),
                ("sample_retained_gene_count", {
                    sample: int((self.cov_long["Sample"] == sample).sum())
                    for sample in self.sample_name
                }),
                ("fit_results", self.fit_table.to_dict(orient="records")),
                ("output_files", self.output_files),
            ]
        )
        with open(summary_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
