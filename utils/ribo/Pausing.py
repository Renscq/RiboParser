#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.005
# Function: Calculate sample-specific codon pausing scores from RPF density data.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Gene-level, gene-codon, codon-level, outlier, summary, and figure files.

"""Core functions for RiboParser codon pausing analysis.

The implementation supports compact JSONL and legacy TXT density files through
``RPFs.RPFData``. Pausing scores are calculated independently for each sample,
and high-expression transcripts are selected using sample-specific CDS counts.

For local-background analysis, the focal codon is excluded from the denominator:

    pausing_score[i] = density[i] / mean(neighbor densities excluding i)

This avoids the downward bias in the legacy centered rolling mean, which included
the focal pileup in its own background.
"""

from __future__ import annotations

import json
import os
import re
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable

from . import RPFs

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Use Arial (with common fallbacks when it is unavailable on the system) for
# every figure produced by this module.
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
    }
)


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
RPM_SCALE = 1_000_000.0
PROGRESS_EVERY = 500
PLOT_CMAP_SEQUENTIAL = "RdBu_r"
PLOT_CMAP_DIVERGING = "RdBu_r"
HEATMAP_ANNOTATION_MAX_SAMPLES = 12


class Pausing(object):
    """Calculate sample-specific codon pausing scores.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Pausing``.
    """

    def __init__(self, args):
        # Input and output.
        self.rpf = args.rpf
        self.output = args.output
        self.transcript = args.list
        self.output_all = args.all

        # Density import and filtering.
        self.site = args.site
        self.frame = args.frame
        self.normal = args.normal
        self.tis = args.tis
        self.tts = args.tts
        self.rpf_num = args.min
        self.background = args.background
        self.exclude_stop = True
        self.individual = args.individual
        self.thread = max(1, int(args.thread))

        # Scaling and plotting.
        self.scale = args.scale
        self.plot_transform = args.plot_transform
        self.rankplot_ncol = max(1, int(getattr(args, "rankplot_ncol", 1)))

        # Outlier filtering.
        self.remove_outlier = args.remove_outlier
        self.outlier_iqr = args.outlier_iqr
        self.outlier_window = args.outlier_window
        self.outlier_local_fold = args.outlier_local_fold

        # Imported data.
        self.rpf_data = None
        self.raw_rpf = None
        self.merged_rpf = None
        self.sample_name = []
        self.sample_num = 0
        self.total_rpf_num = None
        self.file_format = None

        # Derived data.
        self.sample_high_genes = OrderedDict()
        self.sample_gene_counts = pd.DataFrame()
        self.pausing = None
        self.outliers = pd.DataFrame()
        self.gene_pausing = None
        self.gene_codon_pausing = None
        self.codon_pausing = None
        self.sample_summary = pd.DataFrame()
        self.output_files = OrderedDict()

        self.codon_dict, self.codon_table = RPFs.codon_table()
        # Stop codons are never included in pausing calculation or visualization.
        for codon in ("TAA", "TAG", "TGA"):
            self.codon_dict.pop(codon, None)
        self.codon_table = self.codon_table.loc[self.codon_table["Abbr"] != "*", :]

    # ------------------------------------------------------------------
    # Import and validation
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import and prepare RPF density data."""
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
        self.total_rpf_num = self.rpf_data.total_rpf_num.astype(float)
        self.file_format = self.rpf_data.file_format

        self.merged_rpf = self.rpf_data.get_frame(frame=self.frame)
        shift_num = RPFs.set_codon_shift(self.site)
        self.merged_rpf = RPFs.shift_site(self.merged_rpf, self.sample_name, shift_num)

        self._validate_table()
        self.merged_rpf = self.merged_rpf.loc[
            self.merged_rpf["codon"].isin(self.codon_dict.keys()),
            BASE_COLUMNS + self.sample_name,
        ].copy()

        for sample in self.sample_name:
            self.merged_rpf[sample] = pd.to_numeric(
                self.merged_rpf[sample], errors="coerce"
            ).fillna(0.0).astype(float)

        # Select high-expression transcripts from raw sample-specific CDS counts.
        # This keeps --min semantics independent of optional RPM normalization.
        self._select_sample_high_genes()

        if self.normal:
            for sample in self.sample_name:
                total = float(self.total_rpf_num.get(sample, 0.0))
                if total > 0:
                    self.merged_rpf[sample] *= RPM_SCALE / total
                else:
                    self.merged_rpf[sample] = 0.0

        print(
            "Imported RPF density for {samples} sample(s), format={fmt}, rows={rows:,}.".format(
                samples=self.sample_num,
                fmt=self.file_format,
                rows=len(self.merged_rpf),
            ),
            flush=True,
        )

    # Backward-compatible method name.
    read_rpf = import_rpf

    def _validate_table(self) -> None:
        """Validate the imported codon-level density table."""
        if self.merged_rpf is None or self.merged_rpf.empty:
            raise ValueError("RPF density table is empty after import and coordinate filtering.")

        missing = [
            column
            for column in BASE_COLUMNS + self.sample_name
            if column not in self.merged_rpf.columns
        ]
        if missing:
            raise ValueError(
                "RPF density table is missing required column(s): " + ", ".join(missing)
            )

    def _select_sample_high_genes(self) -> None:
        """Select high-expression transcripts independently for each sample."""
        cds = self.merged_rpf.loc[self.merged_rpf["region"] == "cds", :]
        if cds.empty:
            raise ValueError("No CDS codons remain after TIS/TTS filtering.")

        gene_counts = cds.groupby("name", sort=False)[self.sample_name].sum()
        self.sample_gene_counts = gene_counts

        records = []
        for sample in self.sample_name:
            high_genes = gene_counts.index[gene_counts[sample] >= float(self.rpf_num)]
            self.sample_high_genes[sample] = pd.Index(high_genes)
            records.append(
                {
                    "Sample": sample,
                    "HighGeneCount": int(len(high_genes)),
                    "MinimumCDSCount": float(self.rpf_num),
                }
            )
        self.sample_summary = pd.DataFrame.from_records(records)

    # ------------------------------------------------------------------
    # Pausing calculation
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
        """Divide arrays while converting invalid results to NaN."""
        result = np.full(numerator.shape, np.nan, dtype=float)
        valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0)
        np.divide(numerator, denominator, out=result, where=valid)
        return result

    @staticmethod
    def _local_background(values: np.ndarray, window: int) -> np.ndarray:
        """Calculate a symmetric local mean excluding the focal codon.

        Parameters
        ----------
        values : numpy.ndarray
            One-dimensional codon density vector.
        window : int
            Number of neighboring codons on each side.

        Returns
        -------
        numpy.ndarray
            Local background for every position. Boundary positions use only
            available neighbors and are never padded with artificial zeros.
        """
        values = np.asarray(values, dtype=float)
        n = values.size
        if n == 0:
            return np.asarray([], dtype=float)

        valid_values = np.where(np.isfinite(values), values, 0.0)
        valid_counts = np.isfinite(values).astype(float)
        prefix_sum = np.concatenate(([0.0], np.cumsum(valid_values)))
        prefix_count = np.concatenate(([0.0], np.cumsum(valid_counts)))

        index = np.arange(n)
        left = np.maximum(0, index - int(window))
        right = np.minimum(n, index + int(window) + 1)

        neighbor_sum = prefix_sum[right] - prefix_sum[left] - valid_values
        neighbor_count = prefix_count[right] - prefix_count[left] - valid_counts

        background = np.full(n, np.nan, dtype=float)
        np.divide(
            neighbor_sum,
            neighbor_count,
            out=background,
            where=neighbor_count > 0,
        )
        return background

    @staticmethod
    def _robust_upper_cutoff(values: pd.Series, multiplier: float) -> float:
        """Return a conservative robust upper cutoff on the log1p scale."""
        positive = values.loc[values > 0].astype(float)
        if positive.size < 4:
            return float("inf")

        log_values = np.log1p(positive.to_numpy(dtype=float))
        q1 = float(np.percentile(log_values, 25))
        q3 = float(np.percentile(log_values, 75))
        iqr = q3 - q1
        median = float(np.median(log_values))
        mad = float(np.median(np.abs(log_values - median)))
        mad_sigma = 1.4826 * mad

        cutoffs = []
        if iqr > 0:
            cutoffs.append(q3 + float(multiplier) * iqr)
        if mad_sigma > 0:
            cutoffs.append(median + float(multiplier) * mad_sigma)
        if not cutoffs:
            return float("inf")
        return float(np.expm1(max(cutoffs)))

    def _detect_sample_outliers(
        self,
        sample_table: pd.DataFrame,
        sample: str,
    ) -> tuple[pd.Index, pd.DataFrame]:
        """Detect isolated extreme pileups, restricting local scans to candidate genes."""
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

        candidate = sample_table.loc[sample_table[sample] > cutoff, :].copy()
        if candidate.empty:
            return pd.Index([]), pd.DataFrame(columns=columns)

        # Only genes containing a global candidate require local-background scans.
        candidate_genes = set(candidate["name"].astype(str))
        local_background = pd.Series(np.nan, index=candidate.index, dtype=float)
        candidate_by_gene = candidate.groupby("name", sort=False).groups

        restricted = sample_table.loc[sample_table["name"].astype(str).isin(candidate_genes), :]
        for gene_name, gene_df in restricted.groupby("name", sort=False):
            background = self._local_background(
                gene_df[sample].to_numpy(dtype=float),
                self.outlier_window,
            )
            position_lookup = pd.Series(background, index=gene_df.index)
            candidate_index = candidate_by_gene.get(gene_name, [])
            if len(candidate_index) > 0:
                local_background.loc[candidate_index] = position_lookup.loc[candidate_index].to_numpy()

        candidate["LocalBackground"] = local_background.reindex(candidate.index).to_numpy()
        candidate["LocalFold"] = (
            (candidate[sample].astype(float) + 1.0)
            / (candidate["LocalBackground"].fillna(0.0).astype(float) + 1.0)
        )
        candidate = candidate.loc[
            candidate["LocalFold"] >= float(self.outlier_local_fold), :
        ].copy()
        if candidate.empty:
            return pd.Index([]), pd.DataFrame(columns=columns)

        candidate["Sample"] = sample
        candidate["RawDensity"] = candidate[sample].astype(float)
        candidate["OutlierCutoff"] = cutoff
        candidate["OutlierMethod"] = "log1p_robust_cutoff_and_symmetric_local_peak"
        return candidate.index, candidate.reindex(columns=columns)

    def _calculate_sample_pausing(
        self,
        sample: str,
    ) -> tuple[str, pd.DataFrame, pd.DataFrame, dict[str, int]]:
        """Calculate one sample using vectorized global background or preallocated local arrays."""
        high_genes = self.sample_high_genes[sample]
        sample_table = self.merged_rpf.loc[
            self.merged_rpf["name"].isin(high_genes),
            BASE_COLUMNS + [sample],
        ].copy()
        sample_table.reset_index(drop=True, inplace=True)

        outlier_index, outlier_df = self._detect_sample_outliers(sample_table, sample)
        if len(outlier_index) > 0:
            sample_table.loc[outlier_index, sample] = np.nan

        values = sample_table[sample].to_numpy(dtype=float)

        if self.background == 0:
            # Fully vectorized gene/sample CDS mean. This replaces the expensive
            # per-gene dataframe-copy and concat workflow.
            cds_mask = sample_table["region"].eq("cds")
            cds_mean = (
                sample_table.loc[cds_mask, ["name", sample]]
                .groupby("name", sort=False)[sample]
                .mean()
            )
            denominator = sample_table["name"].map(cds_mean).to_numpy(dtype=float)
        else:
            # Allocate one output array and fill gene slices in place.
            denominator = np.full(len(sample_table), np.nan, dtype=float)
            for positions in sample_table.groupby("name", sort=False).indices.values():
                positions = np.asarray(positions, dtype=int)
                denominator[positions] = self._local_background(
                    values[positions],
                    self.background,
                )

        score = self._safe_divide(values, denominator)
        sample_result = sample_table.loc[:, BASE_COLUMNS].copy()
        sample_result["Sample"] = sample
        sample_result["RawDensity"] = values
        sample_result["BackgroundDensity"] = denominator
        sample_result["PausingScore"] = score

        summary = {
            "OutlierCount": int(len(outlier_index)),
            "ValidPausingPositionCount": int(np.isfinite(score).sum()),
            "HighGeneCount": int(len(high_genes)),
        }
        return sample, sample_result, outlier_df, summary

    def calculate_pausing(self) -> None:
        """Calculate pausing scores with optional sample-level multithreading."""
        if self.merged_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        worker_count = min(self.thread, max(1, self.sample_num))
        results: dict[str, tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]] = {}

        if worker_count == 1:
            for sample in self.sample_name:
                sample_name, sample_result, outlier_df, summary = self._calculate_sample_pausing(sample)
                results[sample_name] = (sample_result, outlier_df, summary)
        else:
            print(
                "Calculate pausing scores with {workers} sample workers.".format(
                    workers=worker_count
                ),
                flush=True,
            )
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(self._calculate_sample_pausing, sample): sample
                    for sample in self.sample_name
                }
                for future in as_completed(futures):
                    sample_name, sample_result, outlier_df, summary = future.result()
                    results[sample_name] = (sample_result, outlier_df, summary)
                    print(
                        "Sample {sample}: high genes={genes:,}, valid scores={scores:,}, outliers={outliers:,}.".format(
                            sample=sample_name,
                            genes=summary["HighGeneCount"],
                            scores=summary["ValidPausingPositionCount"],
                            outliers=summary["OutlierCount"],
                        ),
                        flush=True,
                    )

        # Restore the original sample order for deterministic output.
        sample_tables = []
        outlier_tables = []
        summary_updates = []
        for sample in self.sample_name:
            sample_result, outlier_df, summary = results[sample]
            sample_tables.append(sample_result)
            if not outlier_df.empty:
                outlier_tables.append(outlier_df)
            summary_updates.append(
                {
                    "Sample": sample,
                    "OutlierCount": summary["OutlierCount"],
                    "ValidPausingPositionCount": summary["ValidPausingPositionCount"],
                }
            )
            if worker_count == 1:
                print(
                    "Sample {sample}: high genes={genes:,}, valid scores={scores:,}, outliers={outliers:,}.".format(
                        sample=sample,
                        genes=summary["HighGeneCount"],
                        scores=summary["ValidPausingPositionCount"],
                        outliers=summary["OutlierCount"],
                    ),
                    flush=True,
                )

        self.pausing = pd.concat(sample_tables, axis=0, ignore_index=True)
        self.pausing.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.outliers = (
            pd.concat(outlier_tables, axis=0, ignore_index=True)
            if outlier_tables
            else pd.DataFrame(
                columns=BASE_COLUMNS
                + [
                    "Sample",
                    "RawDensity",
                    "LocalBackground",
                    "LocalFold",
                    "OutlierCutoff",
                    "OutlierMethod",
                ]
            )
        )

        update_df = pd.DataFrame.from_records(summary_updates)
        self.sample_summary = self.sample_summary.merge(update_df, on="Sample", how="left")

    # Backward-compatible method name.
    get_pausing_score = calculate_pausing

    # ------------------------------------------------------------------
    # Summary tables
    # ------------------------------------------------------------------

    @staticmethod
    def _flatten_columns(columns: Iterable) -> list[str]:
        """Flatten dataframe MultiIndex columns."""
        flattened = []
        for column in columns:
            if isinstance(column, tuple):
                flattened.append("_".join(str(part) for part in column if str(part)))
            else:
                flattened.append(str(column))
        return flattened

    def summarize_gene_pausing(self) -> None:
        """Summarize CDS pausing scores for each transcript and sample."""
        cds = self.pausing.loc[self.pausing["region"] == "cds", :].copy()
        cds["IsValidCodon"] = cds["RawDensity"].fillna(0.0).gt(0)
        summary = (
            cds.groupby(["name", "Sample"], sort=False)
            .agg(
                CDSCodonCount=("codon", "size"),
                ValidCodonCount=("IsValidCodon", "sum"),
                PausingScoreSum=("PausingScore", "sum"),
                PausingScoreMean=("PausingScore", "mean"),
                PausingScoreMedian=("PausingScore", "median"),
                PausingScoreMax=("PausingScore", "max"),
            )
            .reset_index()
        )
        summary["ValidCodonCount"] = summary["ValidCodonCount"].astype(int)
        self.gene_pausing = summary

    def summarize_gene_codon_pausing(self) -> None:
        """Summarize each codon within each transcript and sample."""
        cds = self.pausing.loc[self.pausing["region"] == "cds", :].copy()
        cds["IsValidCodon"] = cds["RawDensity"].fillna(0.0).gt(0)
        self.gene_codon_pausing = (
            cds.groupby(["name", "codon", "Sample"], sort=False)
            .agg(
                CodonCount=("codon", "size"),
                ValidCodonCount=("IsValidCodon", "sum"),
                PausingScoreSum=("PausingScore", "sum"),
                PausingScoreMean=("PausingScore", "mean"),
            )
            .reset_index()
        )
        self.gene_codon_pausing["ValidCodonCount"] = (
            self.gene_codon_pausing["ValidCodonCount"].astype(int)
        )

    @staticmethod
    def _scale_matrix(matrix: pd.DataFrame, method: str) -> pd.DataFrame:
        """Scale codon pausing scores independently within each sample."""
        scaled = matrix.astype(float).copy()
        if method == "none":
            return scaled
        if method == "minmax":
            maxima = scaled.max(axis=0).replace(0, np.nan)
            return scaled.div(maxima, axis=1).fillna(0.0)
        if method == "zscore":
            means = scaled.mean(axis=0)
            stds = scaled.std(axis=0, ddof=0).replace(0, np.nan)
            return scaled.sub(means, axis=1).div(stds, axis=1).fillna(0.0)
        raise ValueError("Unsupported scale method: {method}".format(method=method))

    def summarize_codon_pausing(self) -> None:
        """Summarize codon pausing scores independently for each sample."""
        cds = self.pausing.loc[self.pausing["region"] == "cds", :].copy()

        total_occurrence = cds.groupby("codon", sort=False).size().rename("TotalCodon")
        raw_count = cds.groupby(["codon", "Sample"], sort=False)["RawDensity"].sum().unstack("Sample")
        valid_count = (
            cds.assign(IsValid=cds["RawDensity"].fillna(0.0) > 0)
            .groupby(["codon", "Sample"], sort=False)["IsValid"]
            .sum()
            .unstack("Sample")
        )

        total_ps = cds.groupby(["codon", "Sample"], sort=False)["PausingScore"].mean().unstack("Sample")

        if self.individual:
            valid_input = cds.loc[cds["RawDensity"].fillna(0.0) > 0, :]
        else:
            # In long format each row belongs to one sample; therefore the
            # cross-sample legacy criterion is represented by a position key.
            key_columns = ["name", "now_nt", "codon"]
            positive_keys = cds.loc[cds["RawDensity"].fillna(0.0) > 0, key_columns].drop_duplicates()
            valid_input = cds.merge(positive_keys, on=key_columns, how="inner")

        valid_ps = (
            valid_input.groupby(["codon", "Sample"], sort=False)["PausingScore"]
            .mean()
            .unstack("Sample")
        )

        codon_index = self.codon_table.index.intersection(total_occurrence.index)
        total_ps = total_ps.reindex(index=codon_index, columns=self.sample_name)
        valid_ps = valid_ps.reindex(index=codon_index, columns=self.sample_name)
        raw_count = raw_count.reindex(index=codon_index, columns=self.sample_name).fillna(0.0)
        valid_count = valid_count.reindex(index=codon_index, columns=self.sample_name).fillna(0).astype(int)

        relative_total = self._scale_matrix(total_ps, self.scale)
        relative_valid = self._scale_matrix(valid_ps, self.scale)

        result = self.codon_table.reindex(codon_index).copy()
        result["TotalCodon"] = total_occurrence.reindex(codon_index).fillna(0).astype(int)

        for sample in self.sample_name:
            result[f"{sample}_valid_codon"] = valid_count[sample]
            result[f"{sample}_rpf_count"] = raw_count[sample]
            result[f"{sample}_absolute_total_ps"] = total_ps[sample]
            result[f"{sample}_absolute_valid_ps"] = valid_ps[sample]
            result[f"{sample}_relative_total_ps"] = relative_total[sample]
            result[f"{sample}_relative_valid_ps"] = relative_valid[sample]

        result.index.name = "Codon"
        result = result.reset_index().sort_values(["Abbr", "Codon"], kind="stable")
        self.codon_pausing = result.reset_index(drop=True)

    # Backward-compatible output wrappers.
    def output_cds_pausing(self) -> None:
        self.summarize_gene_pausing()
        path = self.output + "_cds_pausing_score.txt"
        self.gene_pausing.to_csv(path, sep="\t", index=False)
        self.output_files["gene_pausing_table"] = path

    def output_cds_codon_pausing(self) -> None:
        """Write transcript-codon pausing summaries when ``--all`` is enabled."""
        if not self.output_all:
            return

        self.summarize_gene_codon_pausing()
        path = self.output + "_cds_codon_pausing_score.txt"
        self.gene_codon_pausing.to_csv(path, sep="\t", index=False)
        self.output_files["gene_codon_pausing_table"] = path

    def output_sum_codon_pausing(self) -> None:
        self.summarize_codon_pausing()
        path = self.output + "_sum_codon_pausing_score.txt"
        self.codon_pausing.to_csv(path, sep="\t", index=False)
        self.output_files["codon_pausing_table"] = path

    def output_all_pausing(self) -> None:
        if self.output_all:
            path = self.output + "_all_pausing_score.txt"
            self.pausing.to_csv(path, sep="\t", index=False)
            self.output_files["all_pausing_table"] = path

        if self.remove_outlier:
            path = self.output + "_pausing.outliers.txt"
            self.outliers.to_csv(path, sep="\t", index=False)
            self.output_files["outlier_table"] = path


    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _transform_plot_values(self, values: np.ndarray) -> np.ndarray:
        """Transform non-negative pausing values for plotting only."""
        values = np.asarray(values, dtype=float)
        values = np.where(np.isfinite(values), values, np.nan)
        if self.plot_transform == "none":
            return values
        if self.plot_transform == "sqrt":
            return np.sqrt(np.clip(values, 0.0, None))
        if self.plot_transform in {"log", "log1p"}:
            return np.log1p(np.clip(values, 0.0, None))
        if self.plot_transform == "log2":
            return np.log2(np.clip(values, 0.0, None) + 1.0)
        if self.plot_transform == "log10":
            return np.log10(np.clip(values, 0.0, None) + 1.0)
        raise ValueError("Unsupported plot transform: " + self.plot_transform)

    def _codon_matrix(self, category: str, relative: bool = False) -> pd.DataFrame:
        """Return codon-by-sample matrix used by heatmaps."""
        prefix = "relative" if relative else "absolute"
        columns = [f"{sample}_{prefix}_{category}_ps" for sample in self.sample_name]
        matrix = self.codon_pausing.set_index("Codon")[columns].copy()
        matrix.columns = self.sample_name
        if not relative:
            matrix.loc[:, :] = self._transform_plot_values(matrix.to_numpy(dtype=float))
        return matrix

    @staticmethod
    def _heatmap_limits(matrix_a: pd.DataFrame, matrix_b: pd.DataFrame, diverging: bool) -> tuple[float, float]:
        """Return shared robust color limits for paired heatmaps."""
        values = np.concatenate([matrix_a.to_numpy(dtype=float).ravel(), matrix_b.to_numpy(dtype=float).ravel()])
        values = values[np.isfinite(values)]
        if values.size == 0:
            return (0.0, 1.0)
        if diverging:
            bound = float(np.nanpercentile(np.abs(values), 98))
            bound = bound if bound > 0 else 1.0
            return (-bound, bound)
        upper = float(np.nanpercentile(values, 98))
        return (0.0, upper if upper > 0 else 1.0)

    @staticmethod
    def _adaptive_correlation_limits(values: np.ndarray) -> tuple[float, float]:
        """Return an adaptive correlation range emphasizing sample differences."""
        values = np.asarray(values, dtype=float)
        if values.ndim != 2 or values.shape[0] < 2:
            return -1.0, 1.0

        mask = ~np.eye(values.shape[0], dtype=bool)
        off_diagonal = values[mask]
        off_diagonal = off_diagonal[np.isfinite(off_diagonal)]
        if off_diagonal.size == 0:
            return -1.0, 1.0

        lower = float(np.nanpercentile(off_diagonal, 2.0))
        upper = float(np.nanpercentile(off_diagonal, 98.0))
        spread = max(upper - lower, 0.01)
        margin = max(spread * 0.15, 0.005)
        vmin = max(-1.0, min(lower - margin, 0.95))
        vmin = float(np.floor(vmin * 100.0) / 100.0)
        if vmin >= 0.99:
            vmin = 0.98
        return vmin, 1.0

    def draw_pausing_corr(self) -> None:
        """Draw sample correlation based on absolute valid codon pausing scores."""
        matrix = self._codon_matrix("valid", relative=False)
        corr = matrix.corr(method="pearson")

        txt = self.output + "_pausing_corr.txt"
        corr.to_csv(txt, sep="\t", index=True)
        self.output_files["pausing_correlation_table"] = txt

        size = max(5.5, self.sample_num * 0.48 + 3.0)
        fig, ax = plt.subplots(figsize=(size, size))
        values = corr.to_numpy(dtype=float)
        vmin, vmax = self._adaptive_correlation_limits(values)
        image = ax.imshow(
            values,
            aspect="equal",
            interpolation="nearest",
            cmap="RdBu_r",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_xticks(range(self.sample_num))
        ax.set_yticks(range(self.sample_num))
        ax.set_xticklabels(self.sample_name, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(self.sample_name, fontsize=8)
        ax.set_title("Codon pausing score correlation")
        ax.tick_params(length=0)

        if self.sample_num <= HEATMAP_ANNOTATION_MAX_SAMPLES:
            color_span = max(vmax - vmin, np.finfo(float).eps)
            for row_idx in range(values.shape[0]):
                for col_idx in range(values.shape[1]):
                    value = values[row_idx, col_idx]
                    if not np.isfinite(value):
                        continue
                    color_fraction = (value - vmin) / color_span
                    text_color = "white" if color_fraction < 0.18 or color_fraction > 0.82 else "#222222"
                    ax.text(
                        col_idx,
                        row_idx,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=text_color,
                    )

        cbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.03)
        cbar.set_label(f"Pearson correlation ({vmin:.2f} to {vmax:.2f})")
        fig.tight_layout()

        pdf = self.output + "_pausing_corrplot.pdf"
        png = self.output + "_pausing_corrplot.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["pausing_corrplot_pdf"] = pdf
        self.output_files["pausing_corrplot_png"] = png

    def draw_codon_pausing_heatmap(self) -> None:
        """Draw combined total/valid codon pausing heatmaps."""
        relative = self.scale != "none"
        total_matrix = self._codon_matrix("total", relative=relative)
        valid_matrix = self._codon_matrix("valid", relative=relative)

        labels = (
            self.codon_pausing.set_index("Codon")
            .loc[total_matrix.index, :]
            .apply(lambda row: f"{row.name} [{row['Abbr']}]", axis=1)
            .tolist()
        )
        diverging = relative and self.scale == "zscore"
        cmap = PLOT_CMAP_DIVERGING if diverging else PLOT_CMAP_SEQUENTIAL
        vmin, vmax = self._heatmap_limits(total_matrix, valid_matrix, diverging)

        figure_height = max(9.0, len(total_matrix) * 0.22)
        figure_width = max(10.0, self.sample_num * 0.55 + 7.5)
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(figure_width, figure_height),
            gridspec_kw={"wspace": 0.28},
        )

        images = []
        for ax, matrix, title in zip(
            axes,
            (total_matrix, valid_matrix),
            ("All codon positions", "Positions with detected RPF"),
        ):
            values = matrix.to_numpy(dtype=float)
            image = ax.imshow(
                values,
                aspect="auto",
                interpolation="nearest",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )
            images.append(image)
            ax.set_title(title, fontsize=11)
            ax.set_xlabel("Sample")
            ax.set_xticks(range(self.sample_num))
            ax.set_xticklabels(self.sample_name, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.tick_params(length=0)

            # Add numeric pausing scores only when the sample count is small
            # enough to keep the heatmap readable.
            if self.sample_num <= HEATMAP_ANNOTATION_MAX_SAMPLES:
                color_midpoint = (vmin + vmax) / 2.0
                for row_idx in range(values.shape[0]):
                    for col_idx in range(values.shape[1]):
                        value = values[row_idx, col_idx]
                        if not np.isfinite(value):
                            continue
                        color_span = max(vmax - vmin, np.finfo(float).eps)
                        color_fraction = (value - vmin) / color_span
                        text_color = "white" if color_fraction < 0.18 or color_fraction > 0.82 else "#222222"
                        ax.text(
                            col_idx,
                            row_idx,
                            f"{value:.2f}",
                            ha="center",
                            va="center",
                            fontsize=5.5,
                            color=text_color,
                        )

        axes[0].set_ylabel("Codon [amino acid]")
        axes[1].set_ylabel("")
        label = "Relative pausing score" if relative else "Absolute pausing score"
        if not relative and self.plot_transform != "none":
            label += f" ({self.plot_transform} transformed)"
        cbar = fig.colorbar(images[-1], ax=axes.ravel().tolist(), shrink=0.55, pad=0.02)
        cbar.set_label(label)

        out_pdf = self.output + "_codon_pausing_heatmap.pdf"
        out_png = self.output + "_codon_pausing_heatmap.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["codon_heatmap_pdf"] = out_pdf
        self.output_files["codon_heatmap_png"] = out_png

    def draw_codon_rank_plot(self) -> None:
        """Draw sample-wise codon rank profiles with codons on the x-axis."""
        matrix = self._codon_matrix("valid", relative=False)
        long_df = matrix.reset_index().melt(
            id_vars="Codon", var_name="Sample", value_name="PausingScore"
        )
        long_df = long_df.merge(
            self.codon_pausing[["Codon", "Abbr"]], on="Codon", how="left"
        )

        ncols = max(1, int(self.rankplot_ncol))
        ncols = min(ncols, max(1, self.sample_num))
        nrows = int(np.ceil(self.sample_num / ncols))
        n_codons = int(len(matrix.index))
        tick_fontsize = 14 if ncols == 1 else 12
        label_width = 3 * 0.42 * tick_fontsize / 72.0 + 0.01
        panel_width = max(5.5, 2.0 + n_codons * label_width)
        panel_height = 3.8
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(panel_width * ncols, panel_height * nrows),
            squeeze=False,
        )

        for ax, sample in zip(axes.ravel(), self.sample_name):
            data = long_df.loc[long_df["Sample"] == sample, :].dropna()
            data = data.sort_values("PausingScore", ascending=True).reset_index(drop=True)
            x = np.arange(len(data))
            ax.scatter(
                x,
                data["PausingScore"],
                s=15,
                alpha=0.8,
                edgecolors="none",
                label="_nolegend_",
            )
            top = data.tail(min(6, len(data)))
            ax.scatter(
                top.index.to_numpy(),
                top["PausingScore"],
                s=24,
                alpha=0.95,
                color="#c0392b",
                edgecolors="none",
                zorder=3,
                label="Top 6 paused codons",
            )
            ax.set_xticks(x)
            ax.set_xticklabels(
                data["Codon"] + "[" + data["Abbr"].astype(str) + "]",
                rotation=90,
                ha="center",
                fontsize=tick_fontsize,
            )
            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlim(-0.5, len(data) - 0.5)
            ax.set_title(sample, fontsize=14)
            ax.set_xlabel("Codon [amino acid]", fontsize=14)
            ax.set_ylabel("Absolute valid pausing score", fontsize=14)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(axis="y", linewidth=0.4, alpha=0.25)

        for ax in axes.ravel()[self.sample_num :]:
            ax.set_axis_off()

        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                frameon=False,
                loc="upper center",
                ncol=len(handles),
                fontsize=12,
            )

        out_pdf = self.output + "_codon_pausing_rankplot.pdf"
        out_png = self.output + "_codon_pausing_rankplot.png"
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["codon_rankplot_pdf"] = out_pdf
        self.output_files["codon_rankplot_png"] = out_png

    # Legacy method names now route to the combined figure.
    def draw_codon_total_pausing_heat(self) -> None:
        self.draw_codon_pausing_heatmap()

    def draw_codon_valid_pausing_heat(self) -> None:
        return None

    # ------------------------------------------------------------------
    # Summary and workflow
    # ------------------------------------------------------------------

    def write_summary(self) -> None:
        """Write a machine-readable pausing analysis summary."""
        path = self.output + "_pausing.summary.json"
        self.output_files["summary_json"] = path
        summary = OrderedDict(
            [
                ("tool", "rpf_Pausing"),
                ("version", "0.2.8-dev.005"),
                ("input_rpf", os.path.abspath(self.rpf)),
                ("input_format", self.file_format),
                (
                    "transcript_filter",
                    os.path.abspath(self.transcript) if self.transcript else None,
                ),
                ("sample_count", self.sample_num),
                ("samples", self.sample_name),
                (
                    "parameters",
                    OrderedDict(
                        [
                            ("site", self.site),
                            ("frame", self.frame),
                            ("background", self.background),
                            ("thread", self.thread),
                            ("min", self.rpf_num),
                            ("tis", self.tis),
                            ("tts", self.tts),
                            ("normal", self.normal),
                            ("scale", self.scale),
                            ("exclude_stop", True),
                            ("individual", self.individual),
                            ("plot_transform", self.plot_transform),
                            ("rankplot_ncol", self.rankplot_ncol),
                            ("remove_outlier", self.remove_outlier),
                            ("outlier_iqr", self.outlier_iqr),
                            ("outlier_window", self.outlier_window),
                            ("outlier_local_fold", self.outlier_local_fold),
                        ]
                    ),
                ),
                ("sample_summary", self.sample_summary.to_dict(orient="records")),
                ("output_files", self.output_files),
            ]
        )
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

    def run(self) -> None:
        """Run the full pausing analysis workflow."""
        self.import_rpf()
        self.calculate_pausing()
        self.output_cds_pausing()
        self.output_cds_codon_pausing()
        self.output_sum_codon_pausing()
        self.output_all_pausing()
        self.draw_codon_pausing_heatmap()
        self.draw_codon_rank_plot()
        self.write_summary()
