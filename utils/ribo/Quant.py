#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.005
# Function: Provide core functions for RPF gene-level quantification.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Region-level RPF count, RPM, RPKM, TPM tables, QC plots, outlier tables, and summary JSON.

"""Core quantification utilities for RiboParser RPF density files.

This module contains the reusable implementation used by ``rpf_Quant``. The
command-line parser is intentionally kept in ``rpf_Quant.py`` so that the
quantification workflow can be reused by downstream scripts and tests.

The reader delegates density import to ``utils.ribo.RPFs.RPFData`` so both
compact JSONL files from updated ``rpf_Density`` / ``rpf_Merge`` and legacy TXT
codon-level density files are supported through the same workflow.
"""

from __future__ import annotations

import json
import os
from argparse import Namespace
from collections import OrderedDict
from typing import Any

from . import RPFs

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import seaborn as sns
except ImportError:  # pragma: no cover - optional plotting dependency
    sns = None

BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
REGION_CHOICES = ["cds", "5utr", "3utr", "all"]
REGION_ORDER = ["5utr", "cds", "3utr"]
REGION_LABELS = OrderedDict(
    [
        ("5utr", "5' UTR"),
        ("cds", "CDS"),
        ("3utr", "3' UTR"),
    ]
)
REGION_OUTPUT_NAMES = {
    "5utr": "utr5",
    "cds": "cds",
    "3utr": "utr3",
}
REGION_COLORS = OrderedDict(
    [
        ("5utr", "#E64B35"),
        ("cds", "#4DBBD5"),
        ("3utr", "#00A087"),
    ]
)
RPM_SCALE = 1_000_000.0
RPKM_SCALE = 1_000_000_000.0
TPM_SCALE = 1_000_000.0
HEATMAP_TOP_GENES = 1500
OUTLIER_IQR_MULTIPLIER = 8.0
OUTLIER_WINDOW = 5
OUTLIER_LOCAL_FOLD = 10.0


class Quant(object):
    """Quantify gene-level RPF abundance from codon-level density tables.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Quant``.
    """

    def __init__(self, args: Namespace):
        self.rpf_file = args.rpf
        self.output_prefix = args.output
        self.transcript = args.transcript

        self.frame = args.frame
        self.site = "P"
        self.tis = int(args.tis)
        self.tts = int(args.tts)
        self.region = args.region
        self.remove_outlier = bool(getattr(args, "remove_outlier", False))

        self.rpf_data = None
        self.raw_rpf = None
        self.merged_rpf = None
        self.total_rpf = None
        self.sample_name: list[str] = []
        self.sample_num = 0
        self.file_format = None

        self.region_results: OrderedDict[str, dict[str, pd.DataFrame]] = OrderedDict()
        self.region_summary = pd.DataFrame()
        self.outliers = pd.DataFrame()
        self.output_files: OrderedDict[str, str] = OrderedDict()

        # Backward-compatible attributes used by older downstream code.
        self.utr5_rpf = None
        self.utr5_rpm = None
        self.utr5_rpkm = None
        self.utr5_tpm = None
        self.cds_rpf = None
        self.cds_rpm = None
        self.cds_rpkm = None
        self.cds_tpm = None
        self.utr3_rpf = None
        self.utr3_rpm = None
        self.utr3_rpkm = None
        self.utr3_tpm = None

    # ------------------------------------------------------------------
    # Input helpers
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import RPF density data through the updated RPFData reader."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf_file,
            sample_name=None,
            gene=self.transcript,
            tis=None,
            tts=None,
        )
        self.raw_rpf = self.rpf_data.raw_rpf
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf = self.rpf_data.total_rpf_num.astype(float)
        self.file_format = self.rpf_data.file_format

        self.merged_rpf = self.rpf_data.shifted_frame(sites=self.site, frame=self.frame)
        self._validate_density_table(self.merged_rpf)
        self._coerce_sample_columns()

        print(
            "Imported RPF density for {sample_num} sample(s), format={fmt}, rows={rows:,}.".format(
                sample_num=self.sample_num,
                fmt=self.file_format,
                rows=len(self.merged_rpf),
            ),
            flush=True,
        )

    def read_rpf(self) -> None:
        """Backward-compatible alias of :meth:`import_rpf`."""
        self.import_rpf()

    def _validate_density_table(self, density_table: pd.DataFrame) -> None:
        """Validate required columns in the imported density table."""
        if density_table is None or density_table.empty:
            raise ValueError("RPF density table is empty after import.")

        missing_base = [column for column in BASE_COLUMNS if column not in density_table.columns]
        if missing_base:
            raise ValueError(
                "RPF density table is missing required base column(s): {columns}".format(
                    columns=", ".join(missing_base)
                )
            )

        missing_samples = [sample for sample in self.sample_name if sample not in density_table.columns]
        if missing_samples:
            raise ValueError(
                "RPF density table is missing sample column(s): {columns}".format(
                    columns=", ".join(missing_samples)
                )
            )

    def _coerce_sample_columns(self) -> None:
        """Convert sample columns to numeric values."""
        for sample in self.sample_name:
            self.merged_rpf[sample] = (
                pd.to_numeric(self.merged_rpf[sample], errors="coerce")
                .fillna(0.0)
                .astype("float64")
            )

    def _selected_regions(self) -> list[str]:
        """Return selected transcript regions in biological order."""
        if self.region == "all":
            return REGION_ORDER.copy()
        return [self.region]

    def _sample_total_vector(self) -> pd.Series:
        """Return total RPF count vector aligned to sample columns."""
        totals = self.total_rpf.reindex(self.sample_name).astype(float)
        return totals.replace(0.0, np.nan)

    # ------------------------------------------------------------------
    # Quantification helpers
    # ------------------------------------------------------------------

    def _region_density_table(self, region: str) -> pd.DataFrame:
        """Return codon-level density table for one region."""
        if self.merged_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        region_df = self.merged_rpf.loc[self.merged_rpf["region"] == region, BASE_COLUMNS + self.sample_name].copy()

        if region == "cds":
            if self.tis > 0:
                region_df = region_df.loc[region_df["from_tis"] >= self.tis, :]
            if self.tts > 0:
                region_df = region_df.loc[region_df["from_tts"] < -self.tts, :]

        clean_df, outlier_df = self._remove_outlier_points(region_df, region)
        if not outlier_df.empty:
            if self.outliers.empty:
                self.outliers = outlier_df.copy()
            else:
                self.outliers = pd.concat([self.outliers, outlier_df], axis=0, ignore_index=True)

        return clean_df

    @staticmethod
    def _robust_upper_cutoff(values: pd.Series, multiplier: float = OUTLIER_IQR_MULTIPLIER) -> float:
        """Return a dynamic robust upper cutoff on the log1p scale."""
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
            cutoffs.append(q3 + multiplier * iqr)
        if mad_sigma > 0:
            cutoffs.append(median + multiplier * mad_sigma)

        if not cutoffs:
            return float("inf")

        return float(np.expm1(max(cutoffs)))

    def _remove_outlier_points(self, region_df: pd.DataFrame, region: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Replace isolated sample-specific codon outliers with zero before quantification."""
        outlier_columns = [
            "Sample",
            "Region",
            "name",
            "now_nt",
            "from_tis",
            "from_tts",
            "codon",
            "RawDensity",
            "LocalBackground",
            "LocalFold",
            "OutlierCutoff",
            "OutlierMethod",
            "OutlierReason",
        ]

        if not self.remove_outlier or region_df.empty:
            return region_df, pd.DataFrame(columns=outlier_columns)

        clean_df = region_df.copy()
        work_df = region_df.copy()
        sort_columns = [column for column in ["name", "from_tis", "from_tts", "now_nt"] if column in work_df.columns]
        if sort_columns:
            work_df = work_df.sort_values(sort_columns, kind="mergesort")
        work_df["_PositionIndex"] = work_df.groupby("name", sort=False).cumcount()

        outlier_tables = []
        for sample in self.sample_name:
            values = pd.to_numeric(work_df[sample], errors="coerce").fillna(0.0).astype(float)
            cutoff = self._robust_upper_cutoff(values)
            if not np.isfinite(cutoff):
                continue

            candidate_mask = values > cutoff
            if not bool(candidate_mask.any()):
                continue

            density_lookup = (
                work_df.loc[:, ["name", "_PositionIndex", sample]]
                .set_index(["name", "_PositionIndex"])[sample]
                .astype(float)
            )

            candidate_df = work_df.loc[candidate_mask, :].copy()
            local_backgrounds = []
            local_folds = []
            remove_indices = []

            for row in candidate_df.itertuples():
                row_index = row.Index
                transcript = candidate_df.at[row_index, "name"]
                position = int(candidate_df.at[row_index, "_PositionIndex"])
                raw_density = float(candidate_df.at[row_index, sample])

                neighbor_values = []
                for offset in range(1, OUTLIER_WINDOW + 1):
                    left_value = density_lookup.get((transcript, position - offset), np.nan)
                    right_value = density_lookup.get((transcript, position + offset), np.nan)
                    if pd.notna(left_value):
                        neighbor_values.append(float(left_value))
                    if pd.notna(right_value):
                        neighbor_values.append(float(right_value))

                local_background = float(np.median(neighbor_values)) if neighbor_values else 0.0
                local_fold = (raw_density + 1.0) / (local_background + 1.0)
                local_backgrounds.append(local_background)
                local_folds.append(local_fold)

                if local_fold >= OUTLIER_LOCAL_FOLD:
                    remove_indices.append(row.Index)

            if not remove_indices:
                continue

            candidate_df["LocalBackground"] = local_backgrounds
            candidate_df["LocalFold"] = local_folds
            outlier_df = candidate_df.loc[remove_indices, :].copy()
            outlier_df["Sample"] = sample
            outlier_df["Region"] = region
            outlier_df["RawDensity"] = outlier_df[sample].astype(float)
            outlier_df["OutlierCutoff"] = cutoff
            outlier_df["OutlierMethod"] = "fast_log1p_robust_cutoff_and_local_peak"
            outlier_df["OutlierReason"] = (
                "RawDensity > dynamic robust cutoff and LocalFold >= "
                + str(float(OUTLIER_LOCAL_FOLD))
            )
            outlier_tables.append(outlier_df.reindex(columns=outlier_columns))
            clean_df.loc[remove_indices, sample] = 0.0

        if not outlier_tables:
            return clean_df, pd.DataFrame(columns=outlier_columns)

        return clean_df, pd.concat(outlier_tables, axis=0, ignore_index=True)

    @staticmethod
    def _region_lengths(region_df: pd.DataFrame) -> pd.Series:
        """Calculate quantified region length in nucleotides for each transcript."""
        length = region_df.groupby("name", sort=False).size().astype(float) * 3.0
        length.name = "length_nt"
        return length

    def _count_region(self, region_df: pd.DataFrame) -> pd.DataFrame:
        """Summarize codon-level density into transcript-level counts."""
        if region_df.empty:
            return pd.DataFrame(columns=self.sample_name, dtype="float64")

        count = region_df.groupby("name", sort=False)[self.sample_name].sum()
        count = count.astype(float)
        count.index.name = "name"
        return count

    def _calculate_rpm(self, count: pd.DataFrame) -> pd.DataFrame:
        """Calculate RPM from transcript-level RPF counts."""
        totals = self._sample_total_vector()
        rpm = count.div(totals, axis=1).mul(RPM_SCALE).replace([np.inf, -np.inf], np.nan).fillna(0.0)
        rpm.index.name = "name"
        return rpm

    def _calculate_rpkm(self, count: pd.DataFrame, length_nt: pd.Series) -> pd.DataFrame:
        """Calculate RPKM from transcript-level RPF counts."""
        totals = self._sample_total_vector()
        length = length_nt.replace(0.0, np.nan)
        rpkm = count.mul(RPKM_SCALE).div(length, axis=0).div(totals, axis=1)
        rpkm = rpkm.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        rpkm.index.name = "name"
        return rpkm

    def _calculate_tpm(self, count: pd.DataFrame, length_nt: pd.Series) -> pd.DataFrame:
        """Calculate TPM from transcript-level RPF counts."""
        length = length_nt.replace(0.0, np.nan)
        rpk = count.mul(1000.0).div(length, axis=0)
        rpk_sum = rpk.sum(axis=0).replace(0.0, np.nan)
        tpm = rpk.div(rpk_sum, axis=1).mul(TPM_SCALE)
        tpm = tpm.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        tpm.index.name = "name"
        return tpm

    @staticmethod
    def _rename_metric_columns(matrix: pd.DataFrame, region: str, metric: str) -> pd.DataFrame:
        """Return an output matrix with historical sample_region_metric column names."""
        output = matrix.copy()
        output.columns = ["{sample}_{region}_{metric}".format(sample=sample, region=region, metric=metric) for sample in output.columns]
        return output

    def _output_metric_table(self, matrix: pd.DataFrame, region: str, metric: str) -> str:
        """Write one quantification metric table."""
        out_region = REGION_OUTPUT_NAMES[region]
        out_file = "{prefix}_{region}_{metric}_quant.txt".format(
            prefix=self.output_prefix,
            region=out_region,
            metric=metric,
        )

        output = self._rename_metric_columns(matrix, region, metric)
        output.index.name = "name"

        if metric == "rpf":
            output = output.round(0).astype("int64")
        else:
            output = output.astype(float).round(2)

        output.to_csv(out_file, sep="\t", index=True)
        self.output_files["{region}_{metric}".format(region=out_region, metric=metric)] = out_file
        return out_file

    def quantify_region(self, region: str) -> dict[str, pd.DataFrame]:
        """Quantify RPF abundance for one transcript region."""
        if region not in REGION_ORDER:
            raise ValueError("Unsupported region: {region}".format(region=region))

        print("Quantify RPFs in {region}.".format(region=REGION_LABELS[region]), flush=True)
        region_df = self._region_density_table(region)
        count = self._count_region(region_df)
        length_nt = self._region_lengths(region_df)

        if count.empty:
            rpm = pd.DataFrame(columns=self.sample_name, dtype="float64")
            rpkm = pd.DataFrame(columns=self.sample_name, dtype="float64")
            tpm = pd.DataFrame(columns=self.sample_name, dtype="float64")
        else:
            rpm = self._calculate_rpm(count)
            rpkm = self._calculate_rpkm(count, length_nt)
            tpm = self._calculate_tpm(count, length_nt)

        result = {
            "count": count,
            "rpm": rpm,
            "rpkm": rpkm,
            "tpm": tpm,
            "length_nt": length_nt.to_frame(),
        }
        self.region_results[region] = result

        self._assign_backward_attributes(region, result)
        self._output_metric_table(count, region, "rpf")
        self._output_metric_table(rpm, region, "rpm")
        self._output_metric_table(rpkm, region, "rpkm")
        self._output_metric_table(tpm, region, "tpm")
        self._output_length_table(length_nt, region)

        return result

    def _output_length_table(self, length_nt: pd.Series, region: str) -> str:
        """Write quantified region length table."""
        out_region = REGION_OUTPUT_NAMES[region]
        out_file = "{prefix}_{region}_length.txt".format(prefix=self.output_prefix, region=out_region)
        length_nt.to_frame().to_csv(out_file, sep="\t", index=True)
        self.output_files["{region}_length".format(region=out_region)] = out_file
        return out_file

    def _assign_backward_attributes(self, region: str, result: dict[str, pd.DataFrame]) -> None:
        """Populate historical attributes for compatibility."""
        if region == "5utr":
            self.utr5_rpf = result["count"]
            self.utr5_rpm = result["rpm"]
            self.utr5_rpkm = result["rpkm"]
            self.utr5_tpm = result["tpm"]
        elif region == "cds":
            self.cds_rpf = result["count"]
            self.cds_rpm = result["rpm"]
            self.cds_rpkm = result["rpkm"]
            self.cds_tpm = result["tpm"]
        elif region == "3utr":
            self.utr3_rpf = result["count"]
            self.utr3_rpm = result["rpm"]
            self.utr3_rpkm = result["rpkm"]
            self.utr3_tpm = result["tpm"]

    def quantify_regions(self) -> None:
        """Quantify all selected regions."""
        summary_records = []

        for region in self._selected_regions():
            result = self.quantify_region(region)
            count = result["count"]
            length_nt = result["length_nt"]
            summary_records.append(
                {
                    "Region": region,
                    "RegionLabel": REGION_LABELS[region],
                    "TranscriptCount": int(count.shape[0]),
                    "MeanLengthNt": float(length_nt["length_nt"].mean()) if not length_nt.empty else 0.0,
                    "MedianLengthNt": float(length_nt["length_nt"].median()) if not length_nt.empty else 0.0,
                    "TotalRPF": float(count.sum().sum()) if not count.empty else 0.0,
                }
            )

        self.region_summary = pd.DataFrame.from_records(summary_records)
        out_file = self.output_prefix + "_quant.region_summary.txt"
        self.region_summary.to_csv(out_file, sep="\t", index=False)
        self.output_files["region_summary"] = out_file

        if self.remove_outlier:
            outlier_file = self.output_prefix + "_quant.outliers.txt"
            if self.outliers.empty:
                self.outliers = pd.DataFrame(
                    columns=[
                        "Sample",
                        "Region",
                        "name",
                        "now_nt",
                        "from_tis",
                        "from_tts",
                        "codon",
                        "RawDensity",
                        "LocalBackground",
                        "LocalFold",
                        "OutlierCutoff",
                        "OutlierMethod",
                        "OutlierReason",
                    ]
                )
            self.outliers.to_csv(outlier_file, sep="\t", index=False)
            self.output_files["outlier_table"] = outlier_file
            print(
                "Removed {count:,} sample-specific codon outlier(s).".format(count=len(self.outliers)),
                flush=True,
            )

    def quant_region(self) -> None:
        """Backward-compatible alias of :meth:`quantify_regions`."""
        self.quantify_regions()

    def output_total_rpf(self) -> None:
        """Write sample total RPF counts."""
        out_file = self.output_prefix + "_total.txt"
        total_rpf = self.total_rpf.copy()
        total_rpf.index.name = "sample"
        total_rpf.name = "rpf_count"
        total_rpf.to_csv(out_file, sep="\t", index=True)
        self.output_files["total_rpf"] = out_file

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_log2(matrix: pd.DataFrame) -> pd.DataFrame:
        """Return log2(x + 1) transformed matrix."""
        return np.log2(matrix.astype(float).clip(lower=0.0) + 1.0)

    @staticmethod
    def _row_zscore(matrix: pd.DataFrame) -> pd.DataFrame:
        """Return row-wise z-score matrix."""
        centered = matrix.sub(matrix.mean(axis=1), axis=0)
        scaled = centered.div(matrix.std(axis=1).replace(0.0, np.nan), axis=0)
        return scaled.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    @staticmethod
    def _auto_heatmap_limits(matrix: pd.DataFrame) -> tuple[float, float]:
        """Return robust symmetric color limits for a z-score heatmap."""
        values = matrix.to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return -1.0, 1.0
        limit = float(np.nanpercentile(np.abs(values), 98))
        if not np.isfinite(limit) or limit <= 0:
            limit = 1.0
        return -limit, limit

    def _cds_rpm_for_plots(self) -> pd.DataFrame | None:
        """Return CDS RPM matrix used by QC plots."""
        if "cds" not in self.region_results:
            return None
        rpm = self.region_results["cds"]["rpm"]
        if rpm is None or rpm.empty:
            return None
        return rpm

    def draw_rpf_barplot(self) -> None:
        """Draw a polished stacked barplot of RPF proportions across transcript regions."""
        if self.merged_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        out_pdf = self.output_prefix + "_region_rpf_proportion_bar_plot.pdf"
        out_png = self.output_prefix + "_region_rpf_proportion_bar_plot.png"
        out_txt = self.output_prefix + "_region_rpf_proportion.txt"

        region_sum = self.merged_rpf.groupby("region", sort=False)[self.sample_name].sum()
        region_sum = region_sum.reindex(REGION_ORDER).fillna(0.0)
        region_percent = region_sum.div(region_sum.sum(axis=0).replace(0.0, np.nan), axis=1).mul(100.0).fillna(0.0)
        region_percent.to_csv(out_txt, sep="\t", index=True)

        plot_df = region_percent.T
        figure_width = max(4.8, min(6.8, len(self.sample_name) * 0.24 + 2.2))
        fig, ax = plt.subplots(figsize=(figure_width, 7))
        bottom = np.zeros(len(plot_df.index), dtype=float)
        x = np.arange(len(plot_df.index))

        for region in REGION_ORDER:
            values = plot_df[region].to_numpy(dtype=float) if region in plot_df.columns else np.zeros(len(plot_df.index))
            ax.bar(
                x,
                values,
                bottom=bottom,
                width=0.82,
                label=REGION_LABELS[region],
                color=REGION_COLORS[region],
                edgecolor="white",
                linewidth=0.35,
            )
            bottom += values

        ax.set_ylabel("RPF proportion (%)")
        ax.set_xlabel("")
        ax.set_ylim(0, 100)
        ax.set_xticks(x)
        ax.set_xticklabels(plot_df.index.tolist(), rotation=90, ha="center", va="top")
        ax.tick_params(axis="x", pad=6, length=0)
        ax.tick_params(axis="y", length=3)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#BBBBBB")
        ax.spines["bottom"].set_color("#BBBBBB")
        ax.legend(
            title="Region",
            frameon=False,
            ncol=len(REGION_ORDER),
            loc="lower center",
            bbox_to_anchor=(0.5, 1.02),
            borderaxespad=0.0,
            columnspacing=1.2,
            handlelength=1.4,
        )
        fig.subplots_adjust(bottom=0.36, left=0.12, right=0.98, top=0.78)
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        self.output_files["region_bar_pdf"] = out_pdf
        self.output_files["region_bar_png"] = out_png
        self.output_files["region_proportion_table"] = out_txt

    def draw_rpf_cdfplot(self) -> None:
        """Draw empirical CDF curves of log2 CDS RPM values for detected genes.

        Each sample's CDF is computed from the genes with RPM > 0 in that
        sample only, so log2 transformation never produces -inf.
        """
        cds_rpm = self._cds_rpm_for_plots()
        if cds_rpm is None:
            print("Skip eCDF plot because CDS quantification is not available.", flush=True)
            return

        out_pdf = self.output_prefix + "_cds_rpm_cdf_plot.pdf"
        out_png = self.output_prefix + "_cds_rpm_cdf_plot.png"

        detected = cds_rpm.loc[(cds_rpm > 0).any(axis=1), :]
        if detected.empty:
            print("Skip eCDF plot because all CDS RPM values are zero.", flush=True)
            return

        figure_width = max(5.4, min(7.0, len(self.sample_name) * 0.18 + 3.8))
        fig, ax = plt.subplots(figsize=(figure_width, 7))
        cmap = plt.get_cmap("tab20")

        for idx, sample in enumerate(self.sample_name):
            if sample not in detected.columns:
                continue
            sample_values = detected[sample].dropna().to_numpy(dtype=float)
            sample_values = sample_values[sample_values > 0.0]
            if sample_values.size == 0:
                continue
            values = np.sort(np.log2(sample_values))
            y = np.arange(1, values.size + 1, dtype=float) / float(values.size)
            ax.step(values, y, where="post", label=sample, linewidth=1.35, color=cmap(idx % 20), alpha=0.95)

        ax.set_xlabel("log2(CDS RPM)")
        ax.set_ylabel("Empirical cumulative fraction")
        ax.set_ylim(0, 1.01)
        ax.grid(color="#E5E5E5", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#BBBBBB")
        ax.spines["bottom"].set_color("#BBBBBB")
        if len(self.sample_name) <= 24:
            if len(self.sample_name) == 1:
                legend_cols = 1
            elif len(self.sample_name) <= 6:
                legend_cols = 2
            elif len(self.sample_name) < 13:
                legend_cols = 3
            else:
                legend_cols = 4
            ax.legend(
                frameon=False,
                fontsize=7.5,
                ncol=legend_cols,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                borderaxespad=0.0,
                handlelength=1.3,
                columnspacing=0.9,
            )
            fig.subplots_adjust(bottom=0.26, left=0.11, right=0.98, top=0.94)
        else:
            fig.subplots_adjust(bottom=0.14, left=0.11, right=0.98, top=0.94)
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        self.output_files["cdf_pdf"] = out_pdf
        self.output_files["cdf_png"] = out_png

    @staticmethod
    def _place_non_overlapping_labels(
        fig: plt.Figure,
        ax: plt.Axes,
        anchors: list[tuple[float, float]],
        labels: list[str],
        fontsize: float = 7.5,
        max_iter: int = 80,
        step_px: float = 2.0,
        max_offset_px: float = 70.0,
    ) -> None:
        """Place text labels next to scatter anchors and push overlaps apart.

        Positions are adjusted in display (pixel) coordinates and mapped back
        to data coordinates so the layout is independent of figure size. Each
        label is kept within ``max_offset_px`` pixels of its anchor point.
        """
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = ax.transData.inverted()
        texts: list[Any] = []
        anchor_win: list[tuple[float, float]] = []
        for (x, y), label in zip(anchors, labels):
            texts.append(ax.text(x, y, "  " + str(label), fontsize=fontsize, va="center", zorder=4))
            anchor_win.append(tuple(ax.transData.transform((float(x), float(y)))))

        for _ in range(max_iter):
            moved = False
            for i in range(len(texts)):
                for j in range(i + 1, len(texts)):
                    bbox_i = texts[i].get_window_extent(renderer=renderer)
                    bbox_j = texts[j].get_window_extent(renderer=renderer)
                    if not bbox_i.overlaps(bbox_j):
                        continue
                    center_x_i = (bbox_i.x0 + bbox_i.x1) / 2.0
                    center_y_i = (bbox_i.y0 + bbox_i.y1) / 2.0
                    center_x_j = (bbox_j.x0 + bbox_j.x1) / 2.0
                    center_y_j = (bbox_j.y0 + bbox_j.y1) / 2.0
                    delta_x = center_x_j - center_x_i
                    delta_y = center_y_j - center_y_i
                    distance = float(np.hypot(delta_x, delta_y))
                    if distance < 1e-6:
                        delta_x, delta_y = step_px, step_px
                        distance = float(np.hypot(delta_x, delta_y))
                    unit_x = delta_x / distance
                    unit_y = delta_y / distance
                    for index, sign in ((i, -1.0), (j, 1.0)):
                        current = texts[index].get_position()
                        win = ax.transData.transform(current)
                        new_win_x = win[0] + sign * unit_x * step_px
                        new_win_y = win[1] + sign * unit_y * step_px
                        offset = float(np.hypot(new_win_x - anchor_win[index][0], new_win_y - anchor_win[index][1]))
                        if offset > max_offset_px:
                            ratio = max_offset_px / offset
                            new_win_x = anchor_win[index][0] + (new_win_x - anchor_win[index][0]) * ratio
                            new_win_y = anchor_win[index][1] + (new_win_y - anchor_win[index][1]) * ratio
                        texts[index].set_position(inv.transform((new_win_x, new_win_y)))
                    moved = True
            if not moved:
                break

    @staticmethod
    def _calculate_pca_by_svd(matrix: pd.DataFrame, n_components: int = 2) -> tuple[pd.DataFrame, np.ndarray]:
        """Calculate PCA scores and explained variance by centered SVD."""
        if matrix.shape[0] < 2 or matrix.shape[1] < 2:
            raise ValueError("At least two samples and two features are required for PCA.")

        values = matrix.to_numpy(dtype=float)
        values = values - np.nanmean(values, axis=0, keepdims=True)
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)

        feature_sd = values.std(axis=0)
        keep = feature_sd > 0
        if keep.sum() < 2:
            raise ValueError("At least two non-constant features are required for PCA.")
        values = values[:, keep]

        u_matrix, singular_values, _ = np.linalg.svd(values, full_matrices=False)
        component_count = min(n_components, len(singular_values))
        scores = u_matrix[:, :component_count] * singular_values[:component_count]

        variance = singular_values ** 2
        total_variance = variance.sum()
        if total_variance <= 0:
            variance_ratio = np.zeros(component_count, dtype=float)
        else:
            variance_ratio = variance[:component_count] / total_variance * 100.0

        if component_count < n_components:
            pad_width = n_components - component_count
            scores = np.pad(scores, ((0, 0), (0, pad_width)), mode="constant")
            variance_ratio = np.pad(variance_ratio, (0, pad_width), mode="constant")

        score_df = pd.DataFrame(scores[:, :n_components], index=matrix.index, columns=["PC{i}".format(i=i + 1) for i in range(n_components)])
        return score_df, variance_ratio[:n_components]

    def draw_rpf_pcaplot(self) -> None:
        """Draw PCA plot from log2 CDS RPM values using a pure NumPy SVD implementation."""
        cds_rpm = self._cds_rpm_for_plots()
        if cds_rpm is None:
            print("Skip PCA plot because CDS quantification is not available.", flush=True)
            return
        if self.sample_num < 2 or cds_rpm.shape[0] < 2:
            print("Skip PCA plot because at least two samples and two genes are required.", flush=True)
            return

        out_pdf = self.output_prefix + "_cds_rpm_pca_plot.pdf"
        out_png = self.output_prefix + "_cds_rpm_pca_plot.png"
        out_txt = self.output_prefix + "_cds_rpm_pca.txt"

        data_log = self._safe_log2(cds_rpm.loc[(cds_rpm > 0).any(axis=1), :]).T
        try:
            pca_df, variance_ratio = self._calculate_pca_by_svd(data_log, n_components=2)
        except ValueError as exc:
            print("Skip PCA plot because {error}".format(error=exc), flush=True)
            return

        pca_df.index.name = "Sample"
        pca_df["PC1Variance"] = float(variance_ratio[0])
        pca_df["PC2Variance"] = float(variance_ratio[1])
        pca_df.reset_index().to_csv(out_txt, sep="\t", index=False)

        figure_width = max(5.0, min(6.2, self.sample_num * 0.12 + 4.0))
        fig, ax = plt.subplots(figsize=(figure_width, 5.5))
        cmap = plt.get_cmap("tab20")
        anchor_points: list[tuple[float, float]] = []
        sample_labels: list[str] = []
        for idx, sample in enumerate(pca_df.index):
            ax.scatter(
                pca_df.loc[sample, "PC1"],
                pca_df.loc[sample, "PC2"],
                s=78,
                edgecolor="black",
                linewidth=0.45,
                color=cmap(idx % 20),
                label=sample,
                zorder=3,
            )
            anchor_points.append((float(pca_df.loc[sample, "PC1"]), float(pca_df.loc[sample, "PC2"])))
            sample_labels.append(str(sample))

        ax.axhline(0, color="#D0D0D0", linewidth=0.7, zorder=0)
        ax.axvline(0, color="#D0D0D0", linewidth=0.7, zorder=0)
        ax.grid(color="#EAEAEA", linewidth=0.55, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlabel("PC1 ({value:.2f}%)".format(value=float(variance_ratio[0])))
        ax.set_ylabel("PC2 ({value:.2f}%)".format(value=float(variance_ratio[1])))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#BBBBBB")
        ax.spines["bottom"].set_color("#BBBBBB")
        fig.subplots_adjust(bottom=0.14, left=0.14, right=0.96, top=0.93)
        if self.sample_num <= 20:
            self._place_non_overlapping_labels(fig, ax, anchor_points, sample_labels, fontsize=7.5)
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        self.output_files["pca_table"] = out_txt
        self.output_files["pca_pdf"] = out_pdf
        self.output_files["pca_png"] = out_png

    def _order_heatmap_rows(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """Order heatmap rows without SciPy dependency.

        The ordering uses the first two sample-pattern axes from SVD on the
        row-scaled matrix. This gives a stable pattern-oriented order that is
        close to a clustered heatmap display, while avoiding SciPy, which may
        not be available in minimal RiboParser environments.
        """
        if matrix.shape[0] < 2 or matrix.shape[1] < 2:
            return matrix

        values = matrix.to_numpy(dtype=float)
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
        if not np.any(values):
            return matrix

        try:
            centered = values - values.mean(axis=0, keepdims=True)
            u_matrix, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
            score1 = u_matrix[:, 0] * singular_values[0]
            if len(singular_values) > 1:
                score2 = u_matrix[:, 1] * singular_values[1]
            else:
                score2 = np.zeros_like(score1)
            peak_index = np.argmax(values, axis=1)
            order = np.lexsort((score2, peak_index, score1))
            return matrix.iloc[order, :]
        except Exception as exc:  # pragma: no cover - defensive plotting fallback
            print("Warning: heatmap row ordering failed and original order will be used: {error}".format(error=exc), flush=True)
            return matrix

    def draw_rpf_heatmap(self) -> None:
        """Draw an R-like heatmap from row-scaled CDS RPM values without SciPy."""
        cds_rpm = self._cds_rpm_for_plots()
        if cds_rpm is None:
            print("Skip heatmap because CDS quantification is not available.", flush=True)
            return

        out_pdf = self.output_prefix + "_cds_rpm_heatmap.pdf"
        out_png = self.output_prefix + "_cds_rpm_heatmap.png"

        rpm_log = self._safe_log2(cds_rpm.loc[(cds_rpm > 0).any(axis=1), :])
        if rpm_log.empty:
            print("Skip heatmap because all CDS RPM values are zero.", flush=True)
            return

        if rpm_log.shape[0] > HEATMAP_TOP_GENES:
            keep_gene = rpm_log.var(axis=1).sort_values(ascending=False).head(HEATMAP_TOP_GENES).index
            rpm_log = rpm_log.loc[keep_gene, :]

        rpm_scaled = self._row_zscore(rpm_log)
        rpm_scaled = rpm_scaled.loc[(rpm_scaled.abs().sum(axis=1) > 0), :]
        if rpm_scaled.empty:
            print("Skip heatmap because all row-scaled values are zero.", flush=True)
            return

        vmin, vmax = self._auto_heatmap_limits(rpm_scaled)
        figure_width = max(5.0, min(7.0, self.sample_num * 0.20 + 2.6))
        figure_height = min(max(6.5, rpm_scaled.shape[0] * 0.0012 + 2.4), 10.5)

        rpm_scaled = self._order_heatmap_rows(rpm_scaled)

        fig = plt.figure(figsize=(figure_width, figure_height))
        grid = fig.add_gridspec(
            nrows=1,
            ncols=2,
            width_ratios=[1.0, 0.035],
            wspace=0.06,
        )
        ax = fig.add_subplot(grid[0, 0])
        cax = fig.add_subplot(grid[0, 1])

        image = ax.imshow(
            rpm_scaled.to_numpy(dtype=float),
            aspect="auto",
            interpolation="nearest",
            cmap="RdBu_r",
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )

        ax.set_title("CDS expression pattern", fontsize=12, pad=10)
        ax.set_ylabel("Gene")
        ax.set_yticks([])
        ax.set_xticks(np.arange(len(rpm_scaled.columns)))
        ax.set_xticklabels(rpm_scaled.columns.tolist(), rotation=90, ha="center", va="top")
        ax.tick_params(axis="x", pad=7, length=0)
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

        cbar = fig.colorbar(image, cax=cax, orientation="vertical")
        cbar.set_label("Row z-score", fontsize=9)
        cbar.ax.tick_params(labelsize=8, length=2)
        for spine in cax.spines.values():
            spine.set_visible(False)

        fig.subplots_adjust(left=0.10, right=0.90, bottom=0.30, top=0.90)
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        self.output_files["heatmap_pdf"] = out_pdf
        self.output_files["heatmap_png"] = out_png

    def draw_rpf_heatmap2(self) -> None:
        """Backward-compatible alias of :meth:`draw_rpf_heatmap`."""
        self.draw_rpf_heatmap()

    # ------------------------------------------------------------------
    # Summary and workflow
    # ------------------------------------------------------------------

    def write_summary(self) -> None:
        """Write a machine-readable quantification summary JSON."""
        out_file = self.output_prefix + "_quant.summary.json"
        self.output_files["summary_json"] = out_file

        summary = OrderedDict(
            [
                ("tool", "rpf_Quant"),
                ("version", "0.2.8-dev.005"),
                ("input_rpf", os.path.abspath(self.rpf_file)),
                ("input_format", self.file_format),
                ("transcript_filter", os.path.abspath(self.transcript) if self.transcript else None),
                ("sample_count", self.sample_num),
                ("samples", self.sample_name),
                (
                    "parameters",
                    OrderedDict(
                        [
                            ("frame", self.frame),
                            ("site", self.site),
                            ("tis", self.tis),
                            ("tts", self.tts),
                            ("region", self.region),
                            ("remove_outlier", self.remove_outlier),
                            ("outlier_iqr_multiplier", OUTLIER_IQR_MULTIPLIER),
                            ("outlier_window", OUTLIER_WINDOW),
                            ("outlier_local_fold", OUTLIER_LOCAL_FOLD),
                        ]
                    ),
                ),
                ("region_summary", self.region_summary.to_dict(orient="records") if not self.region_summary.empty else []),
                ("outlier_count", int(len(self.outliers)) if self.remove_outlier else 0),
                ("output_files", self.output_files),
            ]
        )

        with open(out_file, "w", encoding="utf-8") as out:
            json.dump(summary, out, ensure_ascii=False, indent=2)
            out.write("\n")

    def run(self) -> None:
        """Run the full quantification workflow."""
        self.import_rpf()
        self.quantify_regions()
        self.output_total_rpf()
        self.draw_rpf_barplot()
        self.draw_rpf_cdfplot()
        self.draw_rpf_pcaplot()
        self.draw_rpf_heatmap()
        self.write_summary()
