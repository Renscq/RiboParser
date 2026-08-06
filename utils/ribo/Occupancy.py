#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.002
# Function: Calculate sample-specific codon occupancy from RPF density data.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Codon occupancy tables, outlier records, summary JSON, and figures.

"""Core functions for RiboParser codon occupancy analysis.

The module reads compact JSONL and legacy TXT density files through
``RPFs.RPFData``. High-expression transcripts and codon occupancy are calculated
independently for each sample.

For transcript ``g`` and codon position ``i`` in sample ``s``::

    normalized_density[g, i, s] = density[g, i, s] / mean_CDS_density[g, s]

Codon occupancy is the mean normalized density across occurrences of each codon.
Stop codons are always excluded from calculation and visualization.
"""

from __future__ import annotations

import json
import os
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed

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
HEATMAP_CMAP = "RdBu_r"
HEATMAP_ANNOTATION_MAX_SAMPLES = 12


class Occupancy(object):
    """Calculate sample-specific codon occupancy.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Occupancy``.
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
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts
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
        self.merged_rpf = None
        self.sample_name = []
        self.sample_num = 0
        self.total_rpf_num = None
        self.file_format = None

        # Derived results.
        self.sample_high_genes = OrderedDict()
        self.sample_gene_counts = pd.DataFrame()
        self.position_occupancy = None
        self.gene_codon_density = None
        self.codon_occupancy_table = None
        self.occupancy_corr = None
        self.outliers = pd.DataFrame()
        self.sample_summary = pd.DataFrame()
        self.output_files = OrderedDict()

        self.codon_dict, self.codon_annotation = RPFs.codon_table()
        for codon in ("TAA", "TAG", "TGA"):
            self.codon_dict.pop(codon, None)
        self.codon_annotation = self.codon_annotation.loc[
            self.codon_annotation["Abbr"] != "*", :
        ]

    # ------------------------------------------------------------------
    # Import and validation
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import RPF density and prepare sample-level codon counts."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf,
            sample_name=None,
            gene=self.transcript,
            tis=self.tis,
            tts=self.tts,
        )
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf_num = self.rpf_data.total_rpf_num.astype(float)
        self.file_format = self.rpf_data.file_format

        self.merged_rpf = self.rpf_data.get_frame(frame=self.frame)
        shift_num = RPFs.set_codon_shift(self.site)
        self.merged_rpf = RPFs.shift_site(
            self.merged_rpf,
            self.sample_name,
            shift_num,
        )

        self._validate_table()
        self.merged_rpf = self.merged_rpf.loc[
            self.merged_rpf["codon"].isin(self.codon_dict.keys()),
            BASE_COLUMNS + self.sample_name,
        ].copy()

        for sample in self.sample_name:
            self.merged_rpf[sample] = pd.to_numeric(
                self.merged_rpf[sample], errors="coerce"
            ).fillna(0.0).astype(float)

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

    def _validate_table(self) -> None:
        """Validate the imported codon-level density table."""
        if self.merged_rpf is None or self.merged_rpf.empty:
            raise ValueError("RPF density table is empty after import and filtering.")
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

        self.sample_gene_counts = cds.groupby("name", sort=False)[self.sample_name].sum()
        records = []
        for sample in self.sample_name:
            high_genes = self.sample_gene_counts.index[
                self.sample_gene_counts[sample] >= float(self.rpf_num)
            ]
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
    # Outlier detection
    # ------------------------------------------------------------------

    @staticmethod
    def _local_background(values: np.ndarray, window: int) -> np.ndarray:
        """Calculate a symmetric local mean excluding the focal codon."""
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
        np.divide(neighbor_sum, neighbor_count, out=background, where=neighbor_count > 0)
        return background

    @staticmethod
    def _robust_upper_cutoff(values: pd.Series, multiplier: float) -> float:
        """Return a conservative robust upper cutoff on the log1p scale."""
        positive = values.loc[values > 0].astype(float)
        if positive.size < 4:
            return float("inf")

        log_values = np.log1p(positive.to_numpy(dtype=float))
        q1, q3 = np.percentile(log_values, [25, 75])
        iqr = float(q3 - q1)
        median = float(np.median(log_values))
        mad = float(np.median(np.abs(log_values - median)))
        mad_sigma = 1.4826 * mad
        cutoffs = []
        if iqr > 0:
            cutoffs.append(float(q3) + float(multiplier) * iqr)
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
        """Detect isolated extreme raw-density pileups for one sample."""
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

        candidate_genes = set(candidate["name"].astype(str))
        local_background = pd.Series(np.nan, index=candidate.index, dtype=float)
        candidate_by_gene = candidate.groupby("name", sort=False).groups
        restricted = sample_table.loc[
            sample_table["name"].astype(str).isin(candidate_genes), :
        ]

        for gene_name, gene_df in restricted.groupby("name", sort=False):
            background = self._local_background(
                gene_df[sample].to_numpy(dtype=float),
                self.outlier_window,
            )
            lookup = pd.Series(background, index=gene_df.index)
            candidate_index = candidate_by_gene.get(gene_name, [])
            if len(candidate_index) > 0:
                local_background.loc[candidate_index] = lookup.loc[candidate_index].to_numpy()

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

    # ------------------------------------------------------------------
    # Occupancy calculation
    # ------------------------------------------------------------------

    def _calculate_sample_occupancy(
        self,
        sample: str,
    ) -> tuple[str, pd.DataFrame, pd.DataFrame, dict[str, int]]:
        """Calculate transcript-normalized codon density for one sample."""
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
        cds_mask = sample_table["region"].eq("cds")
        gene_mean = (
            sample_table.loc[cds_mask, ["name", sample]]
            .groupby("name", sort=False)[sample]
            .mean()
        )
        denominator = sample_table["name"].map(gene_mean).to_numpy(dtype=float)
        occupancy = np.full(values.shape, np.nan, dtype=float)
        valid = np.isfinite(values) & np.isfinite(denominator) & (denominator > 0)
        np.divide(values, denominator, out=occupancy, where=valid)

        result = sample_table.loc[:, BASE_COLUMNS].copy()
        result["Sample"] = sample
        result["RawDensity"] = values
        result["GeneMeanDensity"] = denominator
        result["Occupancy"] = occupancy

        summary = {
            "HighGeneCount": int(len(high_genes)),
            "OutlierCount": int(len(outlier_index)),
            "ValidOccupancyPositionCount": int(np.isfinite(occupancy).sum()),
        }
        return sample, result, outlier_df, summary

    def calculate_occupancy(self) -> None:
        """Calculate occupancy with optional sample-level multithreading."""
        if self.merged_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        worker_count = min(self.thread, max(1, self.sample_num))
        results = {}
        if worker_count == 1:
            for sample in self.sample_name:
                sample_name, table, outlier_df, summary = self._calculate_sample_occupancy(sample)
                results[sample_name] = (table, outlier_df, summary)
        else:
            print(
                f"Calculate codon occupancy with {worker_count} sample workers.",
                flush=True,
            )
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(self._calculate_sample_occupancy, sample): sample
                    for sample in self.sample_name
                }
                for future in as_completed(futures):
                    sample_name, table, outlier_df, summary = future.result()
                    results[sample_name] = (table, outlier_df, summary)

        tables = []
        outlier_tables = []
        summary_updates = []
        for sample in self.sample_name:
            table, outlier_df, summary = results[sample]
            tables.append(table)
            if not outlier_df.empty:
                outlier_tables.append(outlier_df)
            summary_updates.append({"Sample": sample, **summary})
            print(
                "Sample {sample}: high genes={genes:,}, valid positions={positions:,}, outliers={outliers:,}.".format(
                    sample=sample,
                    genes=summary["HighGeneCount"],
                    positions=summary["ValidOccupancyPositionCount"],
                    outliers=summary["OutlierCount"],
                ),
                flush=True,
            )

        self.position_occupancy = pd.concat(tables, axis=0, ignore_index=True)
        if outlier_tables:
            self.outliers = pd.concat(outlier_tables, axis=0, ignore_index=True)
        else:
            self.outliers = pd.DataFrame(
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

        update_df = pd.DataFrame.from_records(summary_updates)
        self.sample_summary = (
            self.sample_summary.drop(columns=["HighGeneCount"], errors="ignore")
            .merge(update_df, on="Sample", how="left")
        )

    # Backward-compatible method name.
    codon_occupancy = calculate_occupancy

    @staticmethod
    def _scale_matrix(matrix: pd.DataFrame, method: str) -> pd.DataFrame:
        """Scale each sample column independently."""
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
        raise ValueError(f"Unsupported scale method: {method}")

    def summarize_occupancy(self) -> None:
        """Summarize transcript-codon and codon-level occupancy."""
        if self.position_occupancy is None:
            raise ValueError("Codon occupancy has not been calculated yet.")

        cds = self.position_occupancy.loc[
            self.position_occupancy["region"] == "cds", :
        ].copy()
        cds["IsValidCodon"] = cds["RawDensity"].fillna(0.0) > 0

        self.gene_codon_density = (
            cds.groupby(["name", "codon", "Sample"], sort=False)
            .agg(
                CodonCount=("Occupancy", "size"),
                ValidCodonCount=("IsValidCodon", "sum"),
                RawDensitySum=("RawDensity", "sum"),
                OccupancySum=("Occupancy", "sum"),
                OccupancyMean=("Occupancy", "mean"),
            )
            .reset_index()
        )
        self.gene_codon_density["ValidCodonCount"] = (
            self.gene_codon_density["ValidCodonCount"].astype(int)
        )

        occurrence = cds.groupby(["codon", "Sample"], sort=False).agg(
            CodonCount=("Occupancy", "size"),
            ValidCodonCount=("IsValidCodon", "sum"),
            Density=("Occupancy", "sum"),
            AbsoluteOccupancy=("Occupancy", "mean"),
        )
        occurrence = occurrence.reset_index()

        absolute = occurrence.pivot(
            index="codon", columns="Sample", values="AbsoluteOccupancy"
        ).reindex(columns=self.sample_name)
        density = occurrence.pivot(
            index="codon", columns="Sample", values="Density"
        ).reindex(columns=self.sample_name)
        codon_count = occurrence.pivot(
            index="codon", columns="Sample", values="CodonCount"
        ).reindex(columns=self.sample_name)
        valid_count = occurrence.pivot(
            index="codon", columns="Sample", values="ValidCodonCount"
        ).reindex(columns=self.sample_name)
        relative = self._scale_matrix(absolute, self.scale)

        codon_index = self.codon_annotation.index.intersection(absolute.index)
        result = self.codon_annotation.reindex(codon_index).copy()
        for sample in self.sample_name:
            result[f"{sample}_codon_count"] = codon_count.reindex(codon_index)[sample].fillna(0).astype(int)
            result[f"{sample}_valid_codon"] = valid_count.reindex(codon_index)[sample].fillna(0).astype(int)
            result[f"{sample}_density"] = density.reindex(codon_index)[sample]
            result[f"{sample}_absolute_occupancy"] = absolute.reindex(codon_index)[sample]
            result[f"{sample}_relative_occupancy"] = relative.reindex(codon_index)[sample]

        result.index.name = "Codon"
        self.codon_occupancy_table = (
            result.reset_index()
            .sort_values(["Abbr", "Codon"], kind="stable")
            .reset_index(drop=True)
        )

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def output_tables(self) -> None:
        """Write standard and optional detailed occupancy tables."""
        self.summarize_occupancy()

        path = self.output + "_codon_occupancy.txt"
        self.codon_occupancy_table.to_csv(path, sep="\t", index=False)
        self.output_files["codon_occupancy_table"] = path

        if self.output_all:
            position_path = self.output + "_rpf_density.txt"
            self.position_occupancy.to_csv(position_path, sep="\t", index=False)
            self.output_files["position_occupancy_table"] = position_path

            gene_codon_path = self.output + "_codon_density.txt"
            self.gene_codon_density.to_csv(gene_codon_path, sep="\t", index=False)
            self.output_files["gene_codon_density_table"] = gene_codon_path

        if self.remove_outlier:
            outlier_path = self.output + "_occupancy.outliers.txt"
            self.outliers.to_csv(outlier_path, sep="\t", index=False)
            self.output_files["outlier_table"] = outlier_path

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _transform_plot_values(self, values: np.ndarray) -> np.ndarray:
        """Transform non-negative occupancy values for plotting only."""
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
        raise ValueError(f"Unsupported plot transform: {self.plot_transform}")

    def _occupancy_matrix(self, relative: bool) -> pd.DataFrame:
        """Return codon-by-sample occupancy matrix."""
        suffix = "relative_occupancy" if relative else "absolute_occupancy"
        columns = [f"{sample}_{suffix}" for sample in self.sample_name]
        matrix = self.codon_occupancy_table.set_index("Codon")[columns].copy()
        matrix.columns = self.sample_name
        if not relative:
            matrix.loc[:, :] = self._transform_plot_values(matrix.to_numpy(dtype=float))
        return matrix

    @staticmethod
    def _heatmap_limits(matrix: pd.DataFrame, diverging: bool) -> tuple[float, float]:
        """Return robust heatmap color limits."""
        values = matrix.to_numpy(dtype=float).ravel()
        values = values[np.isfinite(values)]
        if values.size == 0:
            return (0.0, 1.0)
        if diverging:
            bound = float(np.nanpercentile(np.abs(values), 98))
            bound = bound if bound > 0 else 1.0
            return (-bound, bound)
        upper = float(np.nanpercentile(values, 98))
        return (0.0, upper if upper > 0 else 1.0)

    def draw_occupancy_heatmap(self) -> None:
        """Draw a polished codon occupancy heatmap."""
        relative = self.scale != "none"
        matrix = self._occupancy_matrix(relative=relative)
        labels = (
            self.codon_occupancy_table.set_index("Codon")
            .loc[matrix.index, :]
            .apply(lambda row: f"{row.name} [{row['Abbr']}]", axis=1)
            .tolist()
        )
        diverging = relative and self.scale == "zscore"
        vmin, vmax = self._heatmap_limits(matrix, diverging)

        figure_height = max(9.0, len(matrix) * 0.22)
        figure_width = max(7.5, self.sample_num * 0.62 + 4.5)
        fig, ax = plt.subplots(figsize=(figure_width, figure_height))
        values = matrix.to_numpy(dtype=float)
        image = ax.imshow(
            values,
            aspect="auto",
            interpolation="nearest",
            cmap=HEATMAP_CMAP,
            vmin=vmin,
            vmax=vmax,
        )
        title = "Relative codon occupancy" if relative else "Absolute codon occupancy"
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Sample")
        ax.set_ylabel("Codon")
        ax.set_xticks(range(self.sample_num))
        ax.set_xticklabels(self.sample_name, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.tick_params(length=0)

        if self.sample_num <= HEATMAP_ANNOTATION_MAX_SAMPLES:
            span = max(vmax - vmin, np.finfo(float).eps)
            for row_idx in range(values.shape[0]):
                for col_idx in range(values.shape[1]):
                    value = values[row_idx, col_idx]
                    if not np.isfinite(value):
                        continue
                    fraction = (value - vmin) / span
                    color = "white" if fraction < 0.18 or fraction > 0.82 else "#222222"
                    ax.text(
                        col_idx,
                        row_idx,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=5.5,
                        color=color,
                    )

        cbar = fig.colorbar(image, ax=ax, shrink=0.68, pad=0.02)
        cbar.set_label(title)
        fig.tight_layout()

        pdf = self.output + "_occupancy_heatmap.pdf"
        png = self.output + "_occupancy_heatmap.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["occupancy_heatmap_pdf"] = pdf
        self.output_files["occupancy_heatmap_png"] = png

    @staticmethod
    def _adaptive_correlation_limits(values: np.ndarray) -> tuple[float, float]:
        """Return an adaptive correlation range emphasizing sample differences.

        Diagonal values are excluded when estimating the lower limit. A small
        robust margin is added below the lower off-diagonal correlations, while
        the upper limit remains 1.0. This preserves the interpretation of a
        correlation matrix but avoids compressing highly correlated samples into
        nearly identical colors.
        """
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

    def draw_occupancy_corr(self) -> None:
        """Draw sample correlation based on absolute codon occupancy."""
        matrix = self._occupancy_matrix(relative=False)
        corr = matrix.corr(method="pearson")
        self.occupancy_corr = corr

        txt = self.output + "_occupancy_corr.txt"
        corr.to_csv(txt, sep="\t", index=True)
        self.output_files["occupancy_correlation_table"] = txt

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
        ax.set_title("Codon occupancy correlation")
        ax.tick_params(length=0)

        if self.sample_num <= HEATMAP_ANNOTATION_MAX_SAMPLES:
            for row_idx in range(values.shape[0]):
                for col_idx in range(values.shape[1]):
                    value = values[row_idx, col_idx]
                    if np.isfinite(value):
                        color_fraction = (value - vmin) / (vmax - vmin)
                        color = (
                            "white"
                            if color_fraction < 0.18 or color_fraction > 0.82
                            else "#222222"
                        )
                        ax.text(
                            col_idx,
                            row_idx,
                            f"{value:.2f}",
                            ha="center",
                            va="center",
                            fontsize=7,
                            color=color,
                        )

        cbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.03)
        cbar.set_label(f"Pearson correlation ({vmin:.2f} to {vmax:.2f})")
        fig.tight_layout()

        pdf = self.output + "_occupancy_corrplot.pdf"
        png = self.output + "_occupancy_corrplot.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["occupancy_corrplot_pdf"] = pdf
        self.output_files["occupancy_corrplot_png"] = png

    def draw_occupancy_rankplot(self) -> None:
        """Draw sample-wise codon occupancy ranks with codons on the x-axis."""
        matrix = self._occupancy_matrix(relative=False)
        long_df = matrix.reset_index().melt(
            id_vars="Codon", var_name="Sample", value_name="OccupancyScore"
        )
        long_df = long_df.merge(
            self.codon_occupancy_table[["Codon", "Abbr"]], on="Codon", how="left"
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
            data = data.sort_values("OccupancyScore", ascending=True).reset_index(drop=True)
            x = np.arange(len(data))
            ax.scatter(
                x,
                data["OccupancyScore"],
                s=15,
                alpha=0.8,
                edgecolors="none",
                label="_nolegend_",
            )
            top = data.tail(min(6, len(data)))
            ax.scatter(
                top.index.to_numpy(),
                top["OccupancyScore"],
                s=24,
                alpha=0.95,
                color="#c0392b",
                edgecolors="none",
                zorder=3,
                label="Top 6 occupied codons",
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
            ax.set_ylabel("Absolute occupancy", fontsize=14)
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

        pdf = self.output + "_occupancy_rankplot.pdf"
        png = self.output + "_occupancy_rankplot.png"
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["occupancy_rankplot_pdf"] = pdf
        self.output_files["occupancy_rankplot_png"] = png

    # Backward-compatible wrappers.
    draw_occupancy_heat = draw_occupancy_heatmap
    draw_occupancy_relative_heat = draw_occupancy_heatmap
    draw_occupancy_line = draw_occupancy_rankplot

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def write_summary(self) -> None:
        """Write a machine-readable occupancy summary JSON."""
        path = self.output + "_occupancy.summary.json"
        self.output_files["summary_json"] = path
        summary = OrderedDict(
            [
                ("tool", "rpf_Occupancy"),
                ("version", "0.2.8-dev.002"),
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
                            ("min", self.rpf_num),
                            ("tis", self.tis),
                            ("tts", self.tts),
                            ("normal", self.normal),
                            ("thread", self.thread),
                            ("scale", self.scale),
                            ("plot_transform", self.plot_transform),
                            ("rankplot_ncol", self.rankplot_ncol),
                            ("exclude_stop", True),
                            ("remove_outlier", self.remove_outlier),
                            ("outlier_iqr", self.outlier_iqr),
                            ("outlier_window", self.outlier_window),
                            ("outlier_local_fold", self.outlier_local_fold),
                            ("output_all", self.output_all),
                        ]
                    ),
                ),
                ("sample_summary", self.sample_summary.to_dict(orient="records")),
                ("output_files", self.output_files),
            ]
        )
        with open(path, "w", encoding="utf-8") as out:
            json.dump(summary, out, ensure_ascii=False, indent=2)
            out.write("\n")

    def run(self) -> None:
        """Run the complete occupancy workflow."""
        self.import_rpf()
        self.calculate_occupancy()
        self.output_tables()
        self.draw_occupancy_corr()
        self.draw_occupancy_heatmap()
        self.draw_occupancy_rankplot()
        self.write_summary()
