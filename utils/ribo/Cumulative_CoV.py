#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.002
# Function: Calculate sample-specific cumulative coefficient of variation along CDS regions.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Cumulative CoV tables, summary JSON, outlier table, and meta summaries and per-transcript figures.

"""Core functions for cumulative coefficient-of-variation analysis."""

from __future__ import annotations

import json
import math
import os
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import RPFs

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
RPM_SCALE = 1_000_000.0
PLOT_DPI = 300


class CumulativeCoV:
    """Calculate cumulative CDS CoV independently for each sample."""

    def __init__(self, args):
        # Input and output.
        self.rpf = args.rpf
        self.output = args.output
        self.transcript = args.list
        self.output_all = args.all

        # Density filtering.
        self.site = args.site
        self.frame = args.frame
        self.normal = args.normal
        self.tis = args.tis
        self.tts = args.tts
        self.rpf_num = args.min
        self.min_positions = args.min_positions
        self.trim = args.trim
        self.resolution = args.resolution
        self.thread = max(1, int(args.thread))
        self.ddof = int(args.ddof)

        # Outlier filtering.
        self.remove_outlier = args.remove_outlier
        self.outlier_iqr = float(args.outlier_iqr)
        self.outlier_window = int(args.outlier_window)
        self.outlier_local_fold = float(args.outlier_local_fold)

        # Plotting.
        self.plot_stat = args.plot_stat
        self.plot_transform = args.plot_transform
        self.ci_low = float(args.ci_low)
        self.ci_high = float(args.ci_high)
        self.gene_fig = args.gene_fig
        self.gene_plot_dir = f"{self.output}_cumulative_CoV_geneplots"

        # Imported data.
        self.rpf_data = None
        self.sample_name = []
        self.sample_num = 0
        self.total_rpf_num = None
        self.file_format = None
        self.position_table = None

        # Derived results.
        self.sample_high_genes = OrderedDict()
        self.sample_summary = pd.DataFrame()
        self.cumulative_cov = pd.DataFrame()
        self.meta_cov = pd.DataFrame()
        self.outliers = pd.DataFrame()
        self.output_files = OrderedDict()

    # ------------------------------------------------------------------
    # Import and preparation
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import RPF density and construct codon- or nucleotide-level positions."""
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

        if self.resolution == "nucleotide":
            self.position_table = self._build_nucleotide_table()
        else:
            table = self.rpf_data.get_frame(frame=self.frame)
            shift_num = RPFs.set_codon_shift(self.site)
            table = RPFs.shift_site(table, self.sample_name, shift_num)
            self.position_table = self._prepare_codon_table(table)

        self._validate_table()
        self._select_sample_high_genes()

        if self.normal:
            for sample in self.sample_name:
                total = float(self.total_rpf_num.get(sample, 0.0))
                self.position_table[sample] = (
                    self.position_table[sample] * RPM_SCALE / total if total > 0 else 0.0
                )

        print(
            "Imported {samples} sample(s), format={fmt}, resolution={resolution}, rows={rows:,}.".format(
                samples=self.sample_num,
                fmt=self.file_format,
                resolution=self.resolution,
                rows=len(self.position_table),
            ),
            flush=True,
        )

    # Backward-compatible alias.
    retrieve_rpf = import_rpf

    def _build_nucleotide_table(self) -> pd.DataFrame:
        """Expand selected frame densities to ordered nucleotide positions."""
        frames = [0, 1, 2] if self.frame == "all" else [int(self.frame)]
        parts = []
        for frame in frames:
            table = self.rpf_data.get_frame(frame=str(frame)).copy()
            shift_num = RPFs.set_codon_shift(self.site)
            table = RPFs.shift_site(table, self.sample_name, shift_num)
            table = table.loc[table["region"] == "cds", BASE_COLUMNS + self.sample_name].copy()
            table["Frame"] = frame
            parts.append(table)

        if not parts:
            return pd.DataFrame()

        table = pd.concat(parts, ignore_index=True)
        table.sort_values(["name", "now_nt", "Frame"], inplace=True, kind="mergesort")

        # Build contiguous CDS-relative coordinates independently for each transcript.
        codon_order = (
            table[["name", "now_nt"]]
            .drop_duplicates()
            .sort_values(["name", "now_nt"], kind="mergesort")
        )
        codon_order["Codon"] = codon_order.groupby("name", sort=False).cumcount() + 1
        table = table.merge(codon_order, on=["name", "now_nt"], how="left", validate="many_to_one")
        table["Position"] = (table["Codon"] - 1) * 3 + table["Frame"] + 1
        table["Nucleotide"] = table["Position"]
        return table

    def _prepare_codon_table(self, table: pd.DataFrame) -> pd.DataFrame:
        """Prepare ordered codon-level coordinates."""
        table = table.loc[table["region"] == "cds", BASE_COLUMNS + self.sample_name].copy()
        table.sort_values(["name", "now_nt"], inplace=True, kind="mergesort")
        table["Codon"] = table.groupby("name", sort=False).cumcount() + 1
        table["Frame"] = 0
        table["Position"] = table["Codon"]
        table["Nucleotide"] = (table["Codon"] - 1) * 3 + 1
        return table

    def _validate_table(self) -> None:
        """Validate imported position-level table."""
        if self.position_table is None or self.position_table.empty:
            raise ValueError("No CDS positions remain after input filtering.")
        required = ["name", "Position", "Codon", "Frame"] + self.sample_name
        missing = [column for column in required if column not in self.position_table.columns]
        if missing:
            raise ValueError("Density table is missing required columns: " + ", ".join(missing))
        for sample in self.sample_name:
            self.position_table[sample] = pd.to_numeric(
                self.position_table[sample], errors="coerce"
            ).fillna(0.0).astype(float)

    def _select_sample_high_genes(self) -> None:
        """Select transcripts independently from raw sample-specific CDS counts."""
        counts = self.position_table.groupby("name", sort=False)[self.sample_name].sum()
        records = []
        for sample in self.sample_name:
            high = counts.index[counts[sample] >= float(self.rpf_num)]
            self.sample_high_genes[sample] = pd.Index(high)
            records.append(
                {
                    "Sample": sample,
                    "HighTranscriptCount": int(len(high)),
                    "MinimumRPFCount": float(self.rpf_num),
                }
            )
        self.sample_summary = pd.DataFrame(records)

    # ------------------------------------------------------------------
    # Cumulative CoV calculation
    # ------------------------------------------------------------------

    @staticmethod
    def _expanding_cov(values: np.ndarray, ddof: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate expanding mean, standard deviation, and CoV in O(n)."""
        values = values.astype(float, copy=False)
        valid = np.isfinite(values)
        clean = np.where(valid, values, 0.0)
        count = np.cumsum(valid.astype(np.int64))
        total = np.cumsum(clean)
        total_sq = np.cumsum(clean * clean)

        mean = np.full(values.shape, np.nan, dtype=float)
        np.divide(total, count, out=mean, where=count > 0)

        variance = np.full(values.shape, np.nan, dtype=float)
        denominator = count - ddof
        numerator = total_sq - np.divide(total * total, count, out=np.zeros_like(total), where=count > 0)
        np.divide(numerator, denominator, out=variance, where=denominator > 0)
        variance = np.maximum(variance, 0.0, where=np.isfinite(variance), out=variance)
        sd = np.sqrt(variance)

        cov = np.full(values.shape, np.nan, dtype=float)
        np.divide(sd, mean, out=cov, where=np.isfinite(sd) & np.isfinite(mean) & (mean > 0))
        return mean, sd, cov

    def _robust_cutoff(self, values: np.ndarray) -> float:
        """Calculate a robust upper cutoff on log1p-transformed positive values."""
        positive = values[np.isfinite(values) & (values > 0)]
        if positive.size < 8:
            return math.inf
        transformed = np.log1p(positive)
        q1, q3 = np.quantile(transformed, [0.25, 0.75])
        cutoff = q3 + self.outlier_iqr * (q3 - q1)
        return float(np.expm1(cutoff))

    def _detect_outliers(self, gene_table: pd.DataFrame, sample: str) -> tuple[np.ndarray, list[dict]]:
        """Detect extreme local RPF pileups for one transcript and sample."""
        values = gene_table[sample].to_numpy(dtype=float, copy=True)
        if not self.remove_outlier:
            return values, []

        cutoff = self._robust_cutoff(values)
        candidates = np.flatnonzero(values > cutoff)
        if candidates.size == 0:
            return values, []

        records = []
        n = len(values)
        for index in candidates:
            left = max(0, index - self.outlier_window)
            right = min(n, index + self.outlier_window + 1)
            neighbors = np.concatenate((values[left:index], values[index + 1 : right]))
            neighbors = neighbors[np.isfinite(neighbors)]
            local = float(np.mean(neighbors)) if neighbors.size else 0.0
            fold = (values[index] + 1.0) / (local + 1.0)
            if fold >= self.outlier_local_fold:
                row = gene_table.iloc[index]
                records.append(
                    {
                        "name": row["name"],
                        "Sample": sample,
                        "Position": int(row["Position"]),
                        "Codon": int(row["Codon"]),
                        "Frame": int(row["Frame"]),
                        "RawDensity": float(values[index]),
                        "RobustCutoff": cutoff,
                        "LocalBackground": local,
                        "LocalFold": float(fold),
                    }
                )
                values[index] = np.nan
        return values, records

    def _calculate_sample(self, sample: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        """Calculate cumulative CoV for one sample."""
        genes = self.sample_high_genes[sample]
        sample_table = self.position_table.loc[
            self.position_table["name"].isin(genes),
            ["name", "Position", "Nucleotide", "Codon", "Frame", "codon", sample],
        ].copy()
        sample_table.sort_values(["name", "Position"], inplace=True, kind="mergesort")

        results = []
        outliers = []
        retained = 0
        for gene, gene_table in sample_table.groupby("name", sort=False):
            if len(gene_table) < self.min_positions:
                continue
            values, gene_outliers = self._detect_outliers(gene_table, sample)
            mean, sd, cov = self._expanding_cov(values, self.ddof)
            result = gene_table.copy()
            result["Sample"] = sample
            result["RawDensity"] = gene_table[sample].to_numpy(dtype=float)
            result["FilteredDensity"] = values
            result["CumulativeMean"] = mean
            result["CumulativeSD"] = sd
            result["CumulativeCoV"] = cov
            result["ObservedPositionCount"] = np.cumsum(np.isfinite(values).astype(int))
            result.drop(columns=[sample], inplace=True)
            results.append(result)
            outliers.extend(gene_outliers)
            retained += 1

        result_table = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
        outlier_table = pd.DataFrame(outliers)
        summary = {
            "Sample": sample,
            "RetainedTranscriptCount": retained,
            "OutlierCount": int(len(outlier_table)),
        }
        return result_table, outlier_table, summary

    def calculate_cumulative_cov(self) -> None:
        """Calculate cumulative CoV independently for all samples."""
        workers = min(self.thread, max(1, self.sample_num))
        results = {}
        if workers == 1:
            for sample in self.sample_name:
                results[sample] = self._calculate_sample(sample)
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(self._calculate_sample, sample): sample
                    for sample in self.sample_name
                }
                for future in as_completed(futures):
                    sample = futures[future]
                    results[sample] = future.result()

        detail = []
        outliers = []
        summaries = []
        for sample in self.sample_name:
            table, outlier, summary = results[sample]
            if not table.empty:
                detail.append(table)
            if not outlier.empty:
                outliers.append(outlier)
            summaries.append(summary)

        self.cumulative_cov = pd.concat(detail, ignore_index=True) if detail else pd.DataFrame()
        self.outliers = pd.concat(outliers, ignore_index=True) if outliers else pd.DataFrame()
        detail_summary = pd.DataFrame(summaries)
        self.sample_summary = self.sample_summary.merge(detail_summary, on="Sample", how="left")

        if self.cumulative_cov.empty:
            raise ValueError("No transcript passed cumulative CoV filtering.")
        self._summarize_meta_cov()

    # Backward-compatible aliases.
    calc_cov = calculate_cumulative_cov

    def _summarize_meta_cov(self) -> None:
        """Summarize transcript-level cumulative CoV by sample and position."""
        table = self.cumulative_cov.loc[
            self.cumulative_cov["Position"] <= self.trim
        ].copy()
        if table.empty:
            raise ValueError("No cumulative CoV positions remain within --trim.")

        def summarize(group: pd.Series) -> pd.Series:
            values = group.to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                return pd.Series(
                    {
                        "TranscriptCount": 0,
                        "MeanCoV": np.nan,
                        "MedianCoV": np.nan,
                        "SDCoV": np.nan,
                        "LowerCoV": np.nan,
                        "UpperCoV": np.nan,
                    }
                )
            return pd.Series(
                {
                    "TranscriptCount": int(values.size),
                    "MeanCoV": float(np.mean(values)),
                    "MedianCoV": float(np.median(values)),
                    "SDCoV": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
                    "LowerCoV": float(np.quantile(values, self.ci_low)),
                    "UpperCoV": float(np.quantile(values, self.ci_high)),
                }
            )

        meta = (
            table.groupby(["Sample", "Position", "Nucleotide", "Codon", "Frame"], sort=False)[
                "CumulativeCoV"
            ]
            .apply(summarize)
            .unstack()
            .reset_index()
        )
        meta["Meta"] = "TIS"
        self.meta_cov = meta

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def output_results(self) -> None:
        """Write cumulative CoV result tables and summary JSON."""
        meta_file = f"{self.output}_meta{self.trim}_cumulative_CoV.txt"
        self.meta_cov.to_csv(meta_file, sep="\t", index=False)
        self.output_files["meta"] = meta_file

        if self.output_all:
            detail_file = f"{self.output}_cumulative_CoV.txt"
            self.cumulative_cov.to_csv(detail_file, sep="\t", index=False)
            self.output_files["detail"] = detail_file

        if self.remove_outlier and not self.outliers.empty:
            outlier_file = f"{self.output}_cumulative_CoV.outliers.txt"
            self.outliers.to_csv(outlier_file, sep="\t", index=False)
            self.output_files["outliers"] = outlier_file

        summary_file = f"{self.output}_cumulative_CoV.summary.json"
        payload = {
            "module": "rpf_Cumulative_CoV",
            "version": "0.2.8-dev.002",
            "input": self.rpf,
            "input_format": self.file_format,
            "samples": self.sample_name,
            "parameters": {
                "site": self.site,
                "frame": self.frame,
                "resolution": self.resolution,
                "min_rpf": self.rpf_num,
                "min_positions": self.min_positions,
                "trim": self.trim,
                "tis": self.tis,
                "tts": self.tts,
                "normal": self.normal,
                "ddof": self.ddof,
                "thread": self.thread,
                "remove_outlier": self.remove_outlier,
                "outlier_iqr": self.outlier_iqr,
                "outlier_window": self.outlier_window,
                "outlier_local_fold": self.outlier_local_fold,
                "gene_fig": self.gene_fig,
            },
            "sample_summary": self.sample_summary.fillna(0).to_dict(orient="records"),
            "outputs": self.output_files,
        }
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
        self.output_files["summary"] = summary_file

    # Backward-compatible aliases.
    merge_cov_table = lambda self: None
    output_rpf_table = output_results
    rpf_to_rpm = lambda self: None
    melt_rpf_table = lambda self: None

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _transform_plot_values(self, values: np.ndarray) -> np.ndarray:
        """Transform values for visualization only."""
        values = np.asarray(values, dtype=float)
        if self.plot_transform == "none":
            return values
        if self.plot_transform == "sqrt":
            return np.sqrt(np.clip(values, 0, None))
        if self.plot_transform == "log1p":
            return np.log1p(np.clip(values, 0, None))
        if self.plot_transform == "log2":
            return np.log2(np.clip(values, 0, None) + 1.0)
        if self.plot_transform == "log10":
            return np.log10(np.clip(values, 0, None) + 1.0)
        raise ValueError(f"Unknown plot transform: {self.plot_transform}")

    def draw_cumulative_cov(self) -> None:
        """Draw sample cumulative CoV curves with transcript quantile ribbons."""
        fig_width = max(8.5, 6.5 + 0.35 * self.sample_num)
        fig, ax = plt.subplots(figsize=(fig_width, 5.8), dpi=PLOT_DPI)

        center_column = "MedianCoV" if self.plot_stat == "median" else "MeanCoV"
        for sample in self.sample_name:
            data = self.meta_cov.loc[self.meta_cov["Sample"] == sample].sort_values("Position")
            if data.empty:
                continue
            x = data["Position"].to_numpy(dtype=float)
            center = self._transform_plot_values(data[center_column].to_numpy(dtype=float))
            lower = self._transform_plot_values(data["LowerCoV"].to_numpy(dtype=float))
            upper = self._transform_plot_values(data["UpperCoV"].to_numpy(dtype=float))
            line = ax.plot(x, center, linewidth=1.7, label=sample)[0]
            ax.fill_between(x, lower, upper, alpha=0.12, color=line.get_color(), linewidth=0)

        unit = "Nucleotide" if self.resolution == "nucleotide" else "Codon"
        ylabel = f"Cumulative CoV ({self.plot_stat})"
        if self.plot_transform != "none":
            ylabel += f" [{self.plot_transform}]"
        ax.set_xlabel(f"CDS-relative {unit.lower()} position")
        ax.set_ylabel(ylabel)
        ax.set_title("Cumulative coefficient of variation along CDS")
        ax.grid(axis="y", linewidth=0.5, alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if self.sample_num <= 12:
            ax.legend(frameon=False, ncol=min(4, self.sample_num), fontsize=8)
        else:
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        fig.tight_layout()

        pdf = f"{self.output}_cumulative_CoV_curve.pdf"
        png = f"{self.output}_cumulative_CoV_curve.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, bbox_inches="tight")
        plt.close(fig)
        self.output_files["curve_pdf"] = pdf
        self.output_files["curve_png"] = png

    @staticmethod
    def _safe_filename(value: str) -> str:
        """Convert a transcript identifier into a safe file name."""
        safe = str(value).strip()
        for token in (os.sep, os.altsep, ":", "*", "?", '"', "<", ">", "|"):
            if token:
                safe = safe.replace(token, "_")
        return safe or "unnamed_transcript"

    def draw_gene_cumulative_cov(self) -> None:
        """Draw one two-panel cumulative CoV figure for every retained transcript."""
        os.makedirs(self.gene_plot_dir, exist_ok=True)
        genes = self.cumulative_cov["name"].drop_duplicates().tolist()
        generated = 0

        for index, gene in enumerate(genes, start=1):
            gene_table = self.cumulative_cov.loc[
                self.cumulative_cov["name"] == gene
            ].copy()
            if gene_table.empty:
                continue

            fig, axes = plt.subplots(
                2,
                1,
                figsize=(9.0, 6.8),
                dpi=PLOT_DPI,
                sharex=True,
                gridspec_kw={"height_ratios": [1.0, 1.15], "hspace": 0.08},
            )
            density_ax, cov_ax = axes

            plotted_samples = 0
            for sample in self.sample_name:
                data = gene_table.loc[gene_table["Sample"] == sample].sort_values(
                    "Position", kind="mergesort"
                )
                if data.empty:
                    continue

                x = data["Position"].to_numpy(dtype=float)
                density = data["FilteredDensity"].to_numpy(dtype=float)
                cov = self._transform_plot_values(
                    data["CumulativeCoV"].to_numpy(dtype=float)
                )
                line = density_ax.plot(
                    x,
                    density,
                    linewidth=1.15,
                    alpha=0.9,
                    label=sample,
                )[0]
                cov_ax.plot(
                    x,
                    cov,
                    linewidth=1.45,
                    alpha=0.95,
                    color=line.get_color(),
                    label=sample,
                )
                plotted_samples += 1

            if plotted_samples == 0:
                plt.close(fig)
                continue

            density_ylabel = "RPF density"
            if self.normal:
                density_ylabel = "RPF density (RPM)"
            density_ax.set_ylabel(density_ylabel)
            density_ax.set_title(str(gene), loc="left", fontweight="bold")
            density_ax.grid(axis="y", linewidth=0.45, alpha=0.22)

            cov_ylabel = "Cumulative CoV"
            if self.plot_transform != "none":
                cov_ylabel += f" [{self.plot_transform}]"
            cov_ax.set_ylabel(cov_ylabel)
            unit = "nucleotide" if self.resolution == "nucleotide" else "codon"
            cov_ax.set_xlabel(f"CDS-relative {unit} position")
            cov_ax.grid(axis="y", linewidth=0.45, alpha=0.22)

            for axis in axes:
                axis.spines["top"].set_visible(False)
                axis.spines["right"].set_visible(False)

            handles, labels = cov_ax.get_legend_handles_labels()
            if handles:
                ncol = min(4, max(1, len(labels)))
                fig.legend(
                    handles,
                    labels,
                    loc="lower center",
                    bbox_to_anchor=(0.5, -0.01),
                    ncol=ncol,
                    frameon=False,
                    fontsize=7.5,
                )
                fig.subplots_adjust(bottom=0.14)

            safe_gene = self._safe_filename(gene)
            if self.gene_fig in {"png", "both"}:
                out_png = os.path.join(self.gene_plot_dir, f"{safe_gene}_cumulative_CoV.png")
                fig.savefig(out_png, bbox_inches="tight")
            if self.gene_fig in {"pdf", "both"}:
                out_pdf = os.path.join(self.gene_plot_dir, f"{safe_gene}_cumulative_CoV.pdf")
                fig.savefig(out_pdf, bbox_inches="tight")
            plt.close(fig)
            generated += 1

            if index % 100 == 0 or index == len(genes):
                print(
                    f"Generated cumulative CoV figures for {index:,}/{len(genes):,} transcripts.",
                    flush=True,
                )

        self.output_files["gene_plot_directory"] = self.gene_plot_dir
        self.output_files["gene_plot_count"] = generated

    def draw_transcript_count(self) -> None:
        """Draw the number of transcripts contributing at each position."""
        fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=PLOT_DPI)
        for sample in self.sample_name:
            data = self.meta_cov.loc[self.meta_cov["Sample"] == sample].sort_values("Position")
            if data.empty:
                continue
            ax.plot(
                data["Position"],
                data["TranscriptCount"],
                linewidth=1.5,
                label=sample,
            )
        ax.set_xlabel("CDS-relative position")
        ax.set_ylabel("Contributing transcripts")
        ax.set_title("Transcript support across cumulative CoV positions")
        ax.grid(axis="y", linewidth=0.5, alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if self.sample_num <= 12:
            ax.legend(frameon=False, ncol=min(4, self.sample_num), fontsize=8)
        else:
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        fig.tight_layout()

        pdf = f"{self.output}_cumulative_CoV_support.pdf"
        png = f"{self.output}_cumulative_CoV_support.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, bbox_inches="tight")
        plt.close(fig)
        self.output_files["support_pdf"] = pdf
        self.output_files["support_png"] = png
