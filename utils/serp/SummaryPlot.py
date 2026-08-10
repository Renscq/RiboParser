#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-08
# Version: dev001
# Function: Draw global summary figures for SeRP peak-calling results.
# Input: Peak-level and transcript-level SeRP summary tables.
# Output: Peak count, length, enrichment, positional, and heatmap figures.

"""Plotting layer for SeRP peak-result summarization."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


PEAK_COLOR = "#0072B2"
SECONDARY_COLOR = "#D55E00"
GENE_COLOR = "#D9D9D9"


class SeRPSummaryPlot:
    """Draw summary figures from precomputed SeRP peak statistics."""

    def __init__(
        self,
        peak_table: pd.DataFrame,
        transcript_summary: pd.DataFrame,
        normalized_matrix: pd.DataFrame,
        inter_peak_distance: pd.DataFrame,
        output_prefix: str,
        normalized_bins: int,
        max_heatmap_cells: int,
        output_format: str,
        dpi: int,
        font_size: float,
    ) -> None:
        """Initialize plotting data and rendering options."""
        self.peak_table = peak_table
        self.transcript_summary = transcript_summary
        self.normalized_matrix = normalized_matrix
        self.inter_peak_distance = inter_peak_distance
        self.output_prefix = output_prefix
        self.normalized_bins = int(normalized_bins)
        self.max_heatmap_cells = int(max_heatmap_cells)
        self.output_format = output_format
        self.dpi = int(dpi)
        self.font_size = float(font_size)

    def _save_figure(self, fig, suffix: str) -> None:
        """Save one summary figure in the requested format."""
        output_stem = Path(str(self.output_prefix) + suffix)
        if self.output_format in {"png", "both"}:
            fig.savefig(
                Path(str(output_stem) + ".png"),
                dpi=self.dpi,
                bbox_inches="tight",
            )
        if self.output_format in {"pdf", "both"}:
            fig.savefig(Path(str(output_stem) + ".pdf"), bbox_inches="tight")
        plt.close(fig)

    def _setup_axis(self, ax) -> None:
        """Apply a restrained publication-oriented axis style."""
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        ax.xaxis.label.set_size(self.font_size)
        ax.yaxis.label.set_size(self.font_size)

    def _plot_summary_counts(self) -> None:
        """Plot primary peak, transcript, gene, and multi-peak counts."""
        peak = self.peak_table
        tx = self.transcript_summary
        valid_gene = peak["gene_name"].astype(str).str.strip()
        valid_gene = valid_gene.loc[~valid_gene.isin(["", "-", "nan", "None"])]
        labels = ["Peaks", "Transcripts", "Genes", "Multi-peak\ntranscripts"]
        values = [
            len(peak),
            peak["transcripts"].nunique(),
            valid_gene.nunique(),
            int((tx["peak_count"] > 1).sum()),
        ]
        fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=self.dpi)
        bars = ax.bar(labels, values, color=PEAK_COLOR, width=0.72)
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height(),
                "{0:,}".format(int(value)),
                ha="center",
                va="bottom",
                fontsize=max(self.font_size - 1.0, 1.0),
            )
        ax.set_ylabel("Count")
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".summary_counts")

    def _plot_fold_distribution(self) -> None:
        """Plot max-fold and mean-fold distributions when available."""
        available = []
        for column, label, color in [
            ("max_fold", "Max fold", PEAK_COLOR),
            ("mean_fold", "Mean fold", SECONDARY_COLOR),
        ]:
            if column in self.peak_table.columns:
                values = self.peak_table[column].dropna().to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                if len(values):
                    available.append((values, label, color))
        if not available:
            return

        combined = np.concatenate([item[0] for item in available])
        bin_number = max(10, min(60, int(np.sqrt(len(combined))) + 1))
        bins = np.histogram_bin_edges(combined, bins=bin_number)
        fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=self.dpi)
        for values, label, color in available:
            ax.hist(
                values,
                bins=bins,
                histtype="step",
                linewidth=1.4,
                color=color,
                label=label,
            )
        ax.set_xlabel("Peak enrichment fold")
        ax.set_ylabel("Peak count")
        ax.legend(frameon=False, fontsize=max(self.font_size - 1.0, 1.0))
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".peak_fold_distribution")

    def _plot_peak_length(self) -> None:
        """Plot core-length and interval-span distributions."""
        core = self.peak_table["peak_length"].dropna().to_numpy(dtype=float)
        span = self.peak_table["peak_span_length"].dropna().to_numpy(dtype=float)
        combined = np.concatenate([core, span])
        bin_number = max(10, min(60, int(np.sqrt(len(combined))) + 1))
        bins = np.histogram_bin_edges(combined, bins=bin_number)

        fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=self.dpi)
        ax.hist(
            core,
            bins=bins,
            histtype="step",
            linewidth=1.4,
            color=PEAK_COLOR,
            label="Core length",
        )
        ax.hist(
            span,
            bins=bins,
            histtype="step",
            linewidth=1.4,
            color=SECONDARY_COLOR,
            label="Interval span",
        )
        ax.set_xlabel("Peak length (codons)")
        ax.set_ylabel("Peak count")
        ax.legend(frameon=False, fontsize=max(self.font_size - 1.0, 1.0))
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".peak_length_distribution")

    def _plot_peak_count(self) -> None:
        """Plot the number of peaks per transcript."""
        distribution = self.transcript_summary["peak_count"].value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(5.0, 3.8), dpi=self.dpi)
        ax.bar(
            distribution.index.astype(str),
            distribution.to_numpy(),
            color=PEAK_COLOR,
            width=0.8,
        )
        ax.set_xlabel("Peaks per transcript")
        ax.set_ylabel("Transcript count")
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".peaks_per_transcript")

    def _plot_peak_region(self) -> None:
        """Plot broad transcript-region classes of called peaks."""
        distribution = self.peak_table["peak_region_class"].value_counts()
        preferred = ["utr5", "cds", "utr3", "mixed", "unknown"]
        ordered = [value for value in preferred if value in distribution.index]
        ordered.extend([value for value in distribution.index if value not in ordered])
        values = distribution.reindex(ordered)

        fig, ax = plt.subplots(figsize=(5.0, 3.8), dpi=self.dpi)
        ax.bar(values.index, values.to_numpy(), color=PEAK_COLOR, width=0.75)
        ax.set_xlabel("Peak region")
        ax.set_ylabel("Peak count")
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".peak_region_distribution")

    def _plot_peak_center(self) -> None:
        """Plot peak-center positions after transcript-length normalization."""
        values = self.peak_table["peak_center_percent"].to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=self.dpi)
        ax.hist(
            values,
            bins=np.linspace(0.0, 100.0, 21),
            color=PEAK_COLOR,
            edgecolor="white",
            linewidth=0.5,
        )
        ax.set_xlim(0.0, 100.0)
        ax.set_xlabel("Normalized peak center (%)")
        ax.set_ylabel("Peak count")
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".peak_center_distribution")

    def _plot_peak_burden(self) -> None:
        """Plot transcript length against the number of called peaks."""
        table = self.transcript_summary
        fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=self.dpi)
        ax.scatter(
            table["transcript_length_codons"],
            table["peak_count"],
            s=9,
            alpha=0.55,
            color=PEAK_COLOR,
            linewidths=0,
        )
        ax.set_xlabel("Transcript length (codons)")
        ax.set_ylabel("Peak count")
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".peak_burden")

    def _plot_inter_peak_distance(self) -> None:
        """Plot distances between adjacent peaks on multi-peak transcripts."""
        if self.inter_peak_distance.empty:
            return
        values = self.inter_peak_distance["inter_peak_distance_codons"].to_numpy(dtype=float)
        bin_number = max(10, min(60, int(np.sqrt(len(values))) + 1))
        fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=self.dpi)
        ax.hist(values, bins=bin_number, color=PEAK_COLOR, edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Inter-peak distance (codons)")
        ax.set_ylabel("Peak-pair count")
        self._setup_axis(ax)
        fig.tight_layout()
        self._save_figure(fig, ".inter_peak_distance")

    def _plot_normalized_heatmap(self) -> None:
        """Plot transcript-normalized peak coverage using fixed bins."""
        matrix = self.normalized_matrix.iloc[:, 4:].to_numpy(dtype=float)
        cmap = LinearSegmentedColormap.from_list(
            "serp_peak_coverage",
            ["#F2F2F2", PEAK_COLOR],
        )
        height = min(10.0, max(4.0, 3.0 + len(matrix) / 1200.0))
        fig, ax = plt.subplots(figsize=(6.2, height), dpi=self.dpi)
        image = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            rasterized=True,
        )
        tick_positions = np.linspace(0, self.normalized_bins - 1, 5)
        tick_labels = ["0", "25", "50", "75", "100"]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.set_yticks([])
        ax.set_xlabel("Normalized transcript position (%)")
        ax.set_ylabel("Transcripts sorted by peak position")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        colorbar.set_label("Peak coverage within bin", fontsize=self.font_size)
        colorbar.ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        fig.tight_layout()
        self._save_figure(fig, ".normalized_peak_heatmap")

    def _plot_full_length_matrix(self, summary: pd.DataFrame) -> bool:
        """Plot the full-length peak map as an exact codon-resolution heatmap.

        Returns
        -------
        bool
            ``True`` when the matrix renderer was used, otherwise ``False`` when
            the requested matrix would exceed the configured cell limit.
        """
        max_length = int(summary["transcript_length_codons"].max())
        cell_number = int(len(summary) * max_length)
        if cell_number > self.max_heatmap_cells:
            return False

        matrix = np.full((len(summary), max_length), 255, dtype=np.uint8)
        peak_group = {
            transcript: table
            for transcript, table in self.peak_table.groupby("transcripts", sort=False)
        }
        for row_index, row in summary.iterrows():
            length = int(row["transcript_length_codons"])
            matrix[row_index, :length] = 0
            table = peak_group[str(row["transcripts"])]
            for peak_row in table.itertuples():
                start = int(peak_row.peak_local_start)
                end = int(peak_row.peak_local_end)
                matrix[row_index, start : end + 1] = 1

        masked = np.ma.masked_where(matrix == 255, matrix.astype(float))
        cmap = LinearSegmentedColormap.from_list(
            "serp_full_peak",
            [GENE_COLOR, PEAK_COLOR],
        )
        cmap.set_bad("white")
        height = min(10.0, max(4.0, 3.0 + len(summary) / 1200.0))
        fig, ax = plt.subplots(figsize=(7.0, height), dpi=self.dpi)
        ax.imshow(
            masked,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            rasterized=True,
        )
        ax.set_yticks([])
        ax.set_xlabel("Transcript position from 5' end (codons)")
        ax.set_ylabel("Transcripts sorted by length")
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        fig.tight_layout()
        self._save_figure(fig, ".full_length_peak_heatmap")
        return True

    def _plot_full_length_intervals(self, summary: pd.DataFrame) -> None:
        """Render all full-length transcripts without allocating a dense matrix."""
        gene_segments = []
        peak_segments = []
        peak_group = {
            transcript: table
            for transcript, table in self.peak_table.groupby("transcripts", sort=False)
        }

        for row_index, row in summary.iterrows():
            length = int(row["transcript_length_codons"])
            y = float(row_index)
            gene_segments.append([(0.0, y), (float(length), y)])
            table = peak_group[str(row["transcripts"])]
            for peak_row in table.itertuples():
                peak_segments.append(
                    [
                        (float(peak_row.peak_local_start), y),
                        (float(peak_row.peak_local_end + 1), y),
                    ]
                )

        height = min(10.0, max(4.0, 3.0 + len(summary) / 1200.0))
        fig, ax = plt.subplots(figsize=(7.0, height), dpi=self.dpi)
        ax.add_collection(
            LineCollection(
                gene_segments,
                colors=GENE_COLOR,
                linewidths=0.45,
                rasterized=True,
            )
        )
        if peak_segments:
            ax.add_collection(
                LineCollection(
                    peak_segments,
                    colors=PEAK_COLOR,
                    linewidths=0.85,
                    rasterized=True,
                )
            )
        ax.set_xlim(0.0, float(summary["transcript_length_codons"].max()))
        ax.set_ylim(-1.0, float(len(summary)))
        ax.invert_yaxis()
        ax.set_yticks([])
        ax.set_xlabel("Transcript position from 5' end (codons)")
        ax.set_ylabel("Transcripts sorted by length")
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        self._save_figure(fig, ".full_length_peak_heatmap")

    def _plot_full_length_heatmap(self) -> None:
        """Plot all peak-containing transcripts at their actual codon lengths."""
        summary = self.transcript_summary.sort_values(
            ["transcript_length_codons", "weighted_peak_center_percent"],
            ascending=[False, True],
            kind="stable",
        ).reset_index(drop=True)

        if self._plot_full_length_matrix(summary):
            return

        print(
            "Full-length matrix exceeds {0:,} cells; using an exact interval "
            "renderer without dropping transcripts.".format(self.max_heatmap_cells),
            flush=True,
        )
        self._plot_full_length_intervals(summary)


    def draw_all(self) -> None:
        """Generate all default SeRP summary figures."""
        self._plot_summary_counts()
        self._plot_peak_length()
        self._plot_fold_distribution()
        self._plot_peak_count()
        self._plot_peak_region()
        self._plot_peak_center()
        self._plot_peak_burden()
        self._plot_inter_peak_distance()
        self._plot_normalized_heatmap()
        self._plot_full_length_heatmap()