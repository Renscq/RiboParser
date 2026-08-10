#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-09
# Version: dev002
# Function: Draw compact shared, specific, significance, length, and enrichment summary figures.
# Input: Precomputed pairwise SeRP overlap-summary tables.
# Output: Compact publication-oriented SeRP overlap summary figures.

"""Plotting layer for pairwise SeRP overlap-result summarization."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


SAMPLE_COLORS = ["#0072B2", "#D55E00"]
SHARED_COLOR = "#009E73"
SPECIFIC_COLOR = "#CC79A7"
SIGNIFICANCE_COLORS = {
    "significant": "#009E73",
    "not_significant": "#E69F00",
    "untested": "#BDBDBD",
}


class SeRPOverlapSummaryPlot:
    """Draw compact summary figures from pairwise SeRP overlap analysis."""

    def __init__(
        self,
        peak_table: pd.DataFrame,
        category_summary: pd.DataFrame,
        gene_summary: pd.DataFrame,
        cluster_summary: pd.DataFrame,
        relationships: pd.DataFrame,
        sample_order: list[str],
        output_prefix: str,
        top_specific_genes: int,
        output_format: str,
        dpi: int,
        font_size: float,
    ) -> None:
        """Initialize plotting data and rendering options."""
        self.peak_table = peak_table
        self.category_summary = category_summary
        self.gene_summary = gene_summary
        self.cluster_summary = cluster_summary
        self.relationships = relationships
        self.sample_order = list(sample_order)
        self.output_prefix = str(output_prefix)
        self.top_specific_genes = int(top_specific_genes)
        self.output_format = str(output_format)
        self.dpi = int(dpi)
        self.font_size = float(font_size)

    @staticmethod
    def _sanitize_label(label: str) -> str:
        """Convert a condition label to the filename-safe summary key."""
        import re

        cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(label).strip())
        cleaned = cleaned.strip("._-")
        return cleaned if cleaned else "sample"

    def _save_figure(self, fig, suffix: str) -> None:
        """Save one figure in the requested output format."""
        output_stem = Path(self.output_prefix + suffix)
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
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0), length=3.0)
        ax.xaxis.label.set_size(self.font_size)
        ax.yaxis.label.set_size(self.font_size)

    def _create_figure(self, width: float, height: float):
        """Create a compact figure canvas with constrained layout."""
        return plt.subplots(
            figsize=(width, height),
            dpi=self.dpi,
            constrained_layout=True,
        )

    def _draw_violin_box(
        self,
        ax,
        groups: list[np.ndarray],
        labels: list[str],
        colors: list[str],
        ylabel: str,
        baseline: float | None = None,
    ) -> None:
        """Draw a compact violin-plus-boxplot composite for grouped distributions."""
        positions = np.arange(1, len(groups) + 1, dtype=float)
        violin = ax.violinplot(
            groups,
            positions=positions,
            widths=0.82,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body, color in zip(violin["bodies"], colors):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.28)
            body.set_linewidth(0.6)

        box = ax.boxplot(
            groups,
            positions=positions,
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.0},
            whiskerprops={"linewidth": 0.8},
            capprops={"linewidth": 0.8},
            boxprops={"linewidth": 0.8},
        )
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.82)
            patch.set_edgecolor(color)

        if baseline is not None:
            ax.axhline(baseline, color="grey", linewidth=0.6, linestyle="--")
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        self._setup_axis(ax)

    @staticmethod
    def _label_bars(ax, bars) -> None:
        """Add integer count labels above vertical bars."""
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                "{0:,}".format(int(round(height))),
                ha="center",
                va="bottom",
                fontsize=8,
            )

    def _plot_peak_category_counts(self) -> None:
        """Plot shared and condition-specific peak counts for each condition."""
        rows = self.category_summary.loc[
            self.category_summary["category"].isin(["shared", "specific"])
        ]
        if rows.empty:
            return

        x = np.arange(len(self.sample_order), dtype=float)
        width = 0.34
        shared_values = []
        specific_values = []
        for sample in self.sample_order:
            one = rows.loc[rows["sample"].eq(sample)].set_index("category")
            shared_values.append(int(one.loc["shared", "peak_count"]) if "shared" in one.index else 0)
            specific_values.append(int(one.loc["specific", "peak_count"]) if "specific" in one.index else 0)

        fig, ax = self._create_figure(4.8, 3.3)
        bars_shared = ax.bar(
            x - width / 2.0,
            shared_values,
            width=width,
            label="Shared",
            color=SHARED_COLOR,
        )
        bars_specific = ax.bar(
            x + width / 2.0,
            specific_values,
            width=width,
            label="Specific",
            color=SPECIFIC_COLOR,
        )
        self._label_bars(ax, bars_shared)
        self._label_bars(ax, bars_specific)
        ax.set_xticks(x)
        ax.set_xticklabels(self.sample_order)
        ax.set_ylabel("Peak count")
        ax.legend(frameon=False, fontsize=max(self.font_size - 1.0, 1.0))
        self._setup_axis(ax)
        self._save_figure(fig, ".peak_category_counts")

    def _plot_significance_counts(self) -> None:
        """Plot significant, non-significant, and untested peaks by category."""
        groups: list[tuple[str, str]] = []
        for sample in self.sample_order:
            groups.extend([(sample, "shared"), (sample, "specific")])

        labels = ["{0}\n{1}".format(sample, category.title()) for sample, category in groups]
        statuses = ["significant", "not_significant", "untested"]
        values = {status: [] for status in statuses}
        for sample, category in groups:
            subset = self.peak_table.loc[
                self.peak_table["comparison_sample"].astype(str).eq(sample)
                & self.peak_table["peak_category"].eq(category)
            ]
            counts = subset["significance_status"].value_counts()
            for status in statuses:
                values[status].append(int(counts.get(status, 0)))

        x = np.arange(len(groups), dtype=float)
        bottom = np.zeros(len(groups), dtype=float)
        fig, ax = self._create_figure(6.0, 3.5)
        for status in statuses:
            bars = ax.bar(
                x,
                values[status],
                bottom=bottom,
                label=status.replace("_", " ").title(),
                color=SIGNIFICANCE_COLORS[status],
                width=0.7,
            )
            bottom += np.asarray(values[status], dtype=float)
            for bar, count in zip(bars, values[status]):
                if count <= 0:
                    continue
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bar.get_y() + bar.get_height() / 2.0,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=7.5,
                )
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Peak count")
        ax.legend(frameon=False, fontsize=max(self.font_size - 1.0, 1.0))
        self._setup_axis(ax)
        self._save_figure(fig, ".significance_counts")

    def _plot_peak_length(self) -> None:
        """Compare core peak lengths between shared and specific categories."""
        groups = []
        labels = []
        for sample in self.sample_order:
            for category in ["shared", "specific"]:
                values = pd.to_numeric(
                    self.peak_table.loc[
                        self.peak_table["comparison_sample"].astype(str).eq(sample)
                        & self.peak_table["peak_category"].eq(category),
                        "summary_peak_length",
                    ],
                    errors="coerce",
                ).dropna()
                if not values.empty:
                    groups.append(values.to_numpy(dtype=float))
                    labels.append("{0}\n{1}".format(sample, category.title()))
        if not groups:
            return

        colors = [SHARED_COLOR if "Shared" in label else SPECIFIC_COLOR for label in labels]
        fig, ax = self._create_figure(5.8, 3.5)
        self._draw_violin_box(
            ax=ax,
            groups=groups,
            labels=labels,
            colors=colors,
            ylabel="Peak core length (codons)",
        )
        self._save_figure(fig, ".peak_length_by_category")

    def _plot_fold_distribution(self) -> None:
        """Compare maximum enrichment between shared and specific peaks."""
        if "max_fold" not in self.peak_table.columns:
            return
        groups = []
        labels = []
        for sample in self.sample_order:
            for category in ["shared", "specific"]:
                values = pd.to_numeric(
                    self.peak_table.loc[
                        self.peak_table["comparison_sample"].astype(str).eq(sample)
                        & self.peak_table["peak_category"].eq(category),
                        "max_fold",
                    ],
                    errors="coerce",
                ).dropna()
                values = values.loc[values > 0]
                if not values.empty:
                    groups.append(np.log2(values.to_numpy(dtype=float)))
                    labels.append("{0}\n{1}".format(sample, category.title()))
        if not groups:
            return

        colors = [SHARED_COLOR if "Shared" in label else SPECIFIC_COLOR for label in labels]
        fig, ax = self._create_figure(5.8, 3.5)
        self._draw_violin_box(
            ax=ax,
            groups=groups,
            labels=labels,
            colors=colors,
            ylabel="log2(max fold enrichment)",
            baseline=0.0,
        )
        self._save_figure(fig, ".max_fold_by_category")

    def _plot_cluster_types(self) -> None:
        """Plot 1:1, 1:N, N:1, and N:N shared-cluster counts."""
        if self.cluster_summary.empty or "relationship_class" not in self.cluster_summary.columns:
            return
        order = ["1:1", "1:N", "N:1", "N:N", "other"]
        counts = self.cluster_summary["relationship_class"].value_counts()
        labels = [item for item in order if counts.get(item, 0) > 0]
        if not labels:
            return
        values = [int(counts[item]) for item in labels]

        fig, ax = self._create_figure(4.8, 3.3)
        bars = ax.bar(labels, values, color=CLUSTER_COLOR, width=0.68)
        self._label_bars(ax, bars)
        ax.set_xlabel("Shared-cluster relationship")
        ax.set_ylabel("Cluster count")
        self._setup_axis(ax)
        self._save_figure(fig, ".shared_cluster_types")

    def _plot_specific_gene_burden(self) -> None:
        """Plot the genes carrying the largest numbers of condition-specific peaks."""
        if self.gene_summary.empty or self.top_specific_genes <= 0:
            return
        key_a = self._sanitize_label(self.sample_order[0])
        key_b = self._sanitize_label(self.sample_order[1])
        column_a = key_a + "_specific_peak_count"
        column_b = key_b + "_specific_peak_count"
        if column_a not in self.gene_summary.columns or column_b not in self.gene_summary.columns:
            return

        table = self.gene_summary[["gene_name", column_a, column_b]].copy()
        table["total_specific"] = table[column_a] + table[column_b]
        table = table.loc[table["total_specific"] > 0]
        table = table.sort_values(
            ["total_specific", column_a, column_b],
            ascending=[False, False, False],
            kind="stable",
        ).head(self.top_specific_genes)
        if table.empty:
            return
        table = table.sort_values("total_specific", ascending=True, kind="stable")

        y = np.arange(len(table), dtype=float)
        fig_height = max(3.4, min(7.2, 2.1 + len(table) * 0.24))
        fig, ax = self._create_figure(6.0, fig_height)
        ax.barh(
            y,
            -table[column_a].to_numpy(dtype=float),
            color=SAMPLE_COLORS[0],
            label=self.sample_order[0],
        )
        ax.barh(
            y,
            table[column_b].to_numpy(dtype=float),
            color=SAMPLE_COLORS[1],
            label=self.sample_order[1],
        )
        ax.axvline(0.0, color="grey", linewidth=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(table["gene_name"].astype(str))
        max_abs = max(
            float(table[column_a].max()),
            float(table[column_b].max()),
            1.0,
        )
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda value, position: "{0:g}".format(abs(value)))
        )
        ax.set_xlim(-max_abs * 1.2, max_abs * 1.2)
        ax.set_xlabel("Condition-specific peak count")
        ax.legend(frameon=False, fontsize=max(self.font_size - 1.0, 1.0))
        self._setup_axis(ax)
        self._save_figure(fig, ".top_specific_genes")

    def _plot_overlap_quality(self) -> None:
        """Plot reciprocal-overlap quality for qualifying shared relationships."""
        if self.relationships.empty or "reciprocal_overlap" not in self.relationships.columns:
            return
        table = self.relationships.copy()
        if "qualifying_overlap" in table.columns:
            mask = table["qualifying_overlap"].astype(str).str.lower().isin(
                ["true", "1", "yes"]
            )
            table = table.loc[mask]
        values = pd.to_numeric(table["reciprocal_overlap"], errors="coerce").dropna()
        if values.empty:
            return
        bins = max(10, min(40, int(np.sqrt(len(values))) + 1))

        fig, ax = self._create_figure(4.8, 3.3)
        ax.hist(
            values.to_numpy(dtype=float),
            bins=bins,
            range=(0.0, 1.0),
            color=SHARED_COLOR,
            edgecolor="white",
            linewidth=0.5,
        )
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Reciprocal overlap")
        ax.set_ylabel("Shared peak-pair count")
        self._setup_axis(ax)
        self._save_figure(fig, ".reciprocal_overlap_distribution")

    def draw_all(self) -> None:
        """Generate all default compact pairwise SeRP overlap-summary figures."""
        self._plot_peak_category_counts()
        self._plot_significance_counts()
        self._plot_peak_length()
        self._plot_fold_distribution()
        self._plot_specific_gene_burden()
        self._plot_overlap_quality()