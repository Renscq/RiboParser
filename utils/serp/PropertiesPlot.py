#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev001
# Function: Draw group-comparison sequence-property figures for SeRP peaks.
# Input: Precomputed full-protein, local-context, codon, amino-acid, and position-profile tables.
# Output: Compact publication-oriented SeRP sequence-property figures.

"""Plotting layer for SeRP shared and condition-specific sequence properties."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from utils.serp.Properties import AA_ORDER


GROUP_COLORS = ["#009E73", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9"]


class SeRPPropertiesPlot:
    """Draw compact sequence-property comparisons for SeRP peak groups."""

    def __init__(
        self,
        analyzer,
        output_prefix: str,
        output_format: str,
        dpi: int,
        font_size: float,
    ) -> None:
        """Initialize plotting tables and rendering options."""
        self.analyzer = analyzer
        self.output_prefix = str(output_prefix)
        self.output_format = str(output_format)
        self.dpi = int(dpi)
        self.font_size = float(font_size)
        self.groups = list(analyzer.plot_groups)
        self.group_colors = {
            group: GROUP_COLORS[index % len(GROUP_COLORS)]
            for index, group in enumerate(self.groups)
        }

    @staticmethod
    def _group_label(group: str) -> str:
        """Create a readable figure label from one analysis-group name."""
        if group == "shared":
            return "Shared"
        if group.endswith("_significant_specific"):
            sample = group[: -len("_significant_specific")]
            return "{0}\nSignificant specific".format(sample)
        if group.endswith("_specific"):
            sample = group[: -len("_specific")]
            return "{0}\nSpecific".format(sample)
        return group

    def _save_figure(self, fig, suffix: str) -> None:
        """Save one figure in the requested format and close it."""
        stem = Path(self.output_prefix + suffix)
        if self.output_format in {"png", "both"}:
            fig.savefig(
                Path(str(stem) + ".png"),
                dpi=self.dpi,
                bbox_inches="tight",
            )
        if self.output_format in {"pdf", "both"}:
            fig.savefig(Path(str(stem) + ".pdf"), bbox_inches="tight")
        plt.close(fig)

    def _setup_axis(self, ax) -> None:
        """Apply a restrained publication-oriented axis style."""
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0), length=3.0)
        ax.xaxis.label.set_size(self.font_size)
        ax.yaxis.label.set_size(self.font_size)
        ax.title.set_size(self.font_size)

    def _draw_violin_box(
        self,
        ax,
        table: pd.DataFrame,
        group_column: str,
        feature: str,
        ylabel: str,
    ) -> None:
        """Draw violin plus boxplots with robust handling of small groups."""
        values_by_group = []
        positions = []
        labels = []
        colors = []
        for position, group in enumerate(self.groups, start=1):
            values = pd.to_numeric(
                table.loc[table[group_column].eq(group), feature],
                errors="coerce",
            ).dropna()
            if values.empty:
                continue
            values_by_group.append(values.to_numpy(dtype=float))
            positions.append(float(position))
            labels.append(self._group_label(group))
            colors.append(self.group_colors[group])

        if not values_by_group:
            ax.set_visible(False)
            return

        violin_values = []
        violin_positions = []
        violin_colors = []
        for values, position, color in zip(values_by_group, positions, colors):
            if len(values) >= 2 and np.nanstd(values) > 0:
                violin_values.append(values)
                violin_positions.append(position)
                violin_colors.append(color)
        if violin_values:
            violin = ax.violinplot(
                violin_values,
                positions=violin_positions,
                widths=0.82,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body, color in zip(violin["bodies"], violin_colors):
                body.set_facecolor(color)
                body.set_edgecolor(color)
                body.set_alpha(0.28)
                body.set_linewidth(0.6)

        box = ax.boxplot(
            values_by_group,
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
            patch.set_edgecolor(color)
            patch.set_alpha(0.82)

        for values, position, color in zip(values_by_group, positions, colors):
            if len(values) == 1:
                ax.scatter(
                    [position],
                    values,
                    s=18,
                    color=color,
                    edgecolor="black",
                    linewidth=0.3,
                    zorder=4,
                )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        self._setup_axis(ax)

    def _plot_property_grid(
        self,
        table: pd.DataFrame,
        group_column: str,
        features: list[tuple[str, str]],
        suffix: str,
    ) -> None:
        """Draw a compact two-by-three grid of violin-plus-box property panels."""
        available = [item for item in features if item[0] in table.columns]
        if not available:
            return
        fig, axes = plt.subplots(
            2,
            3,
            figsize=(9.0, 5.8),
            dpi=self.dpi,
            constrained_layout=True,
        )
        axes_flat = axes.ravel()
        for ax, (feature, ylabel) in zip(axes_flat, available):
            self._draw_violin_box(
                ax=ax,
                table=table,
                group_column=group_column,
                feature=feature,
                ylabel=ylabel,
            )
        for ax in axes_flat[len(available):]:
            ax.set_visible(False)
        self._save_figure(fig, suffix)

    def _plot_full_properties(self) -> None:
        """Plot core full-protein physicochemical and coding-sequence properties."""
        table = self.analyzer.selected_full_properties()
        if table.empty:
            return
        self._plot_property_grid(
            table=table,
            group_column="Group",
            features=[
                ("Length", "Protein length (aa)"),
                ("Gravy", "GRAVY"),
                ("Charge_Density", "Charge density at configured pH"),
                ("Isoelectric_Point", "Isoelectric point"),
                ("Aromaticity", "Aromaticity"),
                ("GC3", "GC3 fraction"),
            ],
            suffix=".protein_properties",
        )

    def _plot_local_properties(self) -> None:
        """Plot max-site-centered local sequence properties."""
        table = self.analyzer.selected_local_properties()
        if table.empty:
            return
        self._plot_property_grid(
            table=table,
            group_column="Plot_Group",
            features=[
                ("Gravy", "Local GRAVY"),
                ("Signed_Charge_Fraction", "Signed charge-residue fraction"),
                ("Positive_Fraction", "Positive-residue fraction"),
                ("Aromatic_Fraction", "Aromatic-residue fraction"),
                ("Proline_Fraction", "Proline fraction"),
                ("GC3", "Local GC3 fraction"),
            ],
            suffix=".local_properties",
        )

    def _plot_aa_composition_heatmap(self) -> None:
        """Plot full-protein amino-acid composition enrichment versus selected-group pool."""
        table = self.analyzer.group_aa_composition.copy()
        table = table.loc[table["Group"].isin(self.groups)]
        if table.empty:
            return
        pivot = table.pivot(index="Group", columns="AA", values="Count").fillna(0.0)
        pivot = pivot.reindex(index=self.groups, columns=AA_ORDER, fill_value=0.0)
        group_total = pivot.sum(axis=1)
        group_fraction = (pivot + 0.5).div(group_total + 0.5 * len(AA_ORDER), axis=0)
        pooled = pivot.sum(axis=0)
        pooled_fraction = (pooled + 0.5) / (pooled.sum() + 0.5 * len(AA_ORDER))
        matrix = np.log2(group_fraction.div(pooled_fraction, axis=1)).to_numpy(dtype=float)
        vmax = max(float(np.nanmax(np.abs(matrix))), 0.5)

        fig, ax = plt.subplots(
            figsize=(8.5, max(2.3, 1.2 + len(self.groups) * 0.62)),
            dpi=self.dpi,
            constrained_layout=True,
        )
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap="coolwarm",
            norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
            interpolation="nearest",
        )
        ax.set_xticks(np.arange(len(AA_ORDER)))
        ax.set_xticklabels(AA_ORDER)
        ax.set_yticks(np.arange(len(self.groups)))
        ax.set_yticklabels([self._group_label(group).replace("\n", " ") for group in self.groups])
        ax.set_xlabel("Amino acid")
        ax.set_ylabel("")
        ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0), length=2.0)
        cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        cbar.set_label("log2(frequency / pooled frequency)", fontsize=self.font_size)
        cbar.ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        self._save_figure(fig, ".aa_composition_heatmap")

    def _plot_rscu_heatmap(self) -> None:
        """Plot pooled group RSCU for synonymous codons, excluding Met and Trp."""
        table = self.analyzer.group_codon_usage.copy()
        table = table.loc[table["Group"].isin(self.groups)]
        table = table.loc[~table["Abbr."].isin(["M", "W"])]
        if table.empty:
            return
        codon_order = (
            table[["Codon", "Abbr."]]
            .drop_duplicates()
            .sort_values(["Abbr.", "Codon"], kind="stable")["Codon"]
            .tolist()
        )
        pivot = table.pivot(index="Group", columns="Codon", values="RSCU")
        pivot = pivot.reindex(index=self.groups, columns=codon_order)
        matrix = pivot.to_numpy(dtype=float)
        finite = matrix[np.isfinite(matrix)]
        if finite.size == 0:
            return
        vmin = min(0.0, float(np.nanmin(finite)))
        vmax = max(2.0, float(np.nanmax(finite)))

        fig, ax = plt.subplots(
            figsize=(13.0, max(2.3, 1.2 + len(self.groups) * 0.62)),
            dpi=self.dpi,
            constrained_layout=True,
        )
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap="coolwarm",
            norm=TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax),
            interpolation="nearest",
        )
        ax.set_xticks(np.arange(len(codon_order)))
        ax.set_xticklabels(codon_order, rotation=90)
        ax.set_yticks(np.arange(len(self.groups)))
        ax.set_yticklabels([self._group_label(group).replace("\n", " ") for group in self.groups])
        ax.set_xlabel("Synonymous codon")
        ax.set_ylabel("")
        ax.tick_params(labelsize=max(self.font_size - 2.0, 1.0), length=2.0)
        cbar = fig.colorbar(image, ax=ax, fraction=0.018, pad=0.015)
        cbar.set_label("RSCU", fontsize=self.font_size)
        cbar.ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        self._save_figure(fig, ".rscu_heatmap")

    @staticmethod
    def _relative_ticks(minimum: int, maximum: int) -> list[int]:
        """Choose regular relative-position ticks while always retaining zero."""
        span = maximum - minimum
        if span >= 30:
            step = 10
        elif span >= 15:
            step = 5
        elif span >= 8:
            step = 2
        else:
            step = 1
        first = int(np.ceil(minimum / float(step)) * step)
        ticks = list(range(first, maximum + 1, step))
        ticks = sorted(set(ticks).union({0}))
        return [value for value in ticks if minimum <= value <= maximum]

    def _plot_position_profiles(self) -> None:
        """Plot position-resolved hydropathy and charged-residue balance around max sites."""
        table = self.analyzer.position_property_profile.copy()
        table = table.loc[table["Group"].isin(self.groups)]
        if table.empty:
            return
        fig, axes = plt.subplots(
            2,
            1,
            figsize=(7.0, 5.2),
            dpi=self.dpi,
            sharex=True,
            constrained_layout=True,
        )
        for group in self.groups:
            one = table.loc[table["Group"].eq(group)].sort_values("Relative_Position")
            if one.empty:
                continue
            x = one["Relative_Position"].to_numpy(dtype=float)
            color = self.group_colors[group]
            label = self._group_label(group).replace("\n", " ")
            axes[0].plot(
                x,
                one["Mean_Hydropathy"].to_numpy(dtype=float),
                linewidth=1.5,
                color=color,
                label=label,
            )
            axes[1].plot(
                x,
                one["Mean_Signed_Charge_Class"].to_numpy(dtype=float),
                linewidth=1.5,
                color=color,
                label=label,
            )
        for ax in axes:
            ax.axvline(0.0, color="grey", linewidth=0.7, linestyle="--")
            self._setup_axis(ax)
        axes[0].set_ylabel("Mean Kyte-Doolittle hydropathy")
        axes[1].set_ylabel("Mean signed charge class")
        axes[1].set_xlabel("Position relative to peak max site (aa)")
        minimum = int(table["Relative_Position"].min())
        maximum = int(table["Relative_Position"].max())
        axes[1].set_xticks(self._relative_ticks(minimum, maximum))
        axes[0].legend(
            frameon=False,
            fontsize=max(self.font_size - 1.0, 1.0),
            ncol=min(3, len(self.groups)),
        )
        self._save_figure(fig, ".local_position_profiles")

    def _plot_position_aa_enrichment(self) -> None:
        """Plot position-specific amino-acid enrichment around peak max sites."""
        table = self.analyzer.position_aa_profile.copy()
        table = table.loc[table["Group"].isin(self.groups)]
        if table.empty:
            return
        positions = list(
            range(
                int(table["Relative_Position"].min()),
                int(table["Relative_Position"].max()) + 1,
            )
        )
        matrices = []
        for group in self.groups:
            one = table.loc[table["Group"].eq(group)]
            pivot = one.pivot(
                index="AA",
                columns="Relative_Position",
                values="Log2_Enrichment_vs_All",
            )
            pivot = pivot.reindex(index=AA_ORDER, columns=positions).fillna(0.0)
            matrices.append(pivot.to_numpy(dtype=float))
        finite_values = np.concatenate(
            [matrix[np.isfinite(matrix)] for matrix in matrices if np.isfinite(matrix).any()]
        )
        if finite_values.size == 0:
            return
        vmax = max(float(np.nanpercentile(np.abs(finite_values), 98)), 0.5)

        fig, axes = plt.subplots(
            len(self.groups),
            1,
            figsize=(8.5, 1.7 + 2.2 * len(self.groups)),
            dpi=self.dpi,
            sharex=True,
            constrained_layout=True,
        )
        if len(self.groups) == 1:
            axes = [axes]
        image = None
        for ax, group, matrix in zip(axes, self.groups, matrices):
            image = ax.imshow(
                matrix,
                aspect="auto",
                cmap="coolwarm",
                norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
                interpolation="nearest",
                extent=[positions[0] - 0.5, positions[-1] + 0.5, len(AA_ORDER) - 0.5, -0.5],
            )
            ax.axvline(0.0, color="black", linewidth=0.6)
            ax.set_yticks(np.arange(len(AA_ORDER)))
            ax.set_yticklabels(AA_ORDER)
            ax.set_ylabel(self._group_label(group).replace("\n", " "))
            ax.tick_params(labelsize=max(self.font_size - 2.0, 1.0), length=2.0)
        axes[-1].set_xlabel("Position relative to peak max site (aa)")
        axes[-1].set_xticks(self._relative_ticks(positions[0], positions[-1]))
        if image is not None:
            cbar = fig.colorbar(image, ax=axes, fraction=0.018, pad=0.015)
            cbar.set_label("log2(AA frequency / pooled frequency)", fontsize=self.font_size)
            cbar.ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0))
        self._save_figure(fig, ".position_aa_enrichment")

    def draw_all(self) -> None:
        """Generate all default SeRP group sequence-property figures."""
        nonempty_groups = [
            group
            for group in self.groups
            if (
                group in set(self.analyzer.group_full_properties.get("Group", pd.Series(dtype=str)))
                or group in set(self.analyzer.peak_workflow.local_table_with_significant_groups().get("Plot_Group", pd.Series(dtype=str)))
            )
        ]
        if len(nonempty_groups) < 2:
            return
        self.groups = nonempty_groups
        self._plot_full_properties()
        self._plot_local_properties()
        self._plot_aa_composition_heatmap()
        self._plot_rscu_heatmap()
        self._plot_position_profiles()
        self._plot_position_aa_enrichment()