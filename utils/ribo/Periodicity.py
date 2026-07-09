#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.7.7
# Function: Provide RiboParser core functions for Periodicity analysis.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Periodicity summary table and periodicity figures.

"""3-nt periodicity analysis utilities for RiboParser.

This module keeps the original ``Periodicity`` workflow and delegates RPF
import to ``utils.ribo.RPFs.import_rpf``. The updated ``RPFs`` module handles
both current JSONL density files and legacy TXT density files, so downstream
periodicity analysis can keep one unified import path.
"""

from __future__ import annotations

import math
from collections import OrderedDict

from . import RPFs

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
FRAME_ORDER = ["0", "1", "2"]
FRAME_COLORS = {
    "0": "#E64B35",
    "1": "#4DBBD5",
    "2": "#00A087",
}


class Periodicity(object):
    """Calculate and plot 3-nt periodicity from RPF density files.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Periodicity``.
    """

    def __init__(self, args):
        # Input and output.
        self.transcript = args.transcript
        self.rpf = args.rpf
        self.output = args.output

        # Filter range around TIS/TTS and high-expression threshold.
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts

        # Imported RPF data.
        self.sample_name = None
        self.sample_num = 0
        self.raw_rpf = None
        self.merged_rpf = None
        self.total_rpf_num = None
        self.gene_rpf_sum = None
        self.high_gene = None
        self.high_rpf = None

        # Periodicity result.
        self.period = None

    def import_rpf(self):
        """Import RPF density data through the updated ``RPFs`` module.

        Returns
        -------
        None
            The function updates internal state in place.

        Notes
        -----
        The updated ``RPFs.import_rpf`` supports both current compact JSONL
        density files and legacy codon-level TXT files. The return structure is
        kept compatible with the historical RiboParser workflow.
        """
        rpf_results = RPFs.import_rpf(
            rpf_file=self.rpf,
            sites="P",
            frame="all",
            sample_num=None,
            sample_name=None,
            tis=self.tis,
            tts=self.tts,
            gene=self.transcript,
            rpf_num=self.rpf_num,
        )

        # Historical return order:
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num,
        # gene_rpf_sum, high_gene, high_rpf.
        self.raw_rpf = rpf_results[0].to_pandas()
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        self.high_rpf = rpf_results[7]

        del rpf_results

    @staticmethod
    def _frame_columns(raw_rpf: pd.DataFrame) -> list[str]:
        """Return frame columns from a codon-level RPF table."""
        frame_columns = []
        for column in raw_rpf.columns:
            if column.endswith("_f0") or column.endswith("_f1") or column.endswith("_f2"):
                frame_columns.append(column)
        return frame_columns

    @staticmethod
    def _column_to_sample_frame(column: str) -> tuple[str, str]:
        """Convert ``sample_f0`` to ``(sample, 0)``."""
        sample = column[:-3]
        frame = column[-1]
        return sample, frame

    def _get_periodicity_input(self) -> pd.DataFrame:
        """Return the filtered raw frame table used for periodicity."""
        if self.raw_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        period_rpf = self.raw_rpf

        # Apply the high-expression transcript filter calculated by RPFs.
        # This preserves the original RPFs-based filtering criterion while
        # still allowing frame-level periodicity to be calculated from raw_rpf.
        if self.high_gene is not None and len(self.high_gene) > 0:
            period_rpf = period_rpf.loc[period_rpf["name"].isin(self.high_gene), :].copy()

        if period_rpf.empty:
            raise ValueError("Filtered RPF density table is empty after high-expression filtering.")

        return period_rpf

    def calc_3nt_period(self):
        """Calculate 3-nt periodicity for each sample.

        Returns
        -------
        None
            The function updates ``self.period`` in place.
        """
        period_rpf = self._get_periodicity_input()
        frame_columns = self._frame_columns(period_rpf)

        if not frame_columns:
            raise ValueError("No frame density columns were found in the RPF density table.")

        records = []
        for column in frame_columns:
            sample, frame = self._column_to_sample_frame(column)
            count = period_rpf[column].sum()
            records.append([sample, frame, int(count)])

        period = pd.DataFrame(records, columns=["Sample", "Frame", "Count"])

        # Ensure all samples have frame 0, 1, and 2 rows.
        if self.sample_name is None:
            sample_names = period["Sample"].drop_duplicates().tolist()
        else:
            sample_names = list(self.sample_name)

        full_index = pd.MultiIndex.from_product(
            [sample_names, FRAME_ORDER],
            names=["Sample", "Frame"],
        )
        period = (
            period.set_index(["Sample", "Frame"])
            .reindex(full_index, fill_value=0)
            .reset_index()
        )

        period["Total"] = period.groupby("Sample")["Count"].transform("sum")
        period["Ratio"] = np.where(
            period["Total"] > 0,
            period["Count"] / period["Total"] * 100,
            0,
        )
        period.drop(columns=["Total"], inplace=True)

        self.period = period[["Sample", "Frame", "Count", "Ratio"]]

    def output_meta(self):
        """Write the periodicity summary table."""
        if self.period is None:
            raise ValueError("Periodicity has not been calculated yet.")

        out_txt = self.output + "_periodicity.txt"
        self.period.to_csv(out_txt, sep="\t", index=False)

    def get_sub_plot_num(self):
        """Calculate subplot rows and columns for per-sample panels."""
        factor1 = int(self.sample_num ** 0.5)
        factor2 = int(self.sample_num ** 0.5)
        square = factor1 * factor2

        while not square >= self.sample_num:
            if factor2 > factor1:
                factor1 += 1
            else:
                factor2 += 1
            square = factor1 * factor2

        return factor1, factor2

    def draw_3nt_period_count(self):
        """Draw one count barplot panel for each sample."""
        if self.period is None:
            raise ValueError("Periodicity has not been calculated yet.")

        out_pdf = self.output + "_count_periodicity_plot.pdf"
        out_png = self.output + "_count_periodicity_plot.png"

        nrow, ncol = self.get_sub_plot_num()
        fig, axes = plt.subplots(
            nrows=nrow,
            ncols=ncol,
            figsize=(ncol * 3, nrow * 3),
            sharey=True,
        )

        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        axes = axes.flatten()

        for (sample, group), ax in zip(self.period.groupby("Sample", sort=False), axes):
            group = group.copy()
            group["Frame"] = pd.Categorical(group["Frame"], categories=FRAME_ORDER, ordered=True)
            group.sort_values("Frame", inplace=True)
            frame_colors = [FRAME_COLORS[str(frame)] for frame in group["Frame"].astype(str)]
            bars = group.plot(
                x="Frame",
                y="Count",
                kind="bar",
                ax=ax,
                title=sample,
                legend=False,
                color=frame_colors,
                width=0.75,
            )
            ax.set_ylabel("Count")
            ax.set_xlabel("Frame")

            for bar in bars.patches:
                height = bar.get_height()
                ax.annotate(
                    f"{height:.0f}",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

        for ax in axes[len(self.period["Sample"].drop_duplicates()):]:
            ax.axis("off")

        plt.tight_layout()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png, dpi=300)
        plt.close(fig)

    def draw_3nt_period_ratio(self):
        """Draw one ratio barplot panel for each sample."""
        if self.period is None:
            raise ValueError("Periodicity has not been calculated yet.")

        out_pdf = self.output + "_ratio_periodicity_plot.pdf"
        out_png = self.output + "_ratio_periodicity_plot.png"

        nrow, ncol = self.get_sub_plot_num()
        fig, axes = plt.subplots(
            nrows=nrow,
            ncols=ncol,
            figsize=(ncol * 3, nrow * 3),
            sharey=True,
        )

        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        axes = axes.flatten()

        for (sample, group), ax in zip(self.period.groupby("Sample", sort=False), axes):
            group = group.copy()
            group["Frame"] = pd.Categorical(group["Frame"], categories=FRAME_ORDER, ordered=True)
            group.sort_values("Frame", inplace=True)
            frame_colors = [FRAME_COLORS[str(frame)] for frame in group["Frame"].astype(str)]
            bars = group.plot(
                x="Frame",
                y="Ratio",
                kind="bar",
                ax=ax,
                title=sample,
                legend=False,
                color=frame_colors,
                width=0.75,
            )
            ax.set_ylabel("Ratio (%)")
            ax.set_xlabel("Frame")
            ax.set_ylim(0, max(100, group["Ratio"].max() * 1.15))

            for bar in bars.patches:
                height = bar.get_height()
                ax.annotate(
                    f"{height:.0f}",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

        for ax in axes[len(self.period["Sample"].drop_duplicates()):]:
            ax.axis("off")

        plt.tight_layout()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png, dpi=300)
        plt.close(fig)

    def draw_3nt_period_stacked(self):
        """Draw a two-panel stacked barplot for count and ratio.

        The left panel shows RPF counts by sample, stacked by frame 0/1/2. The
        right panel shows the corresponding frame ratios by sample.
        """
        if self.period is None:
            raise ValueError("Periodicity has not been calculated yet.")

        out_pdf = self.output + "_stacked_periodicity_plot.pdf"
        out_png = self.output + "_stacked_periodicity_plot.png"

        count_df = self.period.pivot(index="Sample", columns="Frame", values="Count")
        ratio_df = self.period.pivot(index="Sample", columns="Frame", values="Ratio")

        count_df = count_df.reindex(columns=FRAME_ORDER, fill_value=0)
        ratio_df = ratio_df.reindex(columns=FRAME_ORDER, fill_value=0)

        sample_count = max(1, len(count_df.index))
        figure_width = min(max(8.5, sample_count * 0.22 + 5.5), 16.0)
        figure_height = 5.8
        colors = [FRAME_COLORS[frame] for frame in FRAME_ORDER]

        fig, axes = plt.subplots(
            nrows=1,
            ncols=2,
            figsize=(figure_width, figure_height),
            gridspec_kw={"wspace": 0.22},
        )

        count_df.plot(
            kind="bar",
            stacked=True,
            ax=axes[0],
            color=colors,
            width=0.82,
            edgecolor="none",
            legend=False,
        )
        axes[0].set_title("RPF count", fontsize=11)
        axes[0].set_xlabel("Sample")
        axes[0].set_ylabel("RPF count")
        axes[0].tick_params(axis="x", labelrotation=90, labelsize=8)
        axes[0].tick_params(axis="y", labelsize=8)
        for label in axes[0].get_xticklabels():
            label.set_ha("center")
        axes[0].spines["top"].set_visible(False)
        axes[0].spines["right"].set_visible(False)

        ratio_df.plot(
            kind="bar",
            stacked=True,
            ax=axes[1],
            color=colors,
            width=0.82,
            edgecolor="none",
            legend=False,
        )
        axes[1].set_title("RPF ratio", fontsize=11)
        axes[1].set_xlabel("Sample")
        axes[1].set_ylabel("RPF ratio (%)")
        axes[1].set_ylim(0, 100)
        axes[1].tick_params(axis="x", labelrotation=90, labelsize=8)
        axes[1].tick_params(axis="y", labelsize=8)
        for label in axes[1].get_xticklabels():
            label.set_ha("center")
        axes[1].spines["top"].set_visible(False)
        axes[1].spines["right"].set_visible(False)

        handles = [
            plt.Rectangle((0, 0), 1, 1, facecolor=FRAME_COLORS[frame], edgecolor="none")
            for frame in FRAME_ORDER
        ]
        labels = [f"Frame {frame}" for frame in FRAME_ORDER]
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=3,
            frameon=False,
            bbox_to_anchor=(0.5, -0.02),
            fontsize=9,
        )

        fig.subplots_adjust(bottom=0.28, left=0.08, right=0.98, top=0.90, wspace=0.22)
        fig.savefig(fname=out_pdf, bbox_inches="tight")
        fig.savefig(fname=out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
