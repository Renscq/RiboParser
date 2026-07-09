#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.007
# Function: Provide core functions for RPF sample correlation analysis.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Gene-level and codon-level correlation tables, heatmaps, and summary JSON.

"""Core correlation utilities for RiboParser RPF density files.

This module contains the reusable implementation used by ``rpf_Corr``. The
command-line parser is intentionally kept in ``rpf_Corr.py`` so that the
correlation workflow can be reused by downstream scripts and tests.
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
import seaborn as sns


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
FRAME_LABELS = OrderedDict(
    [
        ("frame", "All frames"),
        ("f0", "Frame 0"),
        ("f1", "Frame 1"),
        ("f2", "Frame 2"),
    ]
)
FRAME_TO_RPF_ARG = {
    "frame": "all",
    "f0": "0",
    "f1": "1",
    "f2": "2",
}
REGION_CHOICES = ["all", "5utr", "cds", "3utr"]
RPM_SCALE = 1_000_000.0


class RPFCorrelation(object):
    """Calculate and plot sample correlations from RPF density tables.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Corr``.
    """

    def __init__(self, args: Namespace):
        self.rpf = args.rpf
        self.output = args.output
        self.transcript = args.transcript

        self.level = args.level
        self.method = args.method
        self.normal = bool(args.normal)
        self.normal_scale = RPM_SCALE

        self.region = args.region
        self.min_gene_count = float(args.min_gene_count)

        self.cmap = args.cmap
        self.figure_width = args.figure_width
        self.figure_height = args.figure_height

        self.rpf_data = None
        self.sample_name: list[str] = []
        self.sample_num = 0
        self.total_rpf_num = None
        self.file_format = None

        self.raw_tables: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self.density_tables: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self.gene_matrices: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self.rpf_matrices: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self.gene_corr: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self.rpf_corr: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self.summary_records: list[dict[str, Any]] = []
        self.output_files: OrderedDict[str, str] = OrderedDict()

    # ------------------------------------------------------------------
    # Input helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_key(key: str) -> str:
        """Return a short output-safe frame key."""
        if key == "frame":
            return "frame"
        return key

    def import_rpf(self) -> None:
        """Import RPF density data through the updated RPFData reader."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf,
            sample_name=None,
            gene=self.transcript,
            tis=None,
            tts=None,
        )
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf_num = self.rpf_data.total_rpf_num
        self.file_format = self.rpf_data.file_format

        if self.sample_num < 2:
            raise ValueError("At least two samples are required for correlation analysis.")

        print(
            "Imported RPF density for {sample_num} sample(s), format={fmt}.".format(
                sample_num=self.sample_num,
                fmt=self.file_format,
            ),
            flush=True,
        )

    def _sample_total(self, sample: str, raw_table: pd.DataFrame) -> float:
        """Return sample total RPF count for RPM normalization."""
        if self.total_rpf_num is not None and sample in self.total_rpf_num.index:
            return float(self.total_rpf_num.loc[sample])
        if sample in raw_table.columns:
            return float(raw_table[sample].sum())
        return 0.0

    def _normalize_table(self, raw_table: pd.DataFrame) -> pd.DataFrame:
        """Normalize sample columns to RPM if requested."""
        density_table = raw_table.copy()
        for sample in self.sample_name:
            density_table[sample] = (
                pd.to_numeric(density_table[sample], errors="coerce")
                .fillna(0.0)
                .astype("float64")
            )

        if not self.normal:
            return density_table

        for sample in self.sample_name:
            total = self._sample_total(sample, raw_table)
            if total > 0:
                density_table[sample] = density_table[sample] * self.normal_scale / total
            else:
                density_table[sample] = 0.0
        return density_table

    def build_frame_tables(self) -> None:
        """Build all-frame and frame-specific sample density tables."""
        if self.rpf_data is None:
            raise ValueError("RPF data has not been imported yet.")

        for frame_key, rpf_frame in FRAME_TO_RPF_ARG.items():
            print("Build density table for {label}.".format(label=FRAME_LABELS[frame_key]), flush=True)
            raw_table = self.rpf_data.shifted_frame(sites="P", frame=rpf_frame)
            self._validate_density_table(raw_table, frame_key)
            self.raw_tables[frame_key] = raw_table
            self.density_tables[frame_key] = self._normalize_table(raw_table)

    def _validate_density_table(self, density_table: pd.DataFrame, frame_key: str) -> None:
        """Validate required columns in one frame table."""
        missing_base = [column for column in BASE_COLUMNS if column not in density_table.columns]
        if missing_base:
            raise ValueError(
                "Density table for {frame} lacks required base column(s): {columns}".format(
                    frame=frame_key,
                    columns=", ".join(missing_base),
                )
            )

        missing_samples = [sample for sample in self.sample_name if sample not in density_table.columns]
        if missing_samples:
            raise ValueError(
                "Density table for {frame} lacks sample column(s): {columns}".format(
                    frame=frame_key,
                    columns=", ".join(missing_samples),
                )
            )

    # ------------------------------------------------------------------
    # Filtering and matrix construction
    # ------------------------------------------------------------------

    @staticmethod
    def _filter_region(table: pd.DataFrame, region: str) -> pd.DataFrame:
        """Filter a density table by transcript region."""
        if region == "all":
            return table
        return table.loc[table["region"] == region, :]

    def _filter_feature_matrix(
        self,
        raw_matrix: pd.DataFrame,
        expr_matrix: pd.DataFrame,
        level: str,
        frame_key: str,
    ) -> pd.DataFrame:
        """Filter low-information features before correlation calculation."""
        if raw_matrix.empty:
            return expr_matrix.iloc[0:0, :]

        keep = pd.Series(True, index=raw_matrix.index)

        if level == "gene" and self.min_gene_count > 0:
            keep &= raw_matrix.sum(axis=1) >= self.min_gene_count
        elif level == "rpf":
            keep &= raw_matrix.sum(axis=1) > 0
        elif level != "gene":
            raise ValueError("level must be 'gene' or 'rpf'.")

        keep &= expr_matrix.var(axis=1) > 0
        filtered = expr_matrix.loc[keep, :].copy()
        self.summary_records.append(
            {
                "Level": level,
                "Frame": frame_key,
                "Region": self.region,
                "InputFeatureCount": int(raw_matrix.shape[0]),
                "FilteredFeatureCount": int(filtered.shape[0]),
                "Method": self.method,
            }
        )
        return filtered

    def _build_gene_matrix(self, frame_key: str) -> pd.DataFrame:
        """Build a transcript-level expression matrix for one frame."""
        raw_table = self._filter_region(self.raw_tables[frame_key], self.region)
        expr_table = self._filter_region(self.density_tables[frame_key], self.region)

        raw_matrix = raw_table.groupby("name", sort=False)[self.sample_name].sum()
        expr_matrix = expr_table.groupby("name", sort=False)[self.sample_name].sum()
        expr_matrix = expr_matrix.astype("float64").replace([np.inf, -np.inf], np.nan).fillna(0.0)
        return self._filter_feature_matrix(raw_matrix, expr_matrix, "gene", frame_key)

    def _build_rpf_matrix(self, frame_key: str) -> pd.DataFrame:
        """Build a codon-position expression matrix for one frame."""
        raw_table = self._filter_region(self.raw_tables[frame_key], self.region)
        expr_table = self._filter_region(self.density_tables[frame_key], self.region)

        raw_matrix = raw_table.loc[:, self.sample_name].copy()
        expr_matrix = expr_table.loc[:, self.sample_name].copy()
        expr_matrix = expr_matrix.astype("float64").replace([np.inf, -np.inf], np.nan).fillna(0.0)
        return self._filter_feature_matrix(raw_matrix, expr_matrix, "rpf", frame_key)

    def calculate_correlations(self) -> None:
        """Calculate gene-level and/or codon-level sample correlations."""
        run_gene = self.level in {"gene", "both"}
        run_rpf = self.level in {"rpf", "both"}

        for frame_key in FRAME_LABELS.keys():
            if run_gene:
                gene_matrix = self._build_gene_matrix(frame_key)
                self.gene_matrices[frame_key] = gene_matrix
                self.gene_corr[frame_key] = self._calculate_corr(gene_matrix, "gene", frame_key)

            if run_rpf:
                rpf_matrix = self._build_rpf_matrix(frame_key)
                self.rpf_matrices[frame_key] = rpf_matrix
                self.rpf_corr[frame_key] = self._calculate_corr(rpf_matrix, "rpf", frame_key)

    def _calculate_corr(self, matrix: pd.DataFrame, level: str, frame_key: str) -> pd.DataFrame:
        """Calculate a sample correlation matrix from a feature-by-sample matrix."""
        if matrix.shape[0] < 2:
            print(
                "Warning: {level} {frame} has fewer than two filtered features; correlation may be NA.".format(
                    level=level,
                    frame=frame_key,
                ),
                flush=True,
            )

        corr = matrix.corr(method=self.method)
        corr = corr.reindex(index=self.sample_name, columns=self.sample_name)
        return corr

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def output_tables(self) -> None:
        """Write correlation and feature matrices."""
        for level, corr_dict, matrix_dict in (
            ("gene", self.gene_corr, self.gene_matrices),
            ("rpf", self.rpf_corr, self.rpf_matrices),
        ):
            if not corr_dict:
                continue

            for frame_key, corr in corr_dict.items():
                out_key = self._safe_key(frame_key)
                corr_file = f"{self.output}_{level}_corr_{out_key}.txt"
                corr.to_csv(corr_file, sep="\t", index=True)
                self.output_files[f"{level}_{out_key}_corr"] = corr_file


        summary_txt = self.output + "_corr.feature_summary.txt"
        pd.DataFrame.from_records(self.summary_records).to_csv(summary_txt, sep="\t", index=False)
        self.output_files["feature_summary"] = summary_txt

    def write_summary(self) -> None:
        """Write a machine-readable correlation analysis summary."""
        summary_json = self.output + "_corr.summary.json"
        self.output_files["summary_json"] = summary_json
        summary = OrderedDict(
            [
                ("tool", "rpf_Corr"),
                ("version", "0.2.8-dev.007"),
                ("input_rpf", os.path.abspath(self.rpf)),
                ("input_format", self.file_format),
                ("transcript_filter", os.path.abspath(self.transcript) if self.transcript else None),
                ("sample_count", self.sample_num),
                ("samples", self.sample_name),
                (
                    "parameters",
                    OrderedDict(
                        [
                            ("level", self.level),
                            ("method", self.method),
                            ("normal", self.normal),
                            ("region", self.region),
                            ("min_gene_count", self.min_gene_count),
                            ("cmap", self.cmap),
                        ]
                    ),
                ),
                ("feature_summary", self.summary_records),
                ("output_files", self.output_files),
            ]
        )
        with open(summary_json, "w", encoding="utf-8") as out:
            json.dump(summary, out, ensure_ascii=False, indent=2)
            out.write("\n")

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    def _panel_order(self, corr_dict: OrderedDict[str, pd.DataFrame]) -> list[str]:
        """Return the original input sample order for all panels."""
        return [sample for sample in self.sample_name if any(sample in corr.index for corr in corr_dict.values())]

    @staticmethod
    def _auto_color_limits(corr_dict: OrderedDict[str, pd.DataFrame], sample_order: list[str]) -> tuple[float, float]:
        """Return automatic color limits for one analysis level.

        The diagonal is excluded when possible because it is always 1 and can
        compress the useful color range of sample-to-sample correlations.
        """
        values = []
        for corr in corr_dict.values():
            if corr is None or corr.empty:
                continue
            corr = corr.reindex(index=sample_order, columns=sample_order).astype(float)
            arr = corr.to_numpy(dtype=float)
            if arr.shape[0] > 1:
                mask = ~np.eye(arr.shape[0], dtype=bool)
                arr = arr[mask]
            arr = arr[np.isfinite(arr)]
            if arr.size:
                values.append(arr)

        if not values:
            return 0.0, 1.0

        merged = np.concatenate(values)
        if merged.size == 0:
            return 0.0, 1.0

        vmin = float(np.nanmin(merged))
        vmax = float(np.nanmax(merged))

        if not np.isfinite(vmin) or not np.isfinite(vmax):
            return 0.0, 1.0

        vmin = max(-1.0, vmin)
        vmax = min(1.0, vmax)

        if np.isclose(vmin, vmax):
            delta = max(0.02, abs(vmax) * 0.02)
            vmin = max(-1.0, vmin - delta)
            vmax = min(1.0, vmax + delta)

        return vmin, vmax

    def _plot_one_heatmap(
        self,
        ax,
        corr: pd.DataFrame,
        title: str,
        sample_order: list[str],
        show_y_label: bool,
        vmin: float,
        vmax: float,
        cbar_ax=None,
    ):
        """Draw one correlation heatmap panel."""
        corr = corr.reindex(index=sample_order, columns=sample_order)

        heatmap = sns.heatmap(
            corr,
            ax=ax,
            cmap=self.cmap,
            vmin=vmin,
            vmax=vmax,
            square=True,
            linewidths=0.35,
            linecolor="white",
            annot=False,
            cbar=cbar_ax is not None,
            cbar_ax=cbar_ax,
            cbar_kws={"orientation": "horizontal", "label": f"{self.method} correlation"} if cbar_ax is not None else None,
        )
        ax.set_title(title, fontsize=13)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="both", length=0)

        if corr.shape[0] <= 16:
            label_size = 9
        elif corr.shape[0] <= 30:
            label_size = 7
        else:
            label_size = 5

        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=90,
            ha="center",
            va="top",
            fontsize=label_size,
        )
        if show_y_label:
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=label_size)
        else:
            ax.set_yticklabels([])

        if not show_y_label:
            ax.set_yticklabels([])

        return heatmap

    def draw_corr_plot(self, level: str, corr_dict: OrderedDict[str, pd.DataFrame]) -> None:
        """Draw a 2x2 correlation heatmap figure for one analysis level."""
        if not corr_dict:
            return

        sample_order = self._panel_order(corr_dict)
        sample_num = len(sample_order)
        width = self.figure_width or min(max(16.0, sample_num * 0.95 + 7.5), 40.0)
        height = self.figure_height or min(max(15.0, sample_num * 0.90 + 7.5), 42.0)
        color_vmin, color_vmax = self._auto_color_limits(corr_dict, sample_order)

        out_pdf = f"{self.output}_{level}_correlation_plot.pdf"
        out_png = f"{self.output}_{level}_correlation_plot.png"

        fig = plt.figure(figsize=(width, height))
        grid = fig.add_gridspec(
            nrows=3,
            ncols=2,
            height_ratios=[1, 1, 0.055],
            hspace=0.32,
            wspace=0.18,
        )
        axes = [
            fig.add_subplot(grid[0, 0]),
            fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, 0]),
            fig.add_subplot(grid[1, 1]),
        ]
        cbar_ax = fig.add_subplot(grid[2, :])

        titles = {
            "frame": "All frames",
            "f0": "Frame 0",
            "f1": "Frame 1",
            "f2": "Frame 2",
        }
        value_label = "RPM" if self.normal else "raw count"
        suptitle = "{level}-level sample correlation ({region}, {value})".format(
            level="Gene" if level == "gene" else "Codon",
            region=self.region,
            value=value_label,
        )
        fig.suptitle(suptitle, fontsize=15, y=0.988)

        for idx, frame_key in enumerate(FRAME_LABELS.keys()):
            corr = corr_dict.get(frame_key)
            if corr is None:
                axes[idx].set_axis_off()
                continue

            self._plot_one_heatmap(
                ax=axes[idx],
                corr=corr,
                title=titles[frame_key],
                sample_order=sample_order,
                show_y_label=idx in {0, 2},
                vmin=color_vmin,
                vmax=color_vmax,
                cbar_ax=cbar_ax if idx == 3 else None,
            )

        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        self.output_files[f"{level}_correlation_pdf"] = out_pdf
        self.output_files[f"{level}_correlation_png"] = out_png

    # ------------------------------------------------------------------
    # Public workflow
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Run the complete correlation analysis workflow."""
        self.import_rpf()
        self.build_frame_tables()
        self.calculate_correlations()
        self.output_tables()

        if self.gene_corr:
            self.draw_corr_plot("gene", self.gene_corr)
        if self.rpf_corr:
            self.draw_corr_plot("rpf", self.rpf_corr)

        self.write_summary()
