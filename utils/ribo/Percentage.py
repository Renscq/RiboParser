#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate transcript-level CDS RPF coverage percentages.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Gene coverage percentage tables, summary tables, and figures.

"""Core functions for transcript-level RPF coverage analysis.

The reader delegates density import to :class:`utils.ribo.RPFs.RPFData`, so
compact JSONL files and legacy TXT density tables are processed through the same
workflow.

Coverage is calculated independently for each sample as::

    covered valid CDS codons / total valid CDS codons * 100

A valid CDS codon is a retained CDS row after common input filtering. A codon is
covered when its selected-frame P-site density is greater than zero. The minimum
RPF threshold is applied independently to each transcript and sample. Coverage
values that do not pass the sample-specific threshold are stored as missing
values rather than zero, preventing unexpressed transcripts from biasing sample
coverage distributions.
"""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import RPFs


RPM_SCALE = 1_000_000.0
BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
PLOT_COLOR = "#4C78A8"
GRID_ALPHA = 0.25


class Percentage:
    """Calculate transcript-level CDS RPF coverage percentages.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Percent``.
    """

    def __init__(self, args: Namespace):
        # Input and output.
        self.rpf = args.rpf
        self.transcript = getattr(args, "transcript", None)
        self.gene = self.transcript
        self.output = args.output

        # Density import and filtering.
        self.site = "P"
        self.frame = str(args.frame)
        self.rpf_num = float(args.min)
        self.norm = bool(args.normal)
        self.tis = int(getattr(args, "tis", 0))
        self.tts = int(getattr(args, "tts", 0))

        # Figure configuration.
        self.fig_format = str(
            getattr(args, "fig_format", "png")
        ).lower()
        self.dpi = int(getattr(args, "dpi", 300))
        self.font_size = float(getattr(args, "font_size", 11.0))

        # Imported RPF information.
        self.rpf_data: RPFs.RPFData | None = None
        self.file_format: str | None = None
        self.sample_name: list[str] = []
        self.sample_num = 0
        self.total_rpf_num = pd.Series(dtype=float)
        self.merged_rpf: pd.DataFrame | None = None
        self.high_rpf: pd.DataFrame | None = None

        # Transcript information.
        self.gene_table: pd.DataFrame | None = None
        self.high_gene = pd.Index([])
        self.gene_num = 0
        self.gene_num_by_sample = pd.Series(dtype="int64")

        # Calculation results.
        self.valid_codon_num = pd.Series(dtype="int64")
        self.covered_codon_num: pd.DataFrame | None = None
        self.gene_rpf_sum_raw: pd.DataFrame | None = None
        self.gene_rpf_sum: pd.DataFrame | None = None
        self.gene_abundance: pd.DataFrame | None = None
        self.expression_pass: pd.DataFrame | None = None
        self.gene_coverage: pd.DataFrame | None = None
        self.coverage_summary: pd.DataFrame | None = None

    # ------------------------------------------------------------------
    # Input
    # ------------------------------------------------------------------
    def _validate_density_table(self, table: pd.DataFrame) -> None:
        """Validate the imported codon-level density table."""
        if table.empty:
            raise ValueError(
                "The RPF density table is empty after input filtering."
            )

        required_columns = BASE_COLUMNS + self.sample_name
        missing_columns = [
            column
            for column in required_columns
            if column not in table.columns
        ]
        if missing_columns:
            raise ValueError(
                "RPF density table is missing required column(s): "
                + ", ".join(missing_columns)
            )
        if not self.sample_name:
            raise ValueError("No sample density columns were detected.")

    def read_rpf(self) -> None:
        """Import JSONL or TXT RPF density using the shared reader."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf,
            sample_name=None,
            gene=self.transcript,
            tis=self.tis,
            tts=self.tts,
        )
        self.file_format = str(self.rpf_data.file_format)
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf_num = (
            self.rpf_data.total_rpf_num
            .reindex(self.sample_name)
            .fillna(0.0)
            .astype(float)
        )

        merged_rpf = self.rpf_data.shifted_frame(
            sites=self.site,
            frame=self.frame,
        )
        self._validate_density_table(merged_rpf)

        for sample in self.sample_name:
            merged_rpf.loc[:, sample] = pd.to_numeric(
                merged_rpf[sample],
                errors="coerce",
            ).fillna(0.0)

        self.merged_rpf = merged_rpf

        print(
            "Imported RPF density: format={fmt}, samples={samples:,}, "
            "rows={rows:,}.".format(
                fmt=self.file_format,
                samples=self.sample_num,
                rows=len(self.merged_rpf),
            ),
            flush=True,
        )

    def import_gene(self) -> None:
        """Import the optional transcript annotation for compatibility.

        Notes
        -----
        Transcript filtering is already performed during ``RPFData.from_file``.
        This method retains the historical public API and stores the matching
        annotation rows for downstream code that accesses ``gene_table``.
        """
        if self.transcript is None:
            self.gene_table = None
            return

        gene_table = pd.read_csv(
            self.transcript,
            sep="\t",
            header=0,
            index_col=False,
        )
        if gene_table.empty:
            self.gene_table = gene_table
            return

        transcript_column = (
            "transcript_id"
            if "transcript_id" in gene_table.columns
            else gene_table.columns[0]
        )
        gene_table.loc[:, transcript_column] = (
            gene_table[transcript_column].astype(str)
        )

        if self.merged_rpf is not None:
            imported_transcripts = set(
                self.merged_rpf["name"].astype(str).unique()
            )
            gene_table = gene_table.loc[
                gene_table[transcript_column].isin(imported_transcripts)
            ].copy()

        self.gene_table = gene_table

    # ------------------------------------------------------------------
    # Calculation
    # ------------------------------------------------------------------
    def _calculate_abundance(
        self,
        gene_rpf_sum_raw: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return raw counts or RPM-normalized transcript abundance."""
        if not self.norm:
            return gene_rpf_sum_raw.astype(float)

        denominator = (
            self.total_rpf_num
            .reindex(self.sample_name)
            .replace(0.0, np.nan)
        )
        abundance = (
            gene_rpf_sum_raw.astype(float)
            .div(denominator, axis=1)
            .mul(RPM_SCALE)
        )

        zero_total_samples = denominator.index[denominator.isna()].tolist()
        if zero_total_samples:
            print(
                "Warning: zero total RPF count for sample(s): {samples}. "
                "Their RPM values were set to missing.".format(
                    samples=", ".join(zero_total_samples)
                ),
                flush=True,
            )
        return abundance

    def _build_summary(self) -> pd.DataFrame:
        """Build sample-level coverage and abundance summary statistics."""
        if (
            self.gene_coverage is None
            or self.gene_abundance is None
            or self.expression_pass is None
        ):
            raise ValueError("Coverage results are unavailable.")

        rows: list[dict[str, Any]] = []
        for sample in self.sample_name:
            pass_mask = self.expression_pass[sample]
            coverage = self.gene_coverage.loc[pass_mask, sample].dropna()
            abundance = self.gene_abundance.loc[pass_mask, sample].dropna()

            rows.append(
                {
                    "Sample": sample,
                    "RetainedTranscripts": int(pass_mask.sum()),
                    "MeanCoverage": (
                        float(coverage.mean())
                        if not coverage.empty
                        else np.nan
                    ),
                    "MedianCoverage": (
                        float(coverage.median())
                        if not coverage.empty
                        else np.nan
                    ),
                    "Q1Coverage": (
                        float(coverage.quantile(0.25))
                        if not coverage.empty
                        else np.nan
                    ),
                    "Q3Coverage": (
                        float(coverage.quantile(0.75))
                        if not coverage.empty
                        else np.nan
                    ),
                    "MeanAbundance": (
                        float(abundance.mean())
                        if not abundance.empty
                        else np.nan
                    ),
                    "MedianAbundance": (
                        float(abundance.median())
                        if not abundance.empty
                        else np.nan
                    ),
                }
            )

        summary = pd.DataFrame.from_records(rows)
        summary.insert(
            2,
            "AbundanceUnit",
            "RPM" if self.norm else "RPF count",
        )
        return summary

    def calc_density_percent(self) -> None:
        """Calculate sample-specific CDS RPF coverage percentages."""
        if self.merged_rpf is None:
            raise ValueError("RPF density has not been imported.")

        cds_mask = (
            self.merged_rpf["region"]
            .astype(str)
            .str.lower()
            .eq("cds")
        )
        cds_rpf = self.merged_rpf.loc[
            cds_mask,
            BASE_COLUMNS + self.sample_name,
        ].copy()
        if cds_rpf.empty:
            raise ValueError("No CDS codon remains after input filtering.")

        # The denominator is the actual number of retained, analyzable CDS
        # codons. This avoids the historical off-by-one error caused by using
        # max(from_tis) as the CDS length.
        valid_codon_num = (
            cds_rpf.groupby("name", sort=False)
            .size()
            .astype("int64")
        )
        valid_codon_num.name = "ValidCodonCount"

        gene_rpf_sum_raw = (
            cds_rpf.groupby("name", sort=False)[self.sample_name]
            .sum()
            .astype(float)
        )
        covered_codon_num = (
            cds_rpf.loc[:, self.sample_name]
            .gt(0)
            .groupby(cds_rpf["name"], sort=False)
            .sum()
            .astype("int64")
        )

        # Apply the expression threshold independently for every sample.
        expression_pass = gene_rpf_sum_raw.ge(self.rpf_num)
        coverage = (
            covered_codon_num
            .div(valid_codon_num, axis=0)
            .mul(100.0)
            .clip(lower=0.0, upper=100.0)
            .where(expression_pass)
        )

        retained_any = expression_pass.any(axis=1)
        if not retained_any.any():
            raise ValueError(
                "No transcript passed the sample-specific minimum CDS RPF "
                f"threshold: {self.rpf_num:g}."
            )

        retained_genes = expression_pass.index[retained_any]
        self.high_gene = retained_genes
        self.gene_num = int(len(retained_genes))
        self.gene_num_by_sample = (
            expression_pass.loc[retained_genes, :]
            .sum(axis=0)
            .astype("int64")
        )

        self.valid_codon_num = valid_codon_num.loc[retained_genes]
        self.covered_codon_num = covered_codon_num.loc[
            retained_genes,
            self.sample_name,
        ]
        self.gene_rpf_sum_raw = gene_rpf_sum_raw.loc[
            retained_genes,
            self.sample_name,
        ]
        self.gene_rpf_sum = self.gene_rpf_sum_raw
        self.gene_abundance = self._calculate_abundance(
            self.gene_rpf_sum_raw
        )
        self.expression_pass = expression_pass.loc[
            retained_genes,
            self.sample_name,
        ]
        self.gene_coverage = coverage.loc[
            retained_genes,
            self.sample_name,
        ]
        self.gene_coverage.index.name = "name"

        self.high_rpf = cds_rpf.loc[
            cds_rpf["name"].isin(retained_genes)
        ].copy()
        self.coverage_summary = self._build_summary()

        retained_text = ", ".join(
            f"{sample}={int(self.gene_num_by_sample[sample]):,}"
            for sample in self.sample_name
        )
        print(
            "Calculated CDS coverage for {genes:,} unique transcript(s); "
            "sample-specific retained transcripts: {retained}.".format(
                genes=self.gene_num,
                retained=retained_text,
            ),
            flush=True,
        )

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------
    @staticmethod
    def _adaptive_bin_number(
        values: np.ndarray,
        minimum: int = 20,
        maximum: int = 100,
    ) -> int:
        """Return a bounded Freedman-Diaconis histogram bin number."""
        values = values[np.isfinite(values)]
        if values.size < 2 or np.all(values == values[0]):
            return minimum

        edges = np.histogram_bin_edges(values, bins="fd")
        return int(np.clip(len(edges) - 1, minimum, maximum))

    def _save_figure(self, figure: plt.Figure, prefix: str) -> None:
        """Save one figure using the configured output format."""
        if self.fig_format in {"png", "both"}:
            figure.savefig(
                f"{prefix}.png",
                dpi=self.dpi,
                bbox_inches="tight",
            )
        if self.fig_format in {"pdf", "both"}:
            figure.savefig(
                f"{prefix}.pdf",
                bbox_inches="tight",
            )
        plt.close(figure)

    @staticmethod
    def _clean_axis(axis: plt.Axes) -> None:
        """Apply common publication-style axis formatting."""
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(
            axis="y",
            linewidth=0.4,
            alpha=GRID_ALPHA,
        )

    def _sample_plot_data(
        self,
        sample: str,
    ) -> tuple[pd.Series, pd.Series]:
        """Return coverage and abundance values passing the sample threshold."""
        if (
            self.gene_coverage is None
            or self.gene_abundance is None
            or self.expression_pass is None
        ):
            raise ValueError("Coverage results are unavailable.")

        pass_mask = self.expression_pass[sample]
        coverage = self.gene_coverage.loc[pass_mask, sample].dropna()
        abundance = self.gene_abundance.loc[pass_mask, sample].dropna()
        return coverage, abundance

    def draw_rpf_histogram(self) -> None:
        """Draw per-sample coverage and abundance histograms."""
        if self.gene_coverage is None:
            raise ValueError("Coverage percentages have not been calculated.")

        abundance_label = "RPM" if self.norm else "RPF count"

        for sample in self.sample_name:
            coverage, abundance = self._sample_plot_data(sample)
            if coverage.empty:
                print(
                    f"Skip histogram for {sample}: no transcript passed "
                    "the expression threshold.",
                    flush=True,
                )
                continue

            print(
                f"Draw coverage histogram of {sample}.",
                flush=True,
            )

            figure, axes = plt.subplots(
                1,
                2,
                figsize=(8.0, 3.6),
                dpi=150,
            )

            coverage_values = coverage.to_numpy(dtype=float)
            axes[0].hist(
                coverage_values,
                bins=self._adaptive_bin_number(coverage_values),
                range=(0.0, 100.0),
                color=PLOT_COLOR,
                alpha=0.80,
                edgecolor="white",
                linewidth=0.4,
            )
            axes[0].set_xlim(0.0, 100.0)
            axes[0].set_xlabel(
                "CDS codon coverage (%)",
                fontsize=self.font_size,
            )
            axes[0].set_ylabel(
                "Transcript count",
                fontsize=self.font_size,
            )
            axes[0].set_title(
                f"Coverage | n={len(coverage):,}",
                fontsize=self.font_size + 1.0,
                fontweight="normal",
            )
            self._clean_axis(axes[0])

            abundance_values = np.log2(
                abundance.to_numpy(dtype=float) + 1.0
            )
            axes[1].hist(
                abundance_values,
                bins=self._adaptive_bin_number(abundance_values),
                color=PLOT_COLOR,
                alpha=0.80,
                edgecolor="white",
                linewidth=0.4,
            )
            axes[1].set_xlabel(
                f"log2({abundance_label} + 1)",
                fontsize=self.font_size,
            )
            axes[1].set_ylabel(
                "Transcript count",
                fontsize=self.font_size,
            )
            axes[1].set_title(
                f"Abundance | n={len(abundance):,}",
                fontsize=self.font_size + 1.0,
                fontweight="normal",
            )
            self._clean_axis(axes[1])

            figure.suptitle(
                sample,
                fontsize=self.font_size + 2.0,
                fontweight="normal",
            )
            figure.tight_layout()
            self._save_figure(
                figure,
                f"{self.output}_{sample}_coverage_histogram",
            )

    def _boxplot_values(
        self,
    ) -> tuple[list[str], list[np.ndarray], list[np.ndarray]]:
        """Return aligned sample values for coverage and abundance boxplots."""
        labels: list[str] = []
        coverage_values: list[np.ndarray] = []
        abundance_values: list[np.ndarray] = []

        for sample in self.sample_name:
            coverage, abundance = self._sample_plot_data(sample)
            if coverage.empty:
                continue
            labels.append(sample)
            coverage_values.append(coverage.to_numpy(dtype=float))
            abundance_values.append(
                np.log2(abundance.to_numpy(dtype=float) + 1.0)
            )

        return labels, coverage_values, abundance_values

    def draw_rpf_boxplot(self) -> None:
        """Draw sample-level coverage and abundance boxplots."""
        labels, coverage_values, abundance_values = self._boxplot_values()
        if not labels:
            raise ValueError(
                "No sample has transcripts available for boxplot generation."
            )

        print(
            "Draw sample-level RPF coverage boxplots.",
            flush=True,
        )

        figure_width = max(
            7.5,
            min(18.0, 5.5 + 0.55 * len(labels)),
        )
        figure, axes = plt.subplots(
            1,
            2,
            figsize=(figure_width, 4.8),
            dpi=150,
        )

        coverage_box = axes[0].boxplot(
            coverage_values,
            tick_labels=labels,
            patch_artist=True,
            showfliers=False,
            widths=0.62,
            medianprops={"linewidth": 1.0},
            boxprops={"linewidth": 0.6},
            whiskerprops={"linewidth": 0.6},
            capprops={"linewidth": 0.6},
        )
        abundance_box = axes[1].boxplot(
            abundance_values,
            tick_labels=labels,
            patch_artist=True,
            showfliers=False,
            widths=0.62,
            medianprops={"linewidth": 1.0},
            boxprops={"linewidth": 0.6},
            whiskerprops={"linewidth": 0.6},
            capprops={"linewidth": 0.6},
        )

        for boxplot in (coverage_box, abundance_box):
            for patch in boxplot["boxes"]:
                patch.set_facecolor(PLOT_COLOR)
                patch.set_alpha(0.75)

        rotation = 0 if len(labels) <= 6 else 45
        horizontal_alignment = "center" if rotation == 0 else "right"

        axes[0].set_ylim(0.0, 100.0)
        axes[0].set_ylabel(
            "CDS codon coverage (%)",
            fontsize=self.font_size,
        )
        axes[0].set_xlabel("Sample", fontsize=self.font_size)
        axes[0].set_title(
            "Coverage distribution",
            fontsize=self.font_size + 1.0,
            fontweight="normal",
        )

        abundance_label = "RPM" if self.norm else "RPF count"
        axes[1].set_ylabel(
            f"log2({abundance_label} + 1)",
            fontsize=self.font_size,
        )
        axes[1].set_xlabel("Sample", fontsize=self.font_size)
        axes[1].set_title(
            "Abundance distribution",
            fontsize=self.font_size + 1.0,
            fontweight="normal",
        )

        for axis in axes:
            axis.tick_params(
                axis="x",
                labelrotation=rotation,
                labelsize=max(self.font_size - 1.0, 1.0),
            )
            for label in axis.get_xticklabels():
                label.set_horizontalalignment(horizontal_alignment)
            self._clean_axis(axis)

        figure.tight_layout()
        self._save_figure(
            figure,
            f"{self.output}_coverage_boxplot",
        )

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    def output_density_percent(self) -> None:
        """Output transcript coverage and sample summary tables."""
        if self.gene_coverage is None or self.coverage_summary is None:
            raise ValueError("Coverage percentages have not been calculated.")

        coverage_file = f"{self.output}_gene_coverage_percent.txt"
        summary_file = f"{self.output}_gene_coverage_summary.txt"

        self.gene_coverage.round(6).to_csv(
            coverage_file,
            sep="\t",
            index=True,
            header=True,
            na_rep="NA",
        )
        self.coverage_summary.round(6).to_csv(
            summary_file,
            sep="\t",
            index=False,
            header=True,
            na_rep="NA",
        )

        print(f"Output gene coverage: {coverage_file}", flush=True)
        print(f"Output coverage summary: {summary_file}", flush=True)
