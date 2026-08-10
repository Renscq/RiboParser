#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev003
# Function: Draw compact SeRP enrichment profiles for selected genes or transcripts.
# Input: SeRP peak table, per-codon enrichment profile, and target gene/transcript identifiers.
# Output: Selected transcript-level SeRP peak figures with enrichment and core peak regions.

"""Targeted plotting workflow for selective ribosome profiling peaks."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


ENRICHMENT_COLOR = "#0072B2"
CORE_PEAK_COLOR = ENRICHMENT_COLOR
UTR_BACKGROUND_COLOR = "#F2F2F2"
REFERENCE_COLOR = "#7F7F7F"


class SeRPPlot:
    """Draw SeRP peak profiles only for explicitly selected targets."""

    PROFILE_COLUMNS = ["name", "from_tis", "region", "enrich"]
    PEAK_COLUMNS = {
        "transcripts",
        "gene_name",
        "peak_num",
        "peak_start",
        "peak_end",
    }

    def __init__(self, args):
        """Initialize the plotting workflow from parsed arguments."""
        self.input_prefix = args.input
        self.peak_file = self.input_prefix + "_peaks.log"
        self.profile_file = self.input_prefix + "_peaks_ratio.txt"
        self.output_prefix = args.output if args.output else self.input_prefix
        self.target = args.target
        self.target_list = args.target_list
        self.threshold = args.threshold
        self.show_max_site = args.show_max_site
        self.shade_utr = args.shade_utr
        self.output_format = args.output_format
        self.dpi = args.dpi
        self.font_size = args.font_size
        self.y_max = args.y_max

        self.peak_table = pd.DataFrame()
        self.profile_table = pd.DataFrame()
        self.target_ids: set[str] = set()
        self.target_transcripts: set[str] = set()
        self.figure_dir = Path(self.output_prefix + "_figures")

    @staticmethod
    def _safe_file_name(value: str) -> str:
        """Convert a transcript identifier into a filesystem-safe name."""
        safe = str(value)
        for old, new in [("/", "_"), ("\\", "_"), (":", "_"), (" ", "_")]:
            safe = safe.replace(old, new)
        return safe

    @staticmethod
    def _read_target_file(path: str) -> set[str]:
        """Read targets from the first non-empty column of a text table."""
        table = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
        if table.empty:
            return set()

        for column in table.columns:
            values = table[column].dropna().astype(str).str.strip()
            values = values.loc[values != ""]
            if not values.empty:
                return set(values.tolist())
        return set()

    def read_targets(self) -> None:
        """Read direct or file-based target identifiers."""
        if self.target:
            self.target_ids = {self.target.strip()}
        else:
            self.target_ids = self._read_target_file(self.target_list)

        if not self.target_ids:
            raise ValueError("No valid target gene or transcript was found.")

    def read_peak_table(self) -> None:
        """Read peak results and resolve gene names to transcript identifiers."""
        self.peak_table = pd.read_csv(self.peak_file, sep="\t", dtype=str)
        missing = sorted(self.PEAK_COLUMNS.difference(self.peak_table.columns))
        if missing:
            raise ValueError(
                "Peak table is missing required column(s): {columns}".format(
                    columns=", ".join(missing)
                )
            )

        transcripts = set(self.peak_table["transcripts"].dropna().astype(str))
        resolved = set()
        for target in self.target_ids:
            if target in transcripts:
                resolved.add(target)

            matched = self.peak_table.loc[
                self.peak_table["gene_name"].astype(str) == target,
                "transcripts",
            ].dropna()
            resolved.update(matched.astype(str).tolist())

            if target not in transcripts and matched.empty:
                resolved.add(target)

        self.target_transcripts = resolved

    def read_profiles(self, chunksize: int = 500_000) -> None:
        """Stream the profile table and retain only requested transcripts."""
        header = pd.read_csv(self.profile_file, sep="\t", nrows=0)
        missing = [
            column for column in self.PROFILE_COLUMNS if column not in header.columns
        ]
        if missing:
            raise ValueError(
                "Profile table is missing required column(s): {columns}".format(
                    columns=", ".join(missing)
                )
            )

        selected_chunks = []
        for chunk in pd.read_csv(
            self.profile_file,
            sep="\t",
            usecols=self.PROFILE_COLUMNS,
            chunksize=chunksize,
        ):
            selected = chunk.loc[
                chunk["name"].astype(str).isin(self.target_transcripts)
            ].copy()
            if not selected.empty:
                selected_chunks.append(selected)

        if not selected_chunks:
            raise ValueError(
                "None of the selected targets were found in the SeRP profile table: "
                "{targets}".format(targets=", ".join(sorted(self.target_transcripts)))
            )

        self.profile_table = pd.concat(selected_chunks, axis=0, ignore_index=True)
        found = set(self.profile_table["name"].astype(str).unique())
        missing_targets = sorted(self.target_transcripts.difference(found))
        if missing_targets:
            print(
                "Warning: no profile was available for target(s): {targets}".format(
                    targets=", ".join(missing_targets)
                ),
                flush=True,
            )

    def _gene_name(self, transcript: str) -> str:
        """Return the annotated gene name for one transcript when available."""
        matched = self.peak_table.loc[
            self.peak_table["transcripts"].astype(str) == transcript,
            "gene_name",
        ].dropna()
        if matched.empty:
            return "-"
        return str(matched.iloc[0])

    def _peak_regions(self, transcript: str) -> list[dict[str, float | int]]:
        """Return called core peak regions for one transcript."""
        table = self.peak_table.loc[
            self.peak_table["transcripts"].astype(str) == transcript
        ].copy()
        if table.empty:
            return []

        table["peak_num_numeric"] = pd.to_numeric(table["peak_num"], errors="coerce")
        table["peak_start_numeric"] = pd.to_numeric(
            table["peak_start"], errors="coerce"
        )
        table["peak_end_numeric"] = pd.to_numeric(table["peak_end"], errors="coerce")
        table = table.loc[
            table["peak_num_numeric"].gt(0)
            & table["peak_start_numeric"].notna()
            & table["peak_end_numeric"].notna()
        ].copy()
        if table.empty:
            return []

        if "max_site" in table.columns:
            table["max_site_numeric"] = pd.to_numeric(
                table["max_site"], errors="coerce"
            )

        regions: list[dict[str, float | int]] = []
        table = table.sort_values(
            ["peak_start_numeric", "peak_end_numeric"], kind="stable"
        )
        for _, row in table.iterrows():
            core_start = int(row["peak_start_numeric"])
            core_end = int(row["peak_end_numeric"])

            max_site = np.nan
            if "max_site_numeric" in row.index and pd.notna(row["max_site_numeric"]):
                max_site = float(row["max_site_numeric"])

            regions.append(
                {
                    "peak_num": int(row["peak_num_numeric"]),
                    "core_start": core_start,
                    "core_end": core_end,
                    "max_site": max_site,
                }
            )
        return regions

    def _save_figure(self, fig, output_stem: Path) -> None:
        """Save one figure in the requested output format."""
        if self.output_format in {"png", "both"}:
            fig.savefig(
                output_stem.with_suffix(".png"),
                dpi=self.dpi,
                bbox_inches="tight",
            )
        if self.output_format in {"pdf", "both"}:
            fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")

    @staticmethod
    def _shade_utr_regions(ax, table: pd.DataFrame) -> None:
        """Shade 5-prime and 3-prime UTRs with a subtle background."""
        for region_name in ["5utr", "3utr"]:
            positions = table.loc[
                table["region"].astype(str).str.lower().eq(region_name),
                "from_tis",
            ]
            if positions.empty:
                continue
            start = int(positions.min())
            end = int(positions.max())
            ax.axvspan(
                start - 0.5,
                end + 0.5,
                facecolor=UTR_BACKGROUND_COLOR,
                edgecolor="none",
                zorder=-4,
            )

    def _draw_peak_regions(self, ax, regions: list[dict[str, float | int]]) -> None:
        """Shade called core peak intervals behind the enrichment curve."""
        label_used = False
        for region in regions:
            core_start = int(region["core_start"])
            core_end = int(region["core_end"])
            ax.axvspan(
                core_start - 0.5,
                core_end + 0.5,
                facecolor=CORE_PEAK_COLOR,
                edgecolor="none",
                alpha=0.12,
                zorder=-1,
                label=None if label_used else "Core peak",
            )
            label_used = True

    def _mark_peak_maxima(
        self,
        ax,
        table: pd.DataFrame,
        regions: list[dict[str, float | int]],
    ) -> None:
        """Mark the reported maximum-enrichment site for each called peak."""
        if not self.show_max_site or not regions:
            return

        indexed = table.set_index("from_tis", drop=False)
        for region in regions:
            max_site = region["max_site"]
            if not np.isfinite(max_site):
                continue
            position = int(max_site)
            if position not in indexed.index:
                continue

            value = indexed.loc[position, "enrich"]
            if isinstance(value, pd.Series):
                value = value.iloc[0]
            value = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            if not np.isfinite(value):
                continue

            ax.scatter(
                [position],
                [float(value)],
                s=18,
                facecolor=CORE_PEAK_COLOR,
                edgecolor="white",
                linewidth=0.6,
                zorder=5,
            )
            ax.annotate(
                "P{0}".format(int(region["peak_num"])),
                xy=(position, float(value)),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=max(self.font_size - 2.0, 6.0),
                color=CORE_PEAK_COLOR,
            )

    def _draw_transcript(self, transcript: str, table: pd.DataFrame) -> None:
        """Draw one transcript-level enrichment profile with called core peaks."""
        table = table.copy()
        table["from_tis"] = pd.to_numeric(
            table["from_tis"], errors="raise"
        ).astype(int)
        table["enrich"] = pd.to_numeric(table["enrich"], errors="coerce")
        table = table.sort_values("from_tis", kind="stable")

        x = table["from_tis"].to_numpy(dtype=int)
        enrich = table["enrich"].to_numpy(dtype=float)
        finite_enrich = enrich[np.isfinite(enrich)]
        if finite_enrich.size == 0:
            raise ValueError(
                "No finite enrichment values were available for transcript: "
                "{0}".format(transcript)
            )

        regions = self._peak_regions(transcript)
        cds_positions = table.loc[
            table["region"].astype(str).str.lower().eq("cds"), "from_tis"
        ].to_numpy(dtype=int)

        mpl.rcParams["axes.linewidth"] = 0.7
        fig, ax = plt.subplots(figsize=(6.0, 3.35), dpi=self.dpi)

        if self.shade_utr:
            self._shade_utr_regions(ax, table)
        self._draw_peak_regions(ax, regions)

        ax.plot(
            x,
            enrich,
            linewidth=1.05,
            color=ENRICHMENT_COLOR,
            label="Enrichment",
            zorder=3,
        )
        ax.axhline(
            1.0,
            color=REFERENCE_COLOR,
            linewidth=0.65,
            linestyle="--",
            zorder=0,
        )
        if self.threshold is not None:
            ax.axhline(
                self.threshold,
                color=REFERENCE_COLOR,
                linewidth=0.75,
                linestyle=":",
                label="Peak threshold",
                zorder=1,
            )

        if cds_positions.size > 0:
            cds_start = int(cds_positions.min())
            cds_end = int(cds_positions.max())
            ax.axvline(
                cds_start - 0.5,
                color=REFERENCE_COLOR,
                linewidth=0.55,
                linestyle=":",
                zorder=0,
            )
            ax.axvline(
                cds_end + 0.5,
                color=REFERENCE_COLOR,
                linewidth=0.55,
                linestyle=":",
                zorder=0,
            )

        self._mark_peak_maxima(ax, table, regions)

        ax.set_xlabel("Position relative to CDS start (codons)", fontsize=self.font_size)
        ax.set_ylabel("IP / control enrichment", fontsize=self.font_size)
        ax.set_xlim(float(x.min()) - 0.5, float(x.max()) + 0.5)
        ax.set_ylim(bottom=0.0)
        if self.y_max is not None:
            ax.set_ylim(0.0, self.y_max)
        else:
            current_top = ax.get_ylim()[1]
            data_top = float(np.nanmax(finite_enrich))
            threshold_top = self.threshold if self.threshold is not None else 0.0
            desired_top = max(data_top, threshold_top, 1.0) * 1.08
            ax.set_ylim(0.0, max(current_top, desired_top))

        ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.tick_params(
            axis="both",
            labelsize=max(self.font_size - 1.0, 1.0),
            length=3.0,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                handles,
                labels,
                frameon=False,
                fontsize=max(self.font_size - 1.5, 6.0),
                loc="upper right",
                ncol=1,
                handlelength=1.8,
                borderaxespad=0.2,
            )

        gene_name = self._gene_name(transcript)
        if gene_name == "-" or gene_name == transcript:
            title = transcript
        else:
            title = "{gene} | {transcript}".format(
                gene=gene_name,
                transcript=transcript,
            )
        ax.set_title(title, fontsize=self.font_size + 0.5, pad=5.0)
        fig.tight_layout(pad=0.6)

        output_stem = self.figure_dir / self._safe_file_name(transcript)
        self._save_figure(fig, output_stem)
        plt.close(fig)

    def draw_targets(self) -> None:
        """Draw all resolved target transcripts and no others."""
        self.figure_dir.mkdir(parents=True, exist_ok=True)
        drawn = 0
        for transcript, table in self.profile_table.groupby("name", sort=False):
            transcript = str(transcript)
            self._draw_transcript(transcript, table)
            drawn += 1

        print("Figures generated: {num:,}".format(num=drawn), flush=True)
        print("Figure directory: {path}".format(path=self.figure_dir), flush=True)