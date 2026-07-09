#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.004
# Function: Provide RiboParser core functions for Metaplot analysis.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Metaplot summary table, outlier table, sample-level plots, and heatmap figures.

"""Metaplot analysis utilities for RiboParser.

This module calculates RPF metaplots around start and stop codons from the
updated RiboParser density format. The RPF import layer is delegated to
``utils.ribo.RPFs.RPFData`` so both compact JSONL density files and legacy TXT
codon-level density files are supported through one workflow.
"""

from __future__ import annotations

import json
import os
import re
from collections import OrderedDict
from typing import Iterable

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
META_COLORS = {
    "TIS": "#3C5488",
    "TTS": "#E64B35",
}
HEATMAP_CMAP = "Blues"
RPM_SCALE = 1_000_000.0
PLOT_TRANSFORM_LABELS = {
    "none": "",
    "sqrt": "sqrt",
    "log": "log1p",
    "log1p": "log1p",
    "log2": "log2(x + 1)",
    "log10": "log10(x + 1)",
}


class Metaplot(object):
    """Calculate and plot RPF metaplots around TIS and TTS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Metaplot``.
    """

    def __init__(self, args):
        # Input and output.
        self.transcript = args.transcript
        self.rpf = args.rpf
        self.output = args.output

        # Window and filtering parameters.
        self.rpf_num = args.min
        self.utr5 = args.utr5
        self.cds = args.cds
        self.utr3 = args.utr3
        self.norm = args.normal
        self.normal_scale = RPM_SCALE

        # Plot transformation and heatmap scaling parameters.
        self.plot_transform = args.plot_transform
        self.scale = args.scale

        # Outlier filtering parameters.
        self.remove_outlier = args.remove_outlier
        self.outlier_iqr = args.outlier_iqr
        self.outlier_window = args.outlier_window
        self.outlier_local_fold = args.outlier_local_fold

        # Imported RPF data.
        self.rpf_data = None
        self.raw_rpf = None
        self.sample_name = []
        self.sample_num = 0
        self.total_rpf_num = None
        self.file_format = None
        self.window_gene = None

        # Metaplot results.
        self.meta = None
        self.outliers = pd.DataFrame()
        self.sample_summary = pd.DataFrame()
        self.sample_dict = OrderedDict()
        self.output_files = OrderedDict()

    # ------------------------------------------------------------------
    # Import and validation
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import RPF density data through the updated ``RPFs`` module.

        Returns
        -------
        None
            The function updates internal state in place.
        """
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf,
            sample_name=None,
            gene=self.transcript,
            tis=None,
            tts=None,
        )

        self.raw_rpf = self.rpf_data.raw_rpf.to_pandas()
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf_num = self.rpf_data.total_rpf_num
        self.file_format = self.rpf_data.file_format

        self._validate_rpf_table()
        self.window_gene = self._get_complete_window_genes()

        print(
            "Imported RPF density for {sample_num} sample(s), format={fmt}, rows={rows:,}.".format(
                sample_num=self.sample_num,
                fmt=self.file_format,
                rows=len(self.raw_rpf),
            ),
            flush=True,
        )
        print(
            "Transcripts with complete TIS/TTS metaplot windows: {count:,}.".format(
                count=len(self.window_gene),
            ),
            flush=True,
        )

    def _validate_rpf_table(self) -> None:
        """Validate required columns in the imported RPF table."""
        if self.raw_rpf is None or self.raw_rpf.empty:
            raise ValueError("RPF density table is empty after import.")

        missing_base = [column for column in BASE_COLUMNS if column not in self.raw_rpf.columns]
        if missing_base:
            raise ValueError(
                "RPF density table is missing required base column(s): {columns}".format(
                    columns=", ".join(missing_base)
                )
            )

        for sample in self.sample_name:
            missing_frame = [column for column in self._sample_frame_columns(sample) if column not in self.raw_rpf.columns]
            if missing_frame:
                raise ValueError(
                    "RPF density table is missing frame column(s) for sample {sample}: {columns}".format(
                        sample=sample,
                        columns=", ".join(missing_frame),
                    )
                )

    def _window_gene_ranges(self) -> pd.DataFrame:
        """Return per-transcript coordinate ranges used for window checks."""
        if self.raw_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        gene_range = self.raw_rpf.groupby("name", sort=False).agg(
            from_tis_min=("from_tis", "min"),
            from_tis_max=("from_tis", "max"),
            from_tts_min=("from_tts", "min"),
            from_tts_max=("from_tts", "max"),
        )

        gene_range["available_utr5"] = (-gene_range["from_tis_min"]).clip(lower=0).astype(int)
        gene_range["available_start_cds"] = (gene_range["from_tis_max"] + 1).clip(lower=0).astype(int)
        gene_range["available_stop_cds"] = (1 - gene_range["from_tts_min"]).clip(lower=0).astype(int)
        gene_range["available_cds"] = gene_range[["available_start_cds", "available_stop_cds"]].min(axis=1).astype(int)
        gene_range["available_utr3"] = gene_range["from_tts_max"].clip(lower=0).astype(int)
        return gene_range

    @staticmethod
    def _genes_for_window(gene_range: pd.DataFrame, utr5: int, cds: int, utr3: int) -> pd.Index:
        """Return transcripts that fully cover a candidate metaplot window."""
        start_min = -int(utr5)
        start_max = int(cds) - 1
        stop_min = -int(cds) + 1
        stop_max = int(utr3)

        keep = gene_range.loc[
            (gene_range["from_tis_min"] <= start_min)
            & (gene_range["from_tis_max"] >= start_max)
            & (gene_range["from_tts_min"] <= stop_min)
            & (gene_range["from_tts_max"] >= stop_max)
        ]
        return keep.index

    def _candidate_window_lengths(self, gene_range: pd.DataFrame) -> list[tuple[int, int, int]]:
        """Generate data-supported fallback window candidates.

        The first candidate is the user-requested window. Additional candidates
        are generated by proportional shrinking and by transcript-capacity
        quantiles. This avoids hard failure when a species or annotation lacks
        long UTR windows, while keeping the selected fallback as close as
        possible to the requested profile.
        """
        requested = (int(self.utr5), int(self.cds), int(self.utr3))
        candidates: list[tuple[int, int, int]] = [requested]

        # Proportional shrinkage preserves the requested window shape.
        for scale in np.linspace(0.95, 0.05, 19):
            cand = (
                max(0, int(np.floor(requested[0] * scale))),
                max(1, int(np.floor(requested[1] * scale))),
                max(0, int(np.floor(requested[2] * scale))),
            )
            candidates.append(cand)

        # Capacity quantiles adapt to annotations with no or short UTRs.
        capacity_columns = ["available_utr5", "available_cds", "available_utr3"]
        for quantile in (0.90, 0.75, 0.50, 0.25, 0.10, 0.05):
            caps = gene_range.loc[:, capacity_columns].quantile(quantile).fillna(0).astype(int)
            cand = (
                min(requested[0], max(0, int(caps["available_utr5"]))),
                min(requested[1], max(1, int(caps["available_cds"]))),
                min(requested[2], max(0, int(caps["available_utr3"]))),
            )
            candidates.append(cand)

        # Individual transcript capacities guarantee that at least one complete
        # candidate can be found when any transcript has valid TIS/TTS coverage.
        if not gene_range.empty:
            capped = gene_range.assign(
                cand_utr5=gene_range["available_utr5"].clip(upper=requested[0]).astype(int),
                cand_cds=gene_range["available_cds"].clip(upper=requested[1]).astype(int),
                cand_utr3=gene_range["available_utr3"].clip(upper=requested[2]).astype(int),
            )
            capped = capped.loc[capped["cand_cds"] >= 1, ["cand_utr5", "cand_cds", "cand_utr3"]]
            capped = capped.drop_duplicates().sort_values(
                ["cand_cds", "cand_utr5", "cand_utr3"], ascending=False
            )
            for _, row in capped.head(200).iterrows():
                candidates.append((int(row["cand_utr5"]), int(row["cand_cds"]), int(row["cand_utr3"])))

        candidates.append((0, 1, 0))

        # Deduplicate while preserving order.
        unique_candidates = []
        seen = set()
        for cand in candidates:
            cand = (max(0, int(cand[0])), max(1, int(cand[1])), max(0, int(cand[2])))
            if cand not in seen:
                unique_candidates.append(cand)
                seen.add(cand)
        return unique_candidates

    def _get_complete_window_genes(self) -> pd.Index:
        """Return genes covering the requested or auto-adjusted metaplot windows."""
        gene_range = self._window_gene_ranges()
        requested = (int(self.utr5), int(self.cds), int(self.utr3))

        keep = self._genes_for_window(gene_range, *requested)
        if len(keep) > 0:
            return keep

        best_candidate = None
        best_keep: pd.Index | None = None
        best_score = -1.0

        for candidate in self._candidate_window_lengths(gene_range):
            candidate_keep = self._genes_for_window(gene_range, *candidate)
            if len(candidate_keep) == 0:
                continue

            # Prefer windows that retain enough genes while preserving profile
            # length. The logarithm prevents very small windows with many genes
            # from always dominating a biologically useful larger fallback.
            window_size = candidate[0] + candidate[1] * 2 + candidate[2]
            score = window_size * float(np.log1p(len(candidate_keep)))
            if score > best_score:
                best_candidate = candidate
                best_keep = candidate_keep
                best_score = score

        if best_candidate is None or best_keep is None or len(best_keep) == 0:
            raise ValueError(
                "No transcript has valid TIS/TTS coordinates for metaplot analysis after automatic window adjustment."
            )

        self.utr5, self.cds, self.utr3 = best_candidate
        print(
            "Warning: requested metaplot window --utr5 {old_utr5}, --cds {old_cds}, --utr3 {old_utr3} "
            "has no complete transcript. Auto-adjusted to --utr5 {new_utr5}, --cds {new_cds}, --utr3 {new_utr3}.".format(
                old_utr5=requested[0],
                old_cds=requested[1],
                old_utr3=requested[2],
                new_utr5=self.utr5,
                new_cds=self.cds,
                new_utr3=self.utr3,
            ),
            flush=True,
        )
        return best_keep

    # ------------------------------------------------------------------
    # Calculation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _sample_frame_columns(sample: str) -> list[str]:
        """Return frame density columns for one sample."""
        return [f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"]

    @staticmethod
    def _safe_name(name: str) -> str:
        """Return a safe sample name for output file names."""
        return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_") or "sample"

    def _sample_total(self, sample: str) -> float:
        """Return total raw RPF count for one sample."""
        if self.total_rpf_num is not None and sample in self.total_rpf_num.index:
            total = float(self.total_rpf_num.loc[sample])
        else:
            frame_columns = self._sample_frame_columns(sample)
            total = float(self.raw_rpf.loc[:, frame_columns].sum().sum())
        return total

    def _get_high_expression_genes(self, sample: str) -> pd.Index:
        """Return sample-specific high-expression transcripts."""
        frame_columns = self._sample_frame_columns(sample)
        cds_rpf = self.raw_rpf.loc[self.raw_rpf["region"] == "cds", ["name"] + frame_columns].copy()
        cds_rpf["Count"] = cds_rpf.loc[:, frame_columns].sum(axis=1)
        gene_count = cds_rpf.groupby("name", sort=False)["Count"].sum()
        high_gene = gene_count.loc[gene_count >= self.rpf_num].index

        high_gene = high_gene.intersection(self.window_gene)
        return high_gene

    def _long_profile(
        self,
        sample: str,
        meta_name: str,
        coord_column: str,
        codon_positions: Iterable[int],
        genes: pd.Index,
    ) -> pd.DataFrame:
        """Build a nucleotide-level long table for one sample and one meta region."""
        frame_columns = self._sample_frame_columns(sample)
        codon_positions = list(codon_positions)
        sample_total = self._sample_total(sample)

        region_df = self.raw_rpf.loc[
            self.raw_rpf["name"].isin(genes) & self.raw_rpf[coord_column].isin(codon_positions),
            ["name", coord_column] + frame_columns,
        ].copy()

        if region_df.empty:
            return pd.DataFrame(
                columns=[
                    "Sample",
                    "Meta",
                    "Transcript",
                    "Codon",
                    "Frame",
                    "Nucleotide",
                    "RawDensity",
                    "Density",
                ]
            )

        frames = []
        for frame_index, frame_column in enumerate(frame_columns):
            now_frame = region_df.loc[:, ["name", coord_column, frame_column]].copy()
            now_frame.rename(
                columns={
                    "name": "Transcript",
                    coord_column: "Codon",
                    frame_column: "RawDensity",
                },
                inplace=True,
            )
            now_frame["Sample"] = sample
            now_frame["Meta"] = meta_name
            now_frame["Frame"] = frame_index
            now_frame["Nucleotide"] = now_frame["Codon"].astype(int) * 3 + frame_index
            if self.norm and sample_total > 0:
                now_frame["Density"] = now_frame["RawDensity"].astype(float) * self.normal_scale / sample_total
            else:
                now_frame["Density"] = now_frame["RawDensity"].astype(float)
            frames.append(now_frame)

        long_df = pd.concat(frames, axis=0, ignore_index=True)
        long_df = long_df[
            ["Sample", "Meta", "Transcript", "Codon", "Frame", "Nucleotide", "RawDensity", "Density"]
        ]
        return long_df

    @staticmethod
    def _robust_upper_cutoff(values: pd.Series, multiplier: float) -> float:
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

        # The larger cutoff is more conservative and avoids deleting profiles
        # that are globally high but not true isolated pileups.
        return float(np.expm1(max(cutoffs)))

    def _candidate_local_background(
        self,
        long_df: pd.DataFrame,
        candidate_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add local background metrics only for global outlier candidates.

        Local-background estimation is intentionally restricted to points that
        already pass the global robust cutoff. This avoids the expensive
        transcript-by-frame rolling calculation over every codon position and
        makes ``--remove-outlier`` scale with the number of candidate pileups
        rather than with the full metaplot table size.
        """
        if candidate_df.empty:
            candidate_df = candidate_df.copy()
            candidate_df["LocalBackground"] = []
            candidate_df["LocalFold"] = []
            return candidate_df

        window = max(1, int(self.outlier_window))

        density_lookup = (
            long_df.loc[:, ["Transcript", "Frame", "Codon", "RawDensity"]]
            .set_index(["Transcript", "Frame", "Codon"])["RawDensity"]
            .astype(float)
        )

        local_backgrounds = []
        local_folds = []

        for row in candidate_df.itertuples(index=False):
            transcript = getattr(row, "Transcript")
            frame = int(getattr(row, "Frame"))
            codon = int(getattr(row, "Codon"))
            raw_density = float(getattr(row, "RawDensity"))

            neighbor_values = []
            for offset in range(1, window + 1):
                left_key = (transcript, frame, codon - offset)
                right_key = (transcript, frame, codon + offset)

                left_value = density_lookup.get(left_key, np.nan)
                right_value = density_lookup.get(right_key, np.nan)

                if pd.notna(left_value):
                    neighbor_values.append(float(left_value))
                if pd.notna(right_value):
                    neighbor_values.append(float(right_value))

            if neighbor_values:
                local_background = float(np.median(neighbor_values))
            else:
                local_background = 0.0

            local_backgrounds.append(local_background)
            local_folds.append((raw_density + 1.0) / (local_background + 1.0))

        candidate_df = candidate_df.copy()
        candidate_df["LocalBackground"] = local_backgrounds
        candidate_df["LocalFold"] = local_folds
        return candidate_df

    def _remove_outlier_points(self, long_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Remove dynamic, locally isolated extreme gene-position pileups.

        The optimized implementation first identifies global extreme candidates
        on the log1p robust scale, then computes local background only for those
        candidates. This preserves the conservative isolated-pileup logic from
        v0.2.8-dev.003 but avoids transcript-by-frame rolling calculations over
        the full long table.
        """
        outlier_columns = list(long_df.columns) + [
            "LocalBackground",
            "LocalFold",
            "OutlierCutoff",
            "OutlierMethod",
            "OutlierReason",
        ]

        if not self.remove_outlier or long_df.empty:
            return long_df, pd.DataFrame(columns=outlier_columns)

        raw_density = long_df["RawDensity"].astype(float)
        cutoff = self._robust_upper_cutoff(raw_density, float(self.outlier_iqr))

        if not np.isfinite(cutoff):
            return long_df, pd.DataFrame(columns=outlier_columns)

        candidate_mask = raw_density > cutoff
        if not bool(candidate_mask.any()):
            return long_df, pd.DataFrame(columns=outlier_columns)

        candidate_df = long_df.loc[candidate_mask, :].copy()
        candidate_df = self._candidate_local_background(long_df, candidate_df)

        local_fold = (
            candidate_df["LocalFold"]
            .astype(float)
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
        )
        keep_candidate_mask = local_fold >= float(self.outlier_local_fold)

        outlier_df = candidate_df.loc[keep_candidate_mask, :].copy()
        if outlier_df.empty:
            return long_df, pd.DataFrame(columns=outlier_columns)

        outlier_df["OutlierCutoff"] = cutoff
        outlier_df["OutlierMethod"] = "fast_log1p_robust_cutoff_and_candidate_local_peak"
        outlier_df["OutlierReason"] = (
            "RawDensity > dynamic robust cutoff and LocalFold >= "
            + str(float(self.outlier_local_fold))
        )

        outlier_index = outlier_df.index
        clean_df = long_df.drop(index=outlier_index).copy()
        outlier_df = outlier_df.reindex(columns=outlier_columns)
        return clean_df, outlier_df

    @staticmethod
    def _complete_position_frame_grid(
        sample: str,
        meta_name: str,
        codon_positions: Iterable[int],
    ) -> pd.DataFrame:
        """Create a complete codon-frame coordinate grid."""
        records = []
        for codon in codon_positions:
            for frame in range(3):
                records.append(
                    {
                        "Sample": sample,
                        "Meta": meta_name,
                        "Codon": int(codon),
                        "Frame": frame,
                        "Nucleotide": int(codon) * 3 + frame,
                    }
                )
        return pd.DataFrame.from_records(records)

    def _aggregate_profile(
        self,
        clean_df: pd.DataFrame,
        sample: str,
        meta_name: str,
        codon_positions: Iterable[int],
        gene_count: int,
        outlier_count: int,
    ) -> pd.DataFrame:
        """Aggregate gene-position density to a metaplot profile."""
        grid = self._complete_position_frame_grid(sample, meta_name, codon_positions)

        if clean_df.empty:
            aggregated = grid.copy()
            aggregated["Density"] = 0.0
            aggregated["RawDensity"] = 0.0
            aggregated["ObservedGeneCount"] = 0
        else:
            aggregated = (
                clean_df.groupby(["Sample", "Meta", "Codon", "Frame", "Nucleotide"], sort=False)
                .agg(
                    Density=("Density", "mean"),
                    RawDensity=("RawDensity", "mean"),
                    ObservedGeneCount=("Transcript", "nunique"),
                )
                .reset_index()
            )
            aggregated = grid.merge(
                aggregated,
                on=["Sample", "Meta", "Codon", "Frame", "Nucleotide"],
                how="left",
            )
            aggregated["Density"] = aggregated["Density"].fillna(0.0)
            aggregated["RawDensity"] = aggregated["RawDensity"].fillna(0.0)
            aggregated["ObservedGeneCount"] = aggregated["ObservedGeneCount"].fillna(0).astype(int)

        aggregated["TotalGeneCount"] = int(gene_count)
        aggregated["OutlierCount"] = int(outlier_count)
        aggregated.sort_values(["Nucleotide", "Frame"], inplace=True)
        return aggregated[
            [
                "Sample",
                "Meta",
                "Nucleotide",
                "Codon",
                "Frame",
                "Density",
                "RawDensity",
                "ObservedGeneCount",
                "TotalGeneCount",
                "OutlierCount",
            ]
        ]

    # ------------------------------------------------------------------
    # Public calculation workflow
    # ------------------------------------------------------------------

    def calc_metaplot(self) -> None:
        """Calculate metaplot profiles for all samples."""
        if self.raw_rpf is None:
            raise ValueError("RPF density data has not been imported yet.")

        tis_codons = list(range(-int(self.utr5), int(self.cds)))
        tts_codons = list(range(-int(self.cds) + 1, int(self.utr3) + 1))

        meta_tables = []
        outlier_tables = []
        summary_records = []

        for sample in self.sample_name:
            high_gene = self._get_high_expression_genes(sample)
            gene_count = len(high_gene)

            if gene_count == 0:
                print(
                    "Warning: sample {sample} has no transcript after expression and window filtering.".format(
                        sample=sample
                    ),
                    flush=True,
                )

            sample_meta_tables = []
            sample_outlier_count = 0

            for meta_name, coord_column, codon_positions in (
                ("TIS", "from_tis", tis_codons),
                ("TTS", "from_tts", tts_codons),
            ):
                long_df = self._long_profile(
                    sample=sample,
                    meta_name=meta_name,
                    coord_column=coord_column,
                    codon_positions=codon_positions,
                    genes=high_gene,
                )
                clean_df, outlier_df = self._remove_outlier_points(long_df)
                outlier_count = len(outlier_df)
                sample_outlier_count += outlier_count

                aggregated = self._aggregate_profile(
                    clean_df=clean_df,
                    sample=sample,
                    meta_name=meta_name,
                    codon_positions=codon_positions,
                    gene_count=gene_count,
                    outlier_count=outlier_count,
                )
                meta_tables.append(aggregated)
                sample_meta_tables.append(aggregated)

                if not outlier_df.empty:
                    outlier_tables.append(outlier_df)

                summary_records.append(
                    {
                        "Sample": sample,
                        "Meta": meta_name,
                        "HighGeneCount": int(gene_count),
                        "OutlierCount": int(outlier_count),
                    }
                )

            self.sample_dict[sample] = {
                "gene_count": int(gene_count),
                "outlier_count": int(sample_outlier_count),
                "meta": pd.concat(sample_meta_tables, axis=0, ignore_index=True),
            }
            print(
                "Sample {sample}: high genes={genes:,}, removed outlier points={outliers:,}.".format(
                    sample=sample,
                    genes=gene_count,
                    outliers=sample_outlier_count,
                ),
                flush=True,
            )

        if not meta_tables:
            raise ValueError("No metaplot result was generated.")

        self.meta = pd.concat(meta_tables, axis=0, ignore_index=True)
        self.sample_summary = pd.DataFrame.from_records(summary_records)
        if outlier_tables:
            self.outliers = pd.concat(outlier_tables, axis=0, ignore_index=True)
        else:
            self.outliers = pd.DataFrame(
                columns=[
                    "Sample",
                    "Meta",
                    "Transcript",
                    "Codon",
                    "Frame",
                    "Nucleotide",
                    "RawDensity",
                    "Density",
                    "LocalBackground",
                    "LocalFold",
                    "OutlierCutoff",
                    "OutlierMethod",
                    "OutlierReason",
                ]
            )

    # ------------------------------------------------------------------
    # Output tables and summary
    # ------------------------------------------------------------------

    def output_meta(self) -> None:
        """Write metaplot summary tables."""
        if self.meta is None:
            raise ValueError("Metaplot has not been calculated yet.")

        out_txt = self.output + "_tis_tts_metaplot.txt"
        self.meta.to_csv(out_txt, sep="\t", index=False)
        self.output_files["metaplot_table"] = out_txt

        summary_txt = self.output + "_metaplot.sample_summary.txt"
        self.sample_summary.to_csv(summary_txt, sep="\t", index=False)
        self.output_files["sample_summary"] = summary_txt

        if self.remove_outlier:
            outlier_txt = self.output + "_metaplot.outliers.txt"
            self.outliers.to_csv(outlier_txt, sep="\t", index=False)
            self.output_files["outlier_table"] = outlier_txt
            print(
                "Outlier table written: {file} (removed points={count:,}).".format(
                    file=outlier_txt,
                    count=len(self.outliers),
                ),
                flush=True,
            )

    def write_summary(self) -> None:
        """Write a machine-readable metaplot summary JSON."""
        summary_json = self.output + "_metaplot.summary.json"
        self.output_files["summary_json"] = summary_json
        summary = OrderedDict(
            [
                ("tool", "rpf_Metaplot"),
                ("version", "0.2.8-dev.004"),
                ("input_rpf", os.path.abspath(self.rpf)),
                ("input_format", self.file_format),
                ("transcript_filter", os.path.abspath(self.transcript) if self.transcript else None),
                ("sample_count", self.sample_num),
                ("samples", self.sample_name),
                ("parameters", OrderedDict(
                    [
                        ("min", self.rpf_num),
                        ("utr5", self.utr5),
                        ("cds", self.cds),
                        ("utr3", self.utr3),
                        ("normal", self.norm),
                        ("rpm_scale", RPM_SCALE),
                        ("plot_transform", self.plot_transform),
                        ("heatmap_scale", self.scale),
                        ("remove_outlier", self.remove_outlier),
                        ("outlier_iqr", self.outlier_iqr),
                        ("outlier_window", self.outlier_window),
                        ("outlier_local_fold", self.outlier_local_fold),
                    ]
                )),
                ("complete_window_gene_count", int(len(self.window_gene)) if self.window_gene is not None else 0),
                ("sample_summary", self.sample_summary.to_dict(orient="records")),
                ("output_files", self.output_files),
            ]
        )
        with open(summary_json, "w", encoding="utf-8") as out:
            json.dump(summary, out, ensure_ascii=False, indent=2)
            out.write("\n")

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _set_common_axis_style(ax) -> None:
        """Apply common axis style."""
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linewidth=0.4, alpha=0.25)
        ax.axvline(0, color="#222222", linewidth=0.8, linestyle="--", alpha=0.75)

    @staticmethod
    def _nice_tick_step(span: float, max_ticks: int = 8) -> int:
        """Return a regular codon-aware tick interval for a nucleotide axis."""
        span = abs(float(span))
        if span <= 0:
            return 1

        raw_step = span / max(1, max_ticks - 1)

        # Metaplot coordinates are nucleotide offsets but biologically read as
        # codon-scale windows. Prefer 15/30/60 nt style ticks instead of odd
        # values produced by linear index sampling.
        preferred_steps = [3, 6, 15, 30, 60, 90, 150, 300, 600, 900, 1500]
        for step in preferred_steps:
            if step >= raw_step:
                return step

        exponent = int(np.ceil(np.log10(raw_step)))
        return int(10 ** exponent)

    @classmethod
    def _common_tick_step(cls, *position_lists: list[int], max_ticks: int = 8) -> int:
        """Return one shared tick interval for TIS and TTS panels."""
        spans = []
        for positions in position_lists:
            if positions:
                spans.append(max(positions) - min(positions))
        if not spans:
            return 1
        return cls._nice_tick_step(max(spans), max_ticks=max_ticks)

    @staticmethod
    def _regular_ticks(x_min: int, x_max: int, step: int, max_ticks: int = 9) -> list[int]:
        """Return regular integer ticks within one axis range and preserve zero."""
        if x_min > x_max:
            x_min, x_max = x_max, x_min

        step = max(1, int(step))
        lower = int(np.floor(x_min / step) * step)
        upper = int(np.ceil(x_max / step) * step)

        # Include a rounded boundary tick when it is very close to the data
        # boundary. This gives axes such as -30, -15, 0, 15, 30 instead of
        # ending at 29 when the last nucleotide is one base before 30.
        edge_tolerance = step * 0.25
        start = lower if (x_min - lower) <= edge_tolerance else int(np.ceil(x_min / step) * step)
        end = upper if (upper - x_max) <= edge_tolerance else int(np.floor(x_max / step) * step)
        ticks = list(range(start, end + step, step)) if start <= end else []

        if x_min <= 0 <= x_max and 0 not in ticks:
            ticks.append(0)

        ticks = sorted(set(int(tick) for tick in ticks))

        while len(ticks) > max_ticks:
            nonzero_ticks = [tick for tick in ticks if tick != 0]
            if len(nonzero_ticks) <= 1:
                break
            # Keep every second tick while always preserving zero.
            reduced = nonzero_ticks[::2]
            if 0 in ticks:
                reduced.append(0)
            ticks = sorted(set(reduced))

        if not ticks:
            ticks = [x_min]
            if x_min <= 0 <= x_max:
                ticks.append(0)
            if x_max != x_min:
                ticks.append(x_max)
            ticks = sorted(set(ticks))

        return ticks

    @classmethod
    def _set_regular_xticks(cls, ax, positions: list[int], tick_step: int, max_ticks: int = 9) -> None:
        """Set horizontal regular x-axis ticks and keep the zero position."""
        if not positions:
            return
        x_min = int(min(positions))
        x_max = int(max(positions))
        tick_positions = cls._regular_ticks(x_min, x_max, tick_step, max_ticks=max_ticks)
        x_left = min(x_min, min(tick_positions)) if tick_positions else x_min
        x_right = max(x_max, max(tick_positions)) if tick_positions else x_max
        ax.set_xlim(x_left - 0.5, x_right + 0.5)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([str(value) for value in tick_positions], rotation=0, ha="center")

    def _transform_density_values(self, values) -> pd.Series:
        """Transform density values for plotting only."""
        transformed = pd.Series(values, dtype="float64").replace([np.inf, -np.inf], np.nan).fillna(0.0)
        transformed = transformed.clip(lower=0.0)

        if self.plot_transform == "none":
            return transformed
        if self.plot_transform == "sqrt":
            return np.sqrt(transformed)
        if self.plot_transform in {"log", "log1p"}:
            return np.log1p(transformed)
        if self.plot_transform == "log2":
            return np.log2(transformed + 1.0)
        if self.plot_transform == "log10":
            return np.log10(transformed + 1.0)

        raise ValueError("Unsupported plot transform: {transform}".format(transform=self.plot_transform))

    def _density_axis_label(self) -> str:
        """Return y-axis or colorbar label for plotted density values."""
        base = "Mean RPF density"
        if self.norm:
            base = "Mean RPF density (RPM)"

        label = PLOT_TRANSFORM_LABELS.get(self.plot_transform, "")
        if label:
            return "{base}, {label} transformed".format(base=base, label=label)
        return base

    def _plot_one_panel(
        self,
        ax,
        data: pd.DataFrame,
        meta_name: str,
        mode: str,
        gene_count: int,
        tick_step: int,
    ) -> None:
        """Draw one TIS or TTS metaplot panel."""
        data = data.sort_values(["Nucleotide", "Frame"]).copy()
        if data.empty:
            ax.set_axis_off()
            return

        data["PlotDensity"] = self._transform_density_values(data["Density"]).to_numpy(dtype=float)

        if mode == "bar":
            for frame in FRAME_ORDER:
                frame_df = data.loc[data["Frame"].astype(str) == frame, :]
                ax.bar(
                    frame_df["Nucleotide"],
                    frame_df["PlotDensity"],
                    width=0.85,
                    color=FRAME_COLORS[frame],
                    edgecolor="none",
                    label=f"Frame {frame}",
                )
        elif mode == "line":
            ax.plot(
                data["Nucleotide"],
                data["PlotDensity"],
                color=META_COLORS.get(meta_name, "#3C5488"),
                linewidth=1.4,
            )
        else:
            raise ValueError("mode must be 'bar' or 'line'.")

        title_region = "start codon" if meta_name == "TIS" else "stop codon"
        ax.set_title("{meta} region ({genes:,} genes)".format(meta=title_region, genes=gene_count), fontsize=11)
        ax.set_xlabel("Position from {region} (nt)".format(region=title_region))
        ax.set_ylabel(self._density_axis_label())
        self._set_common_axis_style(ax)
        self._set_regular_xticks(
            ax,
            data["Nucleotide"].drop_duplicates().astype(int).tolist(),
            tick_step=tick_step,
            max_ticks=9,
        )

    def _sync_pair_yaxis(self, axes, *tables: pd.DataFrame) -> None:
        """Use the same y-axis range for paired TIS and TTS panels."""
        values = []
        for table in tables:
            if table is None or table.empty or "Density" not in table.columns:
                continue
            values.append(self._transform_density_values(table["Density"]).to_numpy(dtype=float))

        if not values:
            y_top = 1.0
        else:
            merged = np.concatenate(values)
            merged = merged[np.isfinite(merged)]
            y_max = float(np.max(merged)) if merged.size else 0.0
            y_top = y_max * 1.08 if y_max > 0 else 1.0

        for ax in axes:
            ax.set_ylim(0, y_top)

    def draw_metaplot(self, mode: str) -> None:
        """Draw per-sample metaplot figures."""
        if self.meta is None:
            raise ValueError("Metaplot has not been calculated yet.")

        if mode not in {"bar", "line"}:
            raise ValueError("mode must be 'bar' or 'line'.")

        for sample, message in self.sample_dict.items():
            sample_meta = message["meta"]
            gene_count = message["gene_count"]
            safe_sample = self._safe_name(sample)

            out_pdf = f"{self.output}_{safe_sample}_meta_{mode}_plot.pdf"
            out_png = f"{self.output}_{safe_sample}_meta_{mode}_plot.png"

            fig, axes = plt.subplots(
                nrows=1,
                ncols=2,
                figsize=(14, 4.8),
                gridspec_kw={"wspace": 0.22},
            )

            tis_meta = sample_meta.loc[sample_meta["Meta"] == "TIS", :]
            tts_meta = sample_meta.loc[sample_meta["Meta"] == "TTS", :]
            tick_step = self._common_tick_step(
                tis_meta["Nucleotide"].drop_duplicates().astype(int).tolist(),
                tts_meta["Nucleotide"].drop_duplicates().astype(int).tolist(),
                max_ticks=9,
            )
            self._plot_one_panel(axes[0], tis_meta, "TIS", mode, gene_count, tick_step)
            self._plot_one_panel(axes[1], tts_meta, "TTS", mode, gene_count, tick_step)
            self._sync_pair_yaxis(axes, tis_meta, tts_meta)

            if mode == "bar":
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
                    bbox_to_anchor=(0.5, -0.04),
                    fontsize=9,
                )
                fig.subplots_adjust(bottom=0.20, left=0.07, right=0.98, top=0.88)
            else:
                fig.subplots_adjust(bottom=0.14, left=0.07, right=0.98, top=0.88)

            fig.savefig(out_pdf, bbox_inches="tight")
            fig.savefig(out_png, dpi=300, bbox_inches="tight")
            plt.close(fig)

            self.output_files[f"{sample}_meta_{mode}_pdf"] = out_pdf
            self.output_files[f"{sample}_meta_{mode}_png"] = out_png

    def _matrix_for_heatmap(self, meta_name: str) -> pd.DataFrame:
        """Return sample-by-position matrix for heatmap plotting."""
        data = self.meta.loc[self.meta["Meta"] == meta_name, ["Sample", "Nucleotide", "Density"]].copy()
        data["PlotDensity"] = self._transform_density_values(data["Density"]).to_numpy(dtype=float)
        matrix = data.pivot_table(
            index="Sample",
            columns="Nucleotide",
            values="PlotDensity",
            aggfunc="mean",
            fill_value=0,
        )
        matrix = matrix.reindex(index=self.sample_name)
        matrix = matrix.reindex(columns=sorted(matrix.columns))

        if self.scale == "row":
            matrix = self._row_minmax_scale(matrix)

        return matrix

    @staticmethod
    def _row_minmax_scale(matrix: pd.DataFrame) -> pd.DataFrame:
        """Scale each heatmap row to the [0, 1] range."""
        scaled = matrix.astype(float).copy()
        row_min = scaled.min(axis=1)
        row_max = scaled.max(axis=1)
        row_range = (row_max - row_min).replace(0, np.nan)
        scaled = scaled.sub(row_min, axis=0).div(row_range, axis=0).fillna(0.0)
        return scaled

    @classmethod
    def _set_heatmap_xticks(cls, ax, positions: list[int], tick_step: int, max_ticks: int = 9) -> None:
        """Set heatmap x-axis ticks with horizontal labels."""
        if not positions:
            return

        x_min = int(min(positions))
        x_max = int(max(positions))
        position_to_index = {int(position): idx for idx, position in enumerate(positions)}
        tick_values = cls._regular_ticks(x_min, x_max, tick_step, max_ticks=max_ticks)
        tick_values = [value for value in tick_values if value in position_to_index]

        if x_min <= 0 <= x_max and 0 in position_to_index and 0 not in tick_values:
            tick_values.append(0)
            tick_values = sorted(set(tick_values))

        ax.set_xticks([position_to_index[value] for value in tick_values])
        ax.set_xticklabels([str(value) for value in tick_values], rotation=0, ha="center", fontsize=8)

    def _plot_heatmap_panel(self, ax, matrix: pd.DataFrame, title: str, tick_step: int):
        """Draw one heatmap panel."""
        image = ax.imshow(
            matrix.to_numpy(dtype=float),
            aspect="auto",
            interpolation="nearest",
            cmap=HEATMAP_CMAP,
        )
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Relative position (nt)")

        if len(matrix.index) <= 60:
            ax.set_yticks(range(len(matrix.index)))
            ax.set_yticklabels(matrix.index.tolist(), fontsize=7)
        else:
            tick_indices = np.linspace(0, len(matrix.index) - 1, 20, dtype=int)
            ax.set_yticks(tick_indices)
            ax.set_yticklabels([matrix.index[idx] for idx in tick_indices], fontsize=6)

        positions = [int(value) for value in matrix.columns.tolist()]
        self._set_heatmap_xticks(ax, positions, tick_step=tick_step, max_ticks=9)
        return image

    def draw_heatmap(self) -> None:
        """Draw the combined TIS/TTS metaplot heatmap."""
        if self.meta is None:
            raise ValueError("Metaplot has not been calculated yet.")

        tis_matrix = self._matrix_for_heatmap("TIS")
        tts_matrix = self._matrix_for_heatmap("TTS")

        sample_count = max(1, len(self.sample_name))
        figure_height = min(max(4.8, sample_count * 0.22 + 2.0), 18.0)
        figure_width = 14.5

        fig, axes = plt.subplots(
            nrows=1,
            ncols=2,
            figsize=(figure_width, figure_height),
            gridspec_kw={"wspace": 0.18},
        )

        tick_step = self._common_tick_step(
            [int(value) for value in tis_matrix.columns.tolist()],
            [int(value) for value in tts_matrix.columns.tolist()],
            max_ticks=9,
        )
        image1 = self._plot_heatmap_panel(axes[0], tis_matrix, "TIS metaplot", tick_step=tick_step)
        image2 = self._plot_heatmap_panel(axes[1], tts_matrix, "TTS metaplot", tick_step=tick_step)
        axes[0].set_ylabel("Sample")
        axes[1].set_ylabel("")

        max_value = max(float(np.nanmax(tis_matrix.to_numpy())), float(np.nanmax(tts_matrix.to_numpy())))
        image1.set_clim(0, max_value if max_value > 0 else 1)
        image2.set_clim(0, max_value if max_value > 0 else 1)
        cbar = fig.colorbar(image2, ax=axes.ravel().tolist(), shrink=0.65, pad=0.02)
        heatmap_label = "Row-scaled " + self._density_axis_label() if self.scale == "row" else self._density_axis_label()
        cbar.set_label(heatmap_label)

        out_pdf = self.output + "_metaplot_heatmap.pdf"
        out_png = self.output + "_metaplot_heatmap.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        self.output_files["heatmap_pdf"] = out_pdf
        self.output_files["heatmap_png"] = out_png

    def run(self, mode: str) -> None:
        """Run the full Metaplot workflow."""
        self.import_rpf()
        self.calc_metaplot()
        self.output_meta()

        if mode in {"bar", "both"}:
            self.draw_metaplot("bar")
        if mode in {"line", "both"}:
            self.draw_metaplot("line")

        self.draw_heatmap()
        self.write_summary()
