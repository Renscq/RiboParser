#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.006
# Function: Detect transcript-level ribosomal frameshift candidates from RPF density data.
# Input: RPF density file in JSONL or TXT format and an optional transcript filter.
# Output: Frame periodicity, change-point statistics, candidate tables, summaries, and figures.

"""Core functions for transcript-level ribosomal frameshift detection.

The legacy implementation classified a transcript from its whole-CDS frame
proportions. That approach detects frame dominance, but it does not demonstrate
a frame transition. This implementation scans internal CDS change points and
compares the upstream and downstream three-frame count distributions.

For each candidate split, a 2 x 3 likelihood-ratio G-test is calculated. The
within-transcript scan p-value is Bonferroni-corrected, followed by BH-FDR across
transcripts. A transcript is reported as a frameshift candidate only when the
statistical test and the biological frame-transition thresholds are both met.
"""

from __future__ import annotations

import json
import math
import os
import re
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

from . import RPFs

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
FRAME_NAMES = ["frame0", "frame1", "frame2"]
PROGRESS_EVERY = 500

FRAME_COLORS = {
    0: "#56B4E9",
    1: "#009E73",
    2: "#E69F00",
}

CLASSIFICATION_COLORS = OrderedDict([
    ("plus1_shift", "#0072B2"), 
    ("minus1_shift", "#D55E00"),
    ("fuzzy", "#CC79A7"), 
    ("no_shift", "#999999"), 
    ("insufficient", "#E6E6E6"),
])

SCATTER_COLORS = {
    "other": "#B0B0B0",
    "candidate": "#0072B2",
}

FRAME_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "frame_intensity", 
    ["#F7F7F7", "#BDD7E8", "#5B8DB8", "#1E3F5F"]
)


class Shift(object):
    """Detect sample-specific ribosomal frameshift candidates.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Shift``.
    """

    def __init__(self, args):
        # Input and output.
        self.rpf = args.rpf
        self.output = args.output
        self.transcript = args.transcript

        # Density and transcript filtering.
        self.site = args.site
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts
        self.thread = max(1, int(args.thread))

        # Change-point scan.
        self.period = float(args.period) / 100.0
        self.pre_period = float(args.pre_period) / 100.0
        self.min_shift = float(args.min_shift)
        self.min_segment = int(args.min_segment)
        self.min_segment_reads = float(args.min_segment_reads)
        self.scan_step = int(args.scan_step)
        self.alpha = float(args.alpha)

        # Outlier filtering.
        self.remove_outlier = args.remove_outlier
        self.outlier_iqr = float(args.outlier_iqr)
        self.outlier_window = int(args.outlier_window)
        self.outlier_local_fold = float(args.outlier_local_fold)

        # Plotting.
        self.smooth_window = int(args.smooth_window)
        self.plot_bin_size = int(args.plot_bin_size)
        self.gene_plot_mode = args.gene_plot_mode
        self.gene_fig = args.gene_fig
        self.max_gene_figures = int(args.max_gene_figures)

        # Imported data.
        self.rpf_data = None
        self.raw_rpf = None
        self.sample_name: list[str] = []
        self.sample_num = 0
        self.file_format = None

        # Results.
        self.sample_high_genes = OrderedDict()
        self.periodicity = pd.DataFrame()
        self.shift_results = pd.DataFrame()
        self.shift_candidates = pd.DataFrame()
        self.shift_summary = pd.DataFrame()
        self.outliers = pd.DataFrame()
        self.output_files = OrderedDict()

        # Precomputed transcript row boundaries for fast repeated sample scans.
        self.transcript_names: np.ndarray | None = None
        self.transcript_starts: np.ndarray | None = None
        self.transcript_ends: np.ndarray | None = None
        self.meta_table: pd.DataFrame | None = None

    # ------------------------------------------------------------------
    # Import and validation
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import frame-resolved RPF density data."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf,
            sample_name=None,
            gene=self.transcript,
            tis=self.tis,
            tts=self.tts,
        )
        self.raw_rpf = self.rpf_data.raw_rpf.to_pandas()
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.file_format = self.rpf_data.file_format

        required = list(BASE_COLUMNS)
        for sample in self.sample_name:
            required.extend([f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"])
        missing = [column for column in required if column not in self.raw_rpf.columns]
        if missing:
            raise ValueError(
                "RPF density table is missing required column(s): " + ", ".join(missing)
            )

        self.raw_rpf = self.raw_rpf.loc[
            self.raw_rpf["region"].eq("cds"), required
        ].copy()
        if self.raw_rpf.empty:
            raise ValueError("No CDS positions remain after input filtering.")

        frame_columns = [column for column in required if re.search(r"_f[012]$", column)]
        self.raw_rpf[frame_columns] = (
            self.raw_rpf[frame_columns]
            .apply(pd.to_numeric, errors="coerce")
            .fillna(0.0)
            .astype(float)
        )

        # E/P/A site adjustment shifts every frame-resolved density vector by codon.
        shift_num = RPFs.set_codon_shift(self.site)
        if shift_num != 0:
            self.raw_rpf.loc[:, frame_columns] = (
                self.raw_rpf.groupby("name", sort=False)[frame_columns]
                .shift(shift_num)
                .fillna(0.0)
                .to_numpy()
            )

        # Sort once and reuse transcript row boundaries for every sample.
        self.raw_rpf.sort_values(
            ["name", "now_nt", "from_tis"],
            kind="mergesort",
            inplace=True,
            ignore_index=True,
        )
        names = self.raw_rpf["name"].astype(str).to_numpy()
        starts = np.r_[0, np.flatnonzero(names[1:] != names[:-1]) + 1]
        ends = np.r_[starts[1:], len(names)]
        self.transcript_names = names[starts]
        self.transcript_starts = starts.astype(np.int64, copy=False)
        self.transcript_ends = ends.astype(np.int64, copy=False)
        self.meta_table = self.raw_rpf.loc[:, BASE_COLUMNS]

        print(
            "Imported frame-resolved RPF density: samples={samples}, format={fmt}, "
            "CDS rows={rows:,}.".format(
                samples=self.sample_num,
                fmt=self.file_format,
                rows=len(self.raw_rpf),
            ),
            flush=True,
        )

    # ------------------------------------------------------------------
    # Statistical helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _local_background(values: np.ndarray, window: int) -> np.ndarray:
        """Calculate a symmetric local mean excluding the focal position."""
        values = np.asarray(values, dtype=float)
        n = values.size
        if n == 0:
            return np.asarray([], dtype=float)

        finite = np.isfinite(values)
        safe = np.where(finite, values, 0.0)
        prefix_sum = np.concatenate(([0.0], np.cumsum(safe)))
        prefix_n = np.concatenate(([0.0], np.cumsum(finite.astype(float))))
        idx = np.arange(n)
        left = np.maximum(0, idx - window)
        right = np.minimum(n, idx + window + 1)
        local_sum = prefix_sum[right] - prefix_sum[left] - safe
        local_n = prefix_n[right] - prefix_n[left] - finite.astype(float)
        result = np.full(n, np.nan, dtype=float)
        np.divide(local_sum, local_n, out=result, where=local_n > 0)
        return result

    def _detect_outliers(
        self,
        sample: str,
        table: pd.DataFrame,
        frame_matrix: np.ndarray,
    ) -> tuple[np.ndarray, list[dict]]:
        """Detect isolated extreme total-density pileups for one transcript."""
        if not self.remove_outlier or frame_matrix.size == 0:
            return np.zeros(frame_matrix.shape[0], dtype=bool), []

        total = np.nansum(frame_matrix, axis=1)
        positive = total[np.isfinite(total) & (total > 0)]
        if positive.size < 8:
            return np.zeros(total.size, dtype=bool), []

        logged = np.log1p(positive)
        q1, q3 = np.quantile(logged, [0.25, 0.75])
        cutoff = np.expm1(q3 + self.outlier_iqr * (q3 - q1))
        candidate = np.isfinite(total) & (total > cutoff)
        if not candidate.any():
            return np.zeros(total.size, dtype=bool), []

        local = self._local_background(total, self.outlier_window)
        local_fold = (total + 1.0) / (local + 1.0)
        mask = candidate & np.isfinite(local_fold) & (
            local_fold >= self.outlier_local_fold
        )

        records = []
        for position in np.flatnonzero(mask):
            records.append(
                {
                    "Sample": sample,
                    "name": str(table.iloc[position]["name"]),
                    "now_nt": table.iloc[position]["now_nt"],
                    "from_tis": table.iloc[position]["from_tis"],
                    "from_tts": table.iloc[position]["from_tts"],
                    "codon": table.iloc[position]["codon"],
                    "TotalDensity": float(total[position]),
                    "GlobalCutoff": float(cutoff),
                    "LocalBackground": float(local[position]),
                    "LocalFold": float(local_fold[position]),
                }
            )
        return mask, records

    @staticmethod
    def _g_test_2x3(upstream: np.ndarray, downstream: np.ndarray) -> np.ndarray:
        """Calculate 2 x 3 likelihood-ratio G-test p-values.

        The asymptotic degrees of freedom equal 2, whose chi-square survival
        function has the closed form ``exp(-G / 2)``. SciPy is therefore not
        required.
        """
        upstream = np.asarray(upstream, dtype=float)
        downstream = np.asarray(downstream, dtype=float)
        observed = np.stack([upstream, downstream], axis=1)
        row_sum = observed.sum(axis=2, keepdims=True)
        col_sum = observed.sum(axis=1, keepdims=True)
        total = observed.sum(axis=(1, 2), keepdims=True)

        expected = np.divide(
            row_sum * col_sum,
            total,
            out=np.zeros_like(observed),
            where=total > 0,
        )
        valid = (observed > 0) & (expected > 0)
        terms = np.zeros_like(observed)
        terms[valid] = observed[valid] * np.log(observed[valid] / expected[valid])
        statistic = 2.0 * terms.sum(axis=(1, 2))
        pvalue = np.exp(-0.5 * np.maximum(statistic, 0.0))
        invalid = (row_sum[:, 0, 0] <= 0) | (row_sum[:, 1, 0] <= 0)
        pvalue[invalid] = 1.0
        return pvalue

    @staticmethod
    def _bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
        """Adjust p-values by the Benjamini-Hochberg procedure."""
        values = np.asarray(list(pvalues), dtype=float)
        result = np.full(values.shape, np.nan, dtype=float)
        finite_idx = np.flatnonzero(np.isfinite(values))
        if finite_idx.size == 0:
            return result
        finite = values[finite_idx]
        order = np.argsort(finite)
        ranked = finite[order]
        n = ranked.size
        adjusted = ranked * n / np.arange(1, n + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0.0, 1.0)
        restored = np.empty_like(adjusted)
        restored[order] = adjusted
        result[finite_idx] = restored
        return result

    # ------------------------------------------------------------------
    # Frameshift scan
    # ------------------------------------------------------------------

    def _scan_transcript(
        self,
        sample: str,
        transcript: str,
        table: pd.DataFrame,
        matrix: np.ndarray,
    ) -> tuple[dict, list[dict]]:
        """Scan one presorted sample-transcript matrix for a frame change point."""
        matrix = np.asarray(matrix, dtype=float).copy()

        outlier_mask, outlier_records = self._detect_outliers(sample, table, matrix)
        if outlier_mask.any():
            matrix[outlier_mask, :] = np.nan

        safe_matrix = np.where(np.isfinite(matrix), matrix, 0.0)
        total_frames = safe_matrix.sum(axis=0)
        total_count = float(total_frames.sum())
        n_positions = int(matrix.shape[0])
        observed_positions = int(np.isfinite(matrix).any(axis=1).sum())

        base = {
            "Sample": sample,
            "name": transcript,
            "PositionCount": n_positions,
            "ObservedPositionCount": observed_positions,
            "TotalCount": total_count,
            "Frame0Count": float(total_frames[0]),
            "Frame1Count": float(total_frames[1]),
            "Frame2Count": float(total_frames[2]),
            "Frame0Ratio": float(total_frames[0] / total_count) if total_count > 0 else np.nan,
            "Frame1Ratio": float(total_frames[1] / total_count) if total_count > 0 else np.nan,
            "Frame2Ratio": float(total_frames[2] / total_count) if total_count > 0 else np.nan,
            "BestSplitIndex": np.nan,
            "ShiftNowNt": np.nan,
            "ShiftFromTIS": np.nan,
            "ShiftCodon": None,
            "ShiftType": "none",
            "UpstreamCount": np.nan,
            "DownstreamCount": np.nan,
            "UpstreamFrame0Ratio": np.nan,
            "DownstreamTargetRatio": np.nan,
            "TargetUpstreamRatio": np.nan,
            "ShiftEffect": np.nan,
            "RawPValue": np.nan,
            "ScanAdjustedPValue": np.nan,
            "TestCount": 0,
            "CandidateRule": False,
        }

        if total_count < self.rpf_num or n_positions < 2 * self.min_segment + 1:
            base["Status"] = "insufficient"
            return base, outlier_records

        cumulative = np.cumsum(safe_matrix, axis=0)
        split_idx = np.arange(
            self.min_segment,
            n_positions - self.min_segment + 1,
            self.scan_step,
            dtype=int,
        )
        if split_idx.size == 0:
            base["Status"] = "insufficient"
            return base, outlier_records

        upstream = cumulative[split_idx - 1]
        downstream = cumulative[-1] - upstream
        upstream_total = upstream.sum(axis=1)
        downstream_total = downstream.sum(axis=1)

        upstream_ratio = np.divide(
            upstream,
            upstream_total[:, None],
            out=np.zeros_like(upstream),
            where=upstream_total[:, None] > 0,
        )
        downstream_ratio = np.divide(
            downstream,
            downstream_total[:, None],
            out=np.zeros_like(downstream),
            where=downstream_total[:, None] > 0,
        )

        target_frame = np.where(downstream_ratio[:, 1] >= downstream_ratio[:, 2], 1, 2)
        row_idx = np.arange(split_idx.size)
        target_down = downstream_ratio[row_idx, target_frame]
        target_up = upstream_ratio[row_idx, target_frame]
        effect = target_down - target_up
        raw_p = self._g_test_2x3(upstream, downstream)
        adjusted_p = np.minimum(raw_p * split_idx.size, 1.0)

        rules = (
            (upstream_total >= self.min_segment_reads)
            & (downstream_total >= self.min_segment_reads)
            & (upstream_ratio[:, 0] >= self.pre_period)
            & (target_down >= self.period)
            & (effect >= self.min_shift)
        )

        # Select the strongest biologically valid split; if none are valid,
        # retain the minimum adjusted-p split for diagnostic output.
        valid_idx = np.flatnonzero(rules)
        if valid_idx.size:
            candidate_order = np.lexsort((-effect[valid_idx], adjusted_p[valid_idx]))
            best = valid_idx[candidate_order[0]]
        else:
            best = int(np.nanargmin(adjusted_p))

        split = int(split_idx[best])
        target = int(target_frame[best])
        shift_type = "+1" if target == 1 else "-1"
        shift_row = table.iloc[min(split, n_positions - 1)]

        base.update(
            {
                "BestSplitIndex": split,
                "ShiftNowNt": shift_row["now_nt"],
                "ShiftFromTIS": shift_row["from_tis"],
                "ShiftCodon": shift_row["codon"],
                "ShiftType": shift_type,
                "UpstreamCount": float(upstream_total[best]),
                "DownstreamCount": float(downstream_total[best]),
                "UpstreamFrame0Ratio": float(upstream_ratio[best, 0]),
                "DownstreamTargetRatio": float(target_down[best]),
                "TargetUpstreamRatio": float(target_up[best]),
                "ShiftEffect": float(effect[best]),
                "RawPValue": float(raw_p[best]),
                "ScanAdjustedPValue": float(adjusted_p[best]),
                "TestCount": int(split_idx.size),
                "CandidateRule": bool(rules[best]),
                "Status": "tested",
            }
        )
        return base, outlier_records

    def _process_sample(self, sample: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Process one sample using precomputed transcript row boundaries."""
        if (
            self.transcript_names is None
            or self.transcript_starts is None
            or self.transcript_ends is None
            or self.meta_table is None
        ):
            raise RuntimeError("Transcript row boundaries were not initialized.")

        frame_columns = [f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"]
        sample_matrix = self.raw_rpf.loc[:, frame_columns].to_numpy(
            dtype=float, copy=False
        )
        row_total = sample_matrix.sum(axis=1)
        gene_counts = np.add.reduceat(row_total, self.transcript_starts)
        keep = gene_counts >= self.rpf_num
        keep_idx = np.flatnonzero(keep)
        self.sample_high_genes[sample] = pd.Index(self.transcript_names[keep_idx])

        results: list[dict] = []
        outliers: list[dict] = []
        for number, gene_idx in enumerate(keep_idx, start=1):
            start = int(self.transcript_starts[gene_idx])
            end = int(self.transcript_ends[gene_idx])
            transcript = str(self.transcript_names[gene_idx])
            result, outlier_records = self._scan_transcript(
                sample=sample,
                transcript=transcript,
                table=self.meta_table.iloc[start:end],
                matrix=sample_matrix[start:end, :],
            )
            results.append(result)
            outliers.extend(outlier_records)
            if number % PROGRESS_EVERY == 0:
                print(
                    f"Sample {sample}: scanned {number:,}/{len(keep_idx):,} transcripts.",
                    flush=True,
                )

        print(
            f"Sample {sample}: scanned {len(results):,} transcripts.",
            flush=True,
        )
        return pd.DataFrame(results), pd.DataFrame(outliers)

    def scan_frame_shift(self) -> None:
        """Scan all samples and transcripts for frame-transition candidates."""
        workers = min(self.thread, max(1, self.sample_num))
        sample_results = {}
        sample_outliers = {}

        if workers == 1:
            for sample in self.sample_name:
                result, outlier = self._process_sample(sample)
                sample_results[sample] = result
                sample_outliers[sample] = outlier
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(self._process_sample, sample): sample
                    for sample in self.sample_name
                }
                for future in as_completed(futures):
                    sample = futures[future]
                    result, outlier = future.result()
                    sample_results[sample] = result
                    sample_outliers[sample] = outlier

        self.shift_results = pd.concat(
            [sample_results[sample] for sample in self.sample_name],
            ignore_index=True,
        )
        outlier_tables = [
            sample_outliers[sample]
            for sample in self.sample_name
            if not sample_outliers[sample].empty
        ]
        self.outliers = (
            pd.concat(outlier_tables, ignore_index=True)
            if outlier_tables
            else pd.DataFrame()
        )

        # FDR is applied independently within each sample after the within-gene
        # Bonferroni correction for scanned split positions.
        self.shift_results["FDR"] = np.nan
        for sample, index in self.shift_results.groupby("Sample", sort=False).groups.items():
            idx = np.asarray(list(index), dtype=int)
            self.shift_results.loc[idx, "FDR"] = self._bh_adjust(
                self.shift_results.loc[idx, "ScanAdjustedPValue"]
            )

        significant = (
            self.shift_results["CandidateRule"].fillna(False)
            & self.shift_results["FDR"].le(self.alpha)
        )
        self.shift_results["Classification"] = "no_shift"
        self.shift_results.loc[
            self.shift_results["Status"].eq("insufficient"), "Classification"
        ] = "insufficient"
        self.shift_results.loc[
            self.shift_results["Status"].eq("tested")
            & ~self.shift_results["CandidateRule"].fillna(False),
            "Classification",
        ] = "fuzzy"
        self.shift_results.loc[significant, "Classification"] = self.shift_results.loc[
            significant, "ShiftType"
        ].map({"+1": "plus1_shift", "-1": "minus1_shift"})

        self.shift_candidates = self.shift_results.loc[significant, :].copy()
        summary_order = ["plus1_shift", "minus1_shift", "fuzzy", "no_shift", "insufficient"]
        self.shift_summary = (
            self.shift_results.groupby(["Sample", "Classification"], sort=False)
            .size()
            .rename("Count")
            .reset_index()
        )
        full_index = pd.MultiIndex.from_product(
            [self.sample_name, summary_order], names=["Sample", "Classification"]
        )
        self.shift_summary = (
            self.shift_summary.set_index(["Sample", "Classification"])
            .reindex(full_index, fill_value=0)
            .reset_index()
        )

        periodicity_columns = [
            "Sample", "name", "PositionCount", "TotalCount",
            "Frame0Count", "Frame1Count", "Frame2Count",
            "Frame0Ratio", "Frame1Ratio", "Frame2Ratio",
        ]
        self.periodicity = self.shift_results.loc[:, periodicity_columns].copy()

    # Backward-compatible method names.
    calc_3nt_period = scan_frame_shift

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def output_results(self) -> None:
        """Write frameshift result tables and summary metadata."""
        periodicity_file = self.output + "_gene_periodicity.txt"
        result_file = self.output + "_frame_shift_scan.txt"
        candidate_file = self.output + "_frame_shift_candidates.txt"
        summary_file = self.output + "_frame_shift_count.txt"

        self.periodicity.to_csv(periodicity_file, sep="\t", index=False)
        self.shift_results.to_csv(result_file, sep="\t", index=False)
        self.shift_candidates.to_csv(candidate_file, sep="\t", index=False)
        self.shift_summary.to_csv(summary_file, sep="\t", index=False)

        self.output_files["periodicity"] = periodicity_file
        self.output_files["scan"] = result_file
        self.output_files["candidates"] = candidate_file
        self.output_files["count"] = summary_file

        if self.remove_outlier:
            outlier_file = self.output + "_frame_shift.outliers.txt"
            self.outliers.to_csv(outlier_file, sep="\t", index=False)
            self.output_files["outliers"] = outlier_file

        summary_json = self.output + "_frame_shift.summary.json"
        payload = {
            "version": "0.2.8-dev.006",
            "input": os.path.abspath(self.rpf),
            "input_format": self.file_format,
            "samples": self.sample_name,
            "parameters": {
                "site": self.site,
                "min_rpf": self.rpf_num,
                "tis": self.tis,
                "tts": self.tts,
                "pre_period_percent": self.pre_period * 100.0,
                "post_period_percent": self.period * 100.0,
                "min_shift": self.min_shift,
                "min_segment": self.min_segment,
                "min_segment_reads": self.min_segment_reads,
                "scan_step": self.scan_step,
                "alpha": self.alpha,
                "remove_outlier": self.remove_outlier,
                "thread": self.thread,
                "smooth_window": self.smooth_window,
                "plot_bin_size": self.plot_bin_size,
                "gene_plot_mode": self.gene_plot_mode,
            },
            "candidate_count": int(len(self.shift_candidates)),
            "sample_candidate_count": {
                sample: int((self.shift_candidates["Sample"] == sample).sum())
                for sample in self.sample_name
            },
            "outputs": dict(self.output_files),
        }
        with open(summary_json, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
        self.output_files["summary"] = summary_json

    # Historical output methods.
    output_meta = lambda self: None
    filter_frame_shift = lambda self: None
    output_frame_shift = output_results

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_filename(value: str) -> str:
        """Convert an identifier to a filesystem-safe filename."""
        return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._") or "transcript"

    def draw_frame_shift_count(self) -> None:
        """Draw sample-level candidate classification counts."""
        if self.shift_summary.empty:
            return
        pivot = self.shift_summary.pivot(
            index="Sample", columns="Classification", values="Count"
        ).fillna(0)
        order = [
            column for column in CLASSIFICATION_COLORS
            if column in pivot.columns
        ]
        pivot = pivot.reindex(index=self.sample_name, columns=order, fill_value=0)

        width = max(7.5, 0.65 * len(self.sample_name) + 4.5)
        fig, ax = plt.subplots(figsize=(width, 5.8), dpi=300)
        x = np.arange(len(pivot), dtype=float)
        bottom = np.zeros(len(pivot), dtype=float)
        for column in order:
            values = pivot[column].to_numpy(dtype=float)
            bars = ax.bar(
                x, values, bottom=bottom, width=0.72,
                label=column.replace("_", " "),
                color=CLASSIFICATION_COLORS.get(column),
                edgecolor="white", linewidth=0.6,
            )
            bottom += values

        ax.set_xlim(-0.5, len(x) - 0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(pivot.index.to_list(), rotation=45, ha="right")
        ax.set_ylabel("Transcript count")
        ax.set_title("Frame-transition classification")
        ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.35)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        fig.tight_layout()

        pdf = self.output + "_frame_shift_count_plot.pdf"
        png = self.output + "_frame_shift_count_plot.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["count_plot_pdf"] = pdf
        self.output_files["count_plot_png"] = png

    def draw_candidate_scatter(self) -> None:
        """Draw effect-size versus statistical-significance scatter plots."""
        tested = self.shift_results.loc[
            self.shift_results["Status"].eq("tested")
            & self.shift_results["FDR"].notna()
            & self.shift_results["ShiftEffect"].notna(),
            :,
        ].copy()
        if tested.empty:
            return
        tested["NegLog10FDR"] = -np.log10(
            np.clip(tested["FDR"].to_numpy(dtype=float), 1e-300, 1.0)
        )

        samples = self.sample_name
        ncol = min(3, max(1, len(samples)))
        nrow = int(math.ceil(len(samples) / ncol))
        fig, axes = plt.subplots(
            nrow, ncol,
            figsize=(5.0 * ncol, 4.0 * nrow),
            dpi=300,
            squeeze=False,
        )
        for ax, sample in zip(axes.ravel(), samples):
            data = tested.loc[tested["Sample"].eq(sample), :]
            candidate = data["Classification"].isin(["plus1_shift", "minus1_shift"])
            ax.scatter(
                data.loc[~candidate, "ShiftEffect"],
                data.loc[~candidate, "NegLog10FDR"],
                s=12,
                alpha=0.4,
                color=SCATTER_COLORS["other"],
                label="Other",
            )
            ax.scatter(
                data.loc[candidate, "ShiftEffect"],
                data.loc[candidate, "NegLog10FDR"],
                s=22,
                alpha=0.85,
                color=SCATTER_COLORS["candidate"],
                label="Candidate",
            )
            ax.axvline(self.min_shift, linewidth=0.8, linestyle="--")
            ax.axhline(-math.log10(self.alpha), linewidth=0.8, linestyle="--")
            ax.set_title(sample)
            ax.set_xlabel("Frame-shift effect")
            ax.set_ylabel("-log10(FDR)")
            ax.spines[["top", "right"]].set_visible(False)
        for ax in axes.ravel()[len(samples):]:
            ax.set_visible(False)
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, frameon=False, loc="upper center", ncol=2)
        fig.suptitle("Frameshift candidate scan", y=1.01)
        fig.tight_layout()

        pdf = self.output + "_frame_shift_scatter.pdf"
        png = self.output + "_frame_shift_scatter.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["scatter_pdf"] = pdf
        self.output_files["scatter_png"] = png

    @staticmethod
    def _rolling_mean(values: np.ndarray, window: int) -> np.ndarray:
        """Calculate a centered rolling mean for plotting."""
        return (
            pd.Series(values, dtype=float)
            .rolling(window=max(1, window), center=True, min_periods=1)
            .mean()
            .to_numpy()
        )

    @staticmethod
    def _bin_frame_proportions(
        proportion: np.ndarray,
        x: np.ndarray,
        bin_size: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Aggregate frame proportions into consecutive codon bins for plotting."""
        proportion = np.asarray(proportion, dtype=float)
        x = np.asarray(x, dtype=float)
        if proportion.size == 0:
            return np.asarray([], dtype=float), np.empty((0, 3), dtype=float)
        bin_size = max(1, int(bin_size))
        n = proportion.shape[0]
        centers = []
        rows = []
        for start in range(0, n, bin_size):
            stop = min(start + bin_size, n)
            chunk = proportion[start:stop, :]
            if chunk.size == 0:
                continue
            rows.append(np.nanmean(chunk, axis=0))
            centers.append(float(np.nanmean(x[start:stop])))
        return np.asarray(centers, dtype=float), np.vstack(rows)

    @staticmethod
    def _draw_stacked_bar(ax, centers: np.ndarray, matrix: np.ndarray) -> None:
        """Draw a 100% stacked bar plot of frame proportions."""
        if matrix.size == 0:
            return
        width = 0.9 if len(centers) <= 1 else max(0.5, 0.85 * np.min(np.diff(centers)))
        bottom = np.zeros(len(centers), dtype=float)
        for frame in range(3):
            values = matrix[:, frame]
            ax.bar(
                centers, values, bottom=bottom, width=width,
                color=FRAME_COLORS[frame], edgecolor="white", linewidth=0.2,
                label=f"Frame {frame}",
            )
            bottom += values
        ax.set_ylim(0, 1)
        ax.set_ylabel("Frame proportion")
        ax.set_xlabel("CDS codon position")
        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.3)
        ax.spines[["top", "right"]].set_visible(False)

    @staticmethod
    def _draw_heatmap(ax, matrix: np.ndarray, extent: tuple[float, float, float, float]) -> None:
        """Draw a frame-proportion heatmap."""
        if matrix.size == 0:
            return
        im = ax.imshow(
            matrix.T, aspect="auto", interpolation="nearest",
            cmap=FRAME_HEATMAP_CMAP, vmin=0, vmax=1, origin="lower", extent=extent,
        )
        ax.set_yticks([0.5, 1.5, 2.5])
        ax.set_yticklabels(["Frame 0", "Frame 1", "Frame 2"])
        ax.set_xlabel("CDS codon position")
        ax.set_ylabel("")
        ax.spines[["top", "right"]].set_visible(False)
        return im

    def draw_candidate_genes(self) -> None:
        """Draw three vertically aligned candidate-level frame profiles."""
        if self.gene_fig == "none" or self.shift_candidates.empty:
            return

        plot_dir = Path(self.output + "_frame_shift_geneplots")
        plot_dir.mkdir(parents=True, exist_ok=True)
        # Remove stale figures with the same output prefix to avoid mixing layouts.
        for pattern in ("*_frame_shift.png", "*_frame_shift.pdf"):
            for old_file in plot_dir.glob(pattern):
                old_file.unlink(missing_ok=True)
        candidates = self.shift_candidates.sort_values(
            ["FDR", "ShiftEffect"], ascending=[True, False]
        )
        if self.max_gene_figures > 0:
            candidates = candidates.head(self.max_gene_figures)

        for number, row in enumerate(candidates.itertuples(index=False), start=1):
            sample = row.Sample
            transcript = row.name
            columns = BASE_COLUMNS + [
                f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"
            ]
            table = self.raw_rpf.loc[
                self.raw_rpf["name"].eq(transcript), columns
            ].sort_values(["now_nt", "from_tis"], kind="mergesort")
            if table.empty:
                continue

            matrix = table[[f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"]].to_numpy(
                dtype=float, copy=True
            )
            x = np.arange(1, len(table) + 1, dtype=float)
            smooth = np.column_stack(
                [self._rolling_mean(matrix[:, frame], self.smooth_window) for frame in range(3)]
            )
            denom = smooth.sum(axis=1)
            proportion = np.divide(
                smooth,
                denom[:, None],
                out=np.zeros_like(smooth),
                where=denom[:, None] > 0,
            )
            bar_x, bar_prop = self._bin_frame_proportions(
                proportion=proportion,
                x=x,
                bin_size=self.plot_bin_size,
            )

            fig, axes = plt.subplots(
                3, 1,
                figsize=(12.0, 9.6),
                dpi=300,
                sharex=True,
                gridspec_kw={"height_ratios": [1.15, 1.0, 0.72], "hspace": 0.16},
            )
            ax_density, ax_bar, ax_heat = axes
            shift_x = int(row.BestSplitIndex)
            ax_density.text(-0.075, 1.06, "A", transform=ax_density.transAxes,
                            fontsize=12, fontweight="bold", va="top")
            ax_bar.text(-0.075, 1.06, "B", transform=ax_bar.transAxes,
                        fontsize=12, fontweight="bold", va="top")
            ax_heat.text(-0.075, 1.06, "C", transform=ax_heat.transAxes,
                         fontsize=12, fontweight="bold", va="top")

            # Panel 1: smoothed frame density.
            for frame in range(3):
                ax_density.plot(
                    x, smooth[:, frame],
                    linewidth=1.5,
                    color=FRAME_COLORS[frame],
                    label=f"Frame {frame}",
                )
            ax_density.set_ylabel("Smoothed RPF density")
            ax_density.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.25)
            ax_density.spines[["top", "right"]].set_visible(False)
            ax_density.legend(
                frameon=False,
                ncol=3,
                loc="upper center",
                bbox_to_anchor=(0.5, 1.13),
                columnspacing=1.5,
                handlelength=2.2,
            )
            ax_density.set_title(
                f"{transcript} | {sample} | {row.ShiftType} shift | "
                f"FDR={row.FDR:.3g} | effect={row.ShiftEffect:.3f}",
                pad=26,
            )

            # Panel 2: 100% stacked frame-proportion bar plot.
            self._draw_stacked_bar(ax_bar, bar_x, bar_prop)
            ax_bar.set_xlabel("")

            # Panel 3: frame-proportion heatmap.
            extent = (x.min() - 0.5, x.max() + 0.5, 0.0, 3.0)
            im = self._draw_heatmap(ax_heat, proportion, extent)
            ax_heat.set_xlabel("CDS codon position")
            ax_heat.set_title("Frame-proportion heatmap", loc="left", fontsize=10, pad=6)
            cbar = fig.colorbar(im, ax=ax_heat, orientation="horizontal", fraction=0.10, pad=0.28, aspect=35)
            cbar.set_label("Frame proportion")

            for ax in axes:
                ax.axvline(shift_x, linestyle="--", linewidth=1.0, color="#555555")
                ax.set_xlim(1, len(table))

            fig.subplots_adjust(top=0.90, bottom=0.13, left=0.10, right=0.98, hspace=0.22)

            stem = (
                self._safe_filename(str(transcript))
                + "__"
                + self._safe_filename(str(sample))
                + "_frame_shift"
            )
            if self.gene_fig in {"png", "both"}:
                fig.savefig(plot_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
            if self.gene_fig in {"pdf", "both"}:
                fig.savefig(plot_dir / f"{stem}.pdf", bbox_inches="tight")
            plt.close(fig)

            if number % 100 == 0 or number == len(candidates):
                print(
                    f"Generated frameshift figures for {number}/{len(candidates)} candidates.",
                    flush=True,
                )

        self.output_files["gene_plot_directory"] = str(plot_dir)

    def draw_all(self) -> None:
        """Draw all summary and candidate-level figures."""
        self.draw_frame_shift_count()
        self.draw_candidate_scatter()
        self.draw_candidate_genes()
