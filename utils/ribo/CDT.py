#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.002
# Function: Calculate sample-specific codon decoding time from paired RPF and RNA density data.
# Input: RPF and RNA density files in JSONL or TXT format and an optional transcript filter.
# Output: Codon decoding time tables, outlier records, summary JSON, correlation, heatmap, and rank plots.

"""Core functions for RiboParser codon decoding time analysis.

RPF and RNA density files are imported through ``RPFs.RPFData``. Each RPF
sample is paired with one RNA sample, high-expression transcripts are selected
independently for every pair, and codon decoding time is calculated with
vectorized group operations.

For transcript ``g``, codon ``c`` and sample pair ``s``::

    absolute_cdt[c, s] = sum(RPF density at c) / valid RPF-supported positions

    normalized_cdt[c, s] =
        sum(RPF density at c / RNA RPKM of g) / valid RPF-supported positions

RNA RPKM uses the effective CDS interval that actually remains after TIS/TTS
filtering, so the count numerator and length denominator refer to the same
analysis window.
"""

from __future__ import annotations

import json
import os
import re
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import RPFs

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
STOP_CODONS = {"TAA", "TAG", "TGA"}
RPM_SCALE = 1_000_000.0
RPKM_SCALE = 1_000_000_000.0
HEATMAP_ANNOTATION_MAX_SAMPLES = 12
PLOT_CMAP = "RdBu_r"


class CodonDecodingTime(object):
    """Calculate sample-specific codon decoding time.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_CDT``.
    """

    def __init__(self, args):
        self.rpf_file = args.rpf
        self.rna_file = args.rna
        self.transcript = args.list
        self.output = args.output

        self.site = args.site
        self.frame = args.frame
        self.rpf_num = args.min
        self.rna_num = args.min_rna
        self.tis = args.tis
        self.tts = args.tts
        self.scale = args.scale
        self.plot_transform = args.plot_transform
        self.rankplot_ncol = max(1, int(getattr(args, "rankplot_ncol", 1)))
        self.thread = max(1, int(args.thread))
        self.output_all = args.all

        self.remove_outlier = args.remove_outlier
        self.outlier_iqr = args.outlier_iqr
        self.outlier_window = args.outlier_window
        self.outlier_local_fold = args.outlier_local_fold

        self.rpf_data = None
        self.rna_data = None
        self.rpf_table = None
        self.rna_table = None
        self.rpf_samples = []
        self.rna_samples = []
        self.sample_pairs = []
        self.rpf_format = None
        self.rna_format = None

        self.position_cdt = None
        self.gene_codon_cdt = None
        self.codon_cdt = None
        self.outliers = pd.DataFrame()
        self.sample_summary = pd.DataFrame()
        self.output_files = OrderedDict()

        _, codon_table = RPFs.codon_table()
        self.codon_annotation = codon_table.loc[codon_table["Abbr"] != "*", :].copy()

    # ------------------------------------------------------------------
    # Import and validation
    # ------------------------------------------------------------------

    def import_density(self) -> None:
        """Import paired RPF and RNA density files."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.rpf_file,
            sample_name=None,
            gene=self.transcript,
            tis=self.tis,
            tts=self.tts,
        )
        self.rna_data = RPFs.RPFData.from_file(
            rpf_file=self.rna_file,
            sample_name=None,
            gene=self.transcript,
            tis=self.tis,
            tts=self.tts,
        )

        self.rpf_samples = list(self.rpf_data.sample_name)
        self.rna_samples = list(self.rna_data.sample_name)
        self.rpf_format = self.rpf_data.file_format
        self.rna_format = self.rna_data.file_format
        self.sample_pairs = self._pair_samples(self.rpf_samples, self.rna_samples)

        self.rpf_table = self.rpf_data.get_frame(frame=self.frame)
        self.rpf_table = RPFs.shift_site(
            self.rpf_table,
            self.rpf_samples,
            RPFs.set_codon_shift(self.site),
        )
        # RNA-seq density has no ribosomal E/P/A-site interpretation. It is
        # always merged across the selected frame without an additional shift.
        self.rna_table = self.rna_data.get_frame(frame=self.frame)

        self._validate_density_table(self.rpf_table, self.rpf_samples, "RPF")
        self._validate_density_table(self.rna_table, self.rna_samples, "RNA")

        keep_rpf = BASE_COLUMNS + self.rpf_samples
        keep_rna = BASE_COLUMNS + self.rna_samples
        self.rpf_table = self.rpf_table.loc[
            (self.rpf_table["region"] == "cds")
            & (~self.rpf_table["codon"].isin(STOP_CODONS)),
            keep_rpf,
        ].copy()
        self.rna_table = self.rna_table.loc[
            (self.rna_table["region"] == "cds")
            & (~self.rna_table["codon"].isin(STOP_CODONS)),
            keep_rna,
        ].copy()

        for sample in self.rpf_samples:
            self.rpf_table[sample] = pd.to_numeric(self.rpf_table[sample], errors="coerce").fillna(0.0)
        for sample in self.rna_samples:
            self.rna_table[sample] = pd.to_numeric(self.rna_table[sample], errors="coerce").fillna(0.0)

        if self.rpf_table.empty or self.rna_table.empty:
            raise ValueError("No CDS codons remain after filtering RPF and RNA density tables.")

        print(
            "Imported RPF format={rpf_fmt}, RNA format={rna_fmt}, sample pairs={pairs}.".format(
                rpf_fmt=self.rpf_format,
                rna_fmt=self.rna_format,
                pairs=len(self.sample_pairs),
            ),
            flush=True,
        )
        for output_sample, rpf_sample, rna_sample in self.sample_pairs:
            print(
                f"Pair {output_sample}: RPF={rpf_sample}, RNA={rna_sample}.",
                flush=True,
            )

    @staticmethod
    def _validate_density_table(table: pd.DataFrame, samples: list[str], label: str) -> None:
        """Validate one imported codon-level table."""
        missing = [column for column in BASE_COLUMNS + samples if column not in table.columns]
        if missing:
            raise ValueError(f"{label} density table is missing column(s): {', '.join(missing)}")

    @staticmethod
    def _pair_samples(rpf_samples: list[str], rna_samples: list[str]) -> list[tuple[str, str, str]]:
        """Pair RPF and RNA samples by name or, when necessary, by order."""
        if not rpf_samples or not rna_samples:
            raise ValueError("RPF and RNA density files must both contain at least one sample.")

        if set(rpf_samples) == set(rna_samples):
            return [(sample, sample, sample) for sample in rpf_samples]

        if len(rpf_samples) != len(rna_samples):
            raise ValueError(
                "RPF and RNA sample names differ and sample counts are unequal. "
                "Use files with matching sample names or equal ordered samples."
            )

        print(
            "Warning: RPF and RNA sample names differ; samples are paired by column order.",
            flush=True,
        )
        return [
            (rpf_sample, rpf_sample, rna_sample)
            for rpf_sample, rna_sample in zip(rpf_samples, rna_samples)
        ]

    # Backward-compatible legacy methods.
    import_gene = lambda self: None
    merge_rpf_rna = lambda self: None
    normalize_density = lambda self: None

    # ------------------------------------------------------------------
    # Outlier helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _robust_upper_cutoff(values: pd.Series, multiplier: float) -> float:
        """Return a conservative robust upper cutoff on the log1p scale."""
        positive = values.loc[values > 0].astype(float)
        if positive.size < 4:
            return float("inf")
        log_values = np.log1p(positive.to_numpy(dtype=float))
        q1, q3 = np.percentile(log_values, [25, 75])
        iqr = float(q3 - q1)
        median = float(np.median(log_values))
        mad = float(np.median(np.abs(log_values - median))) * 1.4826
        cutoffs = []
        if iqr > 0:
            cutoffs.append(float(q3 + multiplier * iqr))
        if mad > 0:
            cutoffs.append(float(median + multiplier * mad))
        return float(np.expm1(max(cutoffs))) if cutoffs else float("inf")

    @staticmethod
    def _local_background(values: np.ndarray, window: int) -> np.ndarray:
        """Calculate a symmetric local mean excluding each focal codon."""
        values = np.asarray(values, dtype=float)
        n = values.size
        result = np.full(n, np.nan, dtype=float)
        if n == 0:
            return result
        finite = np.isfinite(values)
        clean = np.where(finite, values, 0.0)
        prefix_sum = np.concatenate(([0.0], np.cumsum(clean)))
        prefix_n = np.concatenate(([0], np.cumsum(finite.astype(int))))
        idx = np.arange(n)
        left = np.maximum(0, idx - int(window))
        right = np.minimum(n, idx + int(window) + 1)
        sums = prefix_sum[right] - prefix_sum[left] - np.where(finite, clean, 0.0)
        counts = prefix_n[right] - prefix_n[left] - finite.astype(int)
        np.divide(sums, counts, out=result, where=counts > 0)
        return result

    def _remove_sample_outliers(
        self,
        table: pd.DataFrame,
        sample: str,
        sample_label: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Remove isolated extreme RPF pileups for one sample."""
        work = table.copy()
        if not self.remove_outlier or work.empty:
            return work, pd.DataFrame()

        cutoff = self._robust_upper_cutoff(work[sample], float(self.outlier_iqr))
        if not np.isfinite(cutoff):
            return work, pd.DataFrame()

        candidate_idx = work.index[work[sample].astype(float) > cutoff]
        if len(candidate_idx) == 0:
            return work, pd.DataFrame()

        outlier_records = []
        for gene, candidate_gene in work.loc[candidate_idx].groupby("name", sort=False):
            gene_idx = work.index[work["name"] == gene]
            gene_values = work.loc[gene_idx, sample].to_numpy(dtype=float)
            background = self._local_background(gene_values, self.outlier_window)
            index_to_pos = {idx: pos for pos, idx in enumerate(gene_idx)}
            for idx in candidate_gene.index:
                pos = index_to_pos[idx]
                raw = float(work.at[idx, sample])
                local = float(background[pos]) if np.isfinite(background[pos]) else 0.0
                fold = (raw + 1.0) / (local + 1.0)
                if fold >= float(self.outlier_local_fold):
                    record = work.loc[idx, BASE_COLUMNS].to_dict()
                    record.update(
                        {
                            "Sample": sample_label,
                            "RawDensity": raw,
                            "LocalBackground": local,
                            "LocalFold": fold,
                            "OutlierCutoff": cutoff,
                        }
                    )
                    outlier_records.append(record)
                    work.at[idx, sample] = np.nan

        return work, pd.DataFrame.from_records(outlier_records)

    # ------------------------------------------------------------------
    # CDT calculation
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
        """Divide arrays while returning NaN for invalid denominators."""
        result = np.full(np.asarray(numerator).shape, np.nan, dtype=float)
        valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0)
        np.divide(numerator, denominator, out=result, where=valid)
        return result

    def _calculate_sample_pair(
        self,
        output_sample: str,
        rpf_sample: str,
        rna_sample: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        """Calculate position-level CDT components for one RPF/RNA pair."""
        rpf = self.rpf_table.loc[:, BASE_COLUMNS + [rpf_sample]].copy()
        rna = self.rna_table.loc[:, BASE_COLUMNS + [rna_sample]].copy()
        rpf, outliers = self._remove_sample_outliers(rpf, rpf_sample, output_sample)

        rpf_gene_count = rpf.groupby("name", sort=False)[rpf_sample].sum(min_count=1)
        rna_gene_count = rna.groupby("name", sort=False)[rna_sample].sum(min_count=1)
        rna_gene_codons = rna.groupby("name", sort=False).size().astype(float)

        high_rpf = set(rpf_gene_count.index[rpf_gene_count >= float(self.rpf_num)])
        high_rna = set(rna_gene_count.index[rna_gene_count >= float(self.rna_num)])
        overlap = pd.Index([gene for gene in rpf_gene_count.index if gene in high_rpf and gene in high_rna])

        if len(overlap) == 0:
            empty = pd.DataFrame(columns=BASE_COLUMNS + [
                "Sample", "RPF_Sample", "RNA_Sample", "RawRPF", "RNA_Count",
                "RNA_EffectiveLength", "RNA_RPKM", "NormalizedRPF", "IsValidCodon"
            ])
            summary = {
                "Sample": output_sample,
                "RPF_Sample": rpf_sample,
                "RNA_Sample": rna_sample,
                "HighRPFGenes": len(high_rpf),
                "HighRNAGenes": len(high_rna),
                "OverlapGenes": 0,
                "OutlierCount": len(outliers),
            }
            return empty, outliers, summary

        total_rna = float(np.nansum(rna[rna_sample].to_numpy(dtype=float)))
        effective_length = rna_gene_codons.loc[overlap] * 3.0
        gene_rna_count = rna_gene_count.loc[overlap].astype(float)
        if total_rna > 0:
            gene_rpkm = RPKM_SCALE * gene_rna_count / (total_rna * effective_length)
        else:
            gene_rpkm = pd.Series(np.nan, index=overlap, dtype=float)

        result = rpf.loc[rpf["name"].isin(overlap), BASE_COLUMNS + [rpf_sample]].copy()
        result.rename(columns={rpf_sample: "RawRPF"}, inplace=True)
        result["Sample"] = output_sample
        result["RPF_Sample"] = rpf_sample
        result["RNA_Sample"] = rna_sample
        result["RNA_Count"] = result["name"].map(gene_rna_count)
        result["RNA_EffectiveLength"] = result["name"].map(effective_length)
        result["RNA_RPKM"] = result["name"].map(gene_rpkm)
        result["NormalizedRPF"] = self._safe_divide(
            result["RawRPF"].to_numpy(dtype=float),
            result["RNA_RPKM"].to_numpy(dtype=float),
        )
        result["IsValidCodon"] = result["RawRPF"].fillna(0.0) > 0

        summary = {
            "Sample": output_sample,
            "RPF_Sample": rpf_sample,
            "RNA_Sample": rna_sample,
            "HighRPFGenes": int(len(high_rpf)),
            "HighRNAGenes": int(len(high_rna)),
            "OverlapGenes": int(len(overlap)),
            "OutlierCount": int(len(outliers)),
            "ValidCodonPositions": int(result["IsValidCodon"].sum()),
        }
        return result, outliers, summary

    def calculate_cdt(self) -> None:
        """Calculate CDT independently for all sample pairs."""
        workers = min(self.thread, len(self.sample_pairs))
        results = {}
        if workers <= 1:
            for pair in self.sample_pairs:
                results[pair[0]] = self._calculate_sample_pair(*pair)
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(self._calculate_sample_pair, *pair): pair[0]
                    for pair in self.sample_pairs
                }
                for future in as_completed(futures):
                    results[futures[future]] = future.result()

        position_tables = []
        outlier_tables = []
        summaries = []
        for output_sample, _, _ in self.sample_pairs:
            table, outlier, summary = results[output_sample]
            if not table.empty:
                position_tables.append(table)
            if not outlier.empty:
                outlier_tables.append(outlier)
            summaries.append(summary)

        if not position_tables:
            raise ValueError("No sample pair retained transcripts for CDT calculation.")

        self.position_cdt = pd.concat(position_tables, ignore_index=True)
        self.outliers = pd.concat(outlier_tables, ignore_index=True) if outlier_tables else pd.DataFrame()
        self.sample_summary = pd.DataFrame.from_records(summaries)
        self._summarize_gene_codon()
        self._summarize_codon()

    # Backward-compatible method names.
    calc_rpf_norm = calculate_cdt
    codon_decoding_time = lambda self: None

    def _summarize_gene_codon(self) -> None:
        """Summarize CDT components by transcript, codon and sample."""
        table = self.position_cdt.copy()
        grouped = table.groupby(["name", "codon", "Sample"], sort=False)
        self.gene_codon_cdt = grouped.agg(
            CodonCount=("RawRPF", "size"),
            ValidCodonCount=("IsValidCodon", "sum"),
            RPFCount=("RawRPF", "sum"),
            NormalizedRPFSum=("NormalizedRPF", "sum"),
        ).reset_index()
        valid = self.gene_codon_cdt["ValidCodonCount"].to_numpy(dtype=float)
        self.gene_codon_cdt["AbsoluteCDT"] = self._safe_divide(
            self.gene_codon_cdt["RPFCount"].to_numpy(dtype=float), valid
        )
        self.gene_codon_cdt["NormalizedCDT"] = self._safe_divide(
            self.gene_codon_cdt["NormalizedRPFSum"].to_numpy(dtype=float), valid
        )

    @staticmethod
    def _scale_series(values: pd.Series, method: str) -> pd.Series:
        """Scale one sample-wise codon metric."""
        values = values.astype(float)
        if method == "none":
            return values
        if method == "minmax":
            maximum = values.max(skipna=True)
            return values / maximum if pd.notna(maximum) and maximum != 0 else values * np.nan
        if method == "zscore":
            std = values.std(skipna=True)
            return (values - values.mean(skipna=True)) / std if pd.notna(std) and std != 0 else values * np.nan
        raise ValueError(f"Unsupported scale method: {method}")

    def _summarize_codon(self) -> None:
        """Summarize absolute and RNA-normalized CDT for each codon."""
        table = self.position_cdt.copy()
        grouped = table.groupby(["codon", "Sample"], sort=False)
        long = grouped.agg(
            CodonCount=("RawRPF", "size"),
            ValidCodonCount=("IsValidCodon", "sum"),
            RPFCount=("RawRPF", "sum"),
            NormalizedRPFSum=("NormalizedRPF", "sum"),
        ).reset_index()
        valid = long["ValidCodonCount"].to_numpy(dtype=float)
        long["AbsoluteCDT"] = self._safe_divide(long["RPFCount"].to_numpy(dtype=float), valid)
        long["NormalizedCDT"] = self._safe_divide(long["NormalizedRPFSum"].to_numpy(dtype=float), valid)
        long["RelativeCDT"] = long.groupby("Sample", sort=False)["AbsoluteCDT"].transform(
            lambda x: self._scale_series(x, self.scale)
        )
        long["NormalizedRelativeCDT"] = long.groupby("Sample", sort=False)["NormalizedCDT"].transform(
            lambda x: self._scale_series(x, self.scale)
        )

        annotation = self.codon_annotation.reset_index().rename(columns={"index": "codon"})
        self.codon_cdt = annotation.merge(long, on="codon", how="left")
        self.codon_cdt.sort_values(["Abbr", "codon", "Sample"], inplace=True)

    # ------------------------------------------------------------------
    # Output formatting
    # ------------------------------------------------------------------

    def _format_cdt_wide(self) -> pd.DataFrame:
        """Convert codon-level CDT results from long to wide format."""

        base_columns = [
            "codon",
            "AA",
            "Abbr",
        ]

        value_columns = [
            "CodonCount",
            "ValidCodonCount",
            "RPFCount",
            "NormalizedRPFSum",
            "AbsoluteCDT",
            "NormalizedCDT",
            "RelativeCDT",
            "NormalizedRelativeCDT",
        ]

        wide = self.codon_cdt.pivot_table(
            index=base_columns,
            columns="Sample",
            values=value_columns,
            aggfunc="first",
        )

        wide.columns = [
            "{}_{}".format(sample, metric)
            for metric, sample in wide.columns
        ]

        return wide.reset_index()

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def output_tables(self) -> None:
        """Write codon-level and optional detailed CDT tables."""

        wide_cdt = self._format_cdt_wide()

        out_cdt = self.output + "_codon_decoding_time.txt"
        wide_cdt.to_csv(out_cdt, sep="\t", index=False)
        self.output_files["cdt_table"] = out_cdt

        # Keep long-format CDT table for compatibility and debugging.
        out_cdt_long = self.output + "_codon_decoding_time.long.txt"
        self.codon_cdt.to_csv(out_cdt_long, sep="\t", index=False)
        self.output_files["cdt_long_table"] = out_cdt_long

        if self.output_all:
            out_position = self.output + "_cdt_position.txt"
            out_gene_codon = self.output + "_gene_codon_cdt.txt"
            self.position_cdt.to_csv(out_position, sep="\t", index=False)
            self.gene_codon_cdt.to_csv(out_gene_codon, sep="\t", index=False)
            self.output_files["position_table"] = out_position
            self.output_files["gene_codon_table"] = out_gene_codon

        if self.remove_outlier:
            out_outlier = self.output + "_cdt.outliers.txt"
            self.outliers.to_csv(out_outlier, sep="\t", index=False)
            self.output_files["outlier_table"] = out_outlier

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    def _metric_matrix(self, metric: str) -> pd.DataFrame:
        """Return a codon-by-sample matrix for one CDT metric."""
        matrix = self.codon_cdt.pivot_table(
            index="codon", columns="Sample", values=metric, aggfunc="first"
        )
        codon_order = self.codon_annotation.sort_values(["Abbr", "codon"]).index.tolist()
        sample_order = [pair[0] for pair in self.sample_pairs]
        return matrix.reindex(index=codon_order, columns=sample_order)

    @staticmethod
    def _adaptive_correlation_limits(values: np.ndarray) -> tuple[float, float]:
        """Return a focused color range for highly correlated samples."""
        values = np.asarray(values, dtype=float)
        offdiag = values[~np.eye(values.shape[0], dtype=bool)] if values.ndim == 2 else values
        offdiag = offdiag[np.isfinite(offdiag)]
        if offdiag.size == 0:
            return -1.0, 1.0
        low = float(np.quantile(offdiag, 0.02))
        high = float(np.quantile(offdiag, 0.98))
        span = max(high - low, 0.02)
        return max(-1.0, low - span * 0.15), 1.0

    def draw_cdt_corr(self) -> None:
        """Draw sample correlation based on RNA-normalized absolute CDT."""
        matrix = self._metric_matrix("NormalizedCDT")
        corr = matrix.corr(method="pearson", min_periods=3)
        out_txt = self.output + "_cdt_corr.txt"
        corr.to_csv(out_txt, sep="\t")
        self.output_files["correlation_table"] = out_txt

        values = corr.to_numpy(dtype=float)
        vmin, vmax = self._adaptive_correlation_limits(values)
        sample_count = max(1, corr.shape[0])
        size = max(5.5, sample_count * 0.48 + 3.0)
        fig, ax = plt.subplots(figsize=(size, size))
        image = ax.imshow(
            values,
            aspect="equal",
            interpolation="nearest",
            cmap=PLOT_CMAP,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_xticks(range(sample_count))
        ax.set_yticks(range(sample_count))
        ax.set_xticklabels(corr.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(corr.index, fontsize=8)
        ax.set_title("Codon decoding time correlation")
        ax.tick_params(length=0)

        if sample_count <= HEATMAP_ANNOTATION_MAX_SAMPLES:
            color_span = max(vmax - vmin, np.finfo(float).eps)
            for row_idx in range(values.shape[0]):
                for col_idx in range(values.shape[1]):
                    value = values[row_idx, col_idx]
                    if not np.isfinite(value):
                        continue
                    color_fraction = (value - vmin) / color_span
                    text_color = "white" if color_fraction < 0.18 or color_fraction > 0.82 else "#222222"
                    ax.text(
                        col_idx,
                        row_idx,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=text_color,
                    )

        cbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.03)
        cbar.set_label(f"Pearson correlation ({vmin:.2f} to {vmax:.2f})")
        fig.tight_layout()

        out_pdf = self.output + "_cdt_corrplot.pdf"
        out_png = self.output + "_cdt_corrplot.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["correlation_pdf"] = out_pdf
        self.output_files["correlation_png"] = out_png

    def _transform_plot_values(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """Transform plotted values without modifying output tables."""
        values = matrix.astype(float).clip(lower=0)
        if self.plot_transform == "none":
            return values
        if self.plot_transform == "sqrt":
            return np.sqrt(values)
        if self.plot_transform in {"log", "log1p"}:
            return np.log1p(values)
        if self.plot_transform == "log2":
            return np.log2(values + 1.0)
        if self.plot_transform == "log10":
            return np.log10(values + 1.0)
        raise ValueError(f"Unsupported plot transform: {self.plot_transform}")

    def draw_cdt_heat(self) -> None:
        """Draw absolute and RNA-normalized CDT heatmaps."""
        absolute = self._metric_matrix("AbsoluteCDT")
        normalized = self._metric_matrix("NormalizedCDT")
        absolute_plot = self._transform_plot_values(absolute)
        normalized_plot = self._transform_plot_values(normalized)

        sample_count = max(1, len(self.sample_pairs))
        figure_height = max(9.0, len(absolute_plot) * 0.22)
        figure_width = max(10.0, sample_count * 0.55 + 7.5)
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(figure_width, figure_height),
            gridspec_kw={"wspace": 0.28},
        )

        for ax, matrix, title in (
            (axes[0], absolute_plot, "Absolute CDT"),
            (axes[1], normalized_plot, "RNA-normalized CDT"),
        ):
            values = matrix.to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            vmax = float(np.quantile(finite, 0.98)) if finite.size else 1.0
            vmax = vmax if vmax > 0 else 1.0
            image = ax.imshow(
                values,
                aspect="auto",
                interpolation="nearest",
                cmap=PLOT_CMAP,
                vmin=0,
                vmax=vmax,
            )
            ax.set_title(title, fontsize=11)
            ax.set_xlabel("Sample")
            ax.set_xticks(range(matrix.shape[1]))
            ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
            labels = [f"{codon} [{self.codon_annotation.loc[codon, 'Abbr']}]" for codon in matrix.index]
            ax.set_yticks(range(matrix.shape[0]))
            ax.set_yticklabels(labels, fontsize=7)
            ax.tick_params(length=0)

            # Add numeric CDT values only when the sample count is small enough
            # to keep the heatmap readable.
            if sample_count <= HEATMAP_ANNOTATION_MAX_SAMPLES:
                color_span = max(vmax - 0.0, np.finfo(float).eps)
                for row_idx in range(values.shape[0]):
                    for col_idx in range(values.shape[1]):
                        value = values[row_idx, col_idx]
                        if not np.isfinite(value):
                            continue
                        color_fraction = (value - 0.0) / color_span
                        text_color = "white" if color_fraction < 0.18 or color_fraction > 0.82 else "#222222"
                        ax.text(
                            col_idx,
                            row_idx,
                            f"{value:.2f}",
                            ha="center",
                            va="center",
                            fontsize=5.5,
                            color=text_color,
                        )

            cbar = fig.colorbar(image, ax=ax, shrink=0.55, pad=0.02)
            cbar.set_label(title)

        axes[0].set_ylabel("Codon [amino acid]")
        axes[1].set_ylabel("")

        out_pdf = self.output + "_cdt_heatmap.pdf"
        out_png = self.output + "_cdt_heatmap.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["heatmap_pdf"] = out_pdf
        self.output_files["heatmap_png"] = out_png

    def draw_cdt_rank(self) -> None:
        """Draw sample-wise CDT rank profiles with codons on the x-axis."""
        matrix = self._metric_matrix("NormalizedCDT")
        long_df = matrix.reset_index().melt(
            id_vars="codon", var_name="Sample", value_name="CDT"
        )
        long_df = long_df.merge(
            self.codon_annotation.reset_index().rename(columns={"index": "codon"}),
            on="codon",
            how="left",
        )

        # Each x-axis tick is one codon labeled as "Codon[AminoAcid]", so the
        # panel width must be large enough to host all rotated labels side by
        # side. `rankplot_ncol` controls how many panels are placed per row
        # (e.g. 1 wide panel per row, or 2 narrower panels for a more compact
        # layout). The label-width factor is kept small so that larger tick
        # fonts do not inflate the panel width too much.
        ncols = max(1, int(self.rankplot_ncol))
        ncols = min(ncols, max(1, matrix.shape[1]))
        nrows = int(np.ceil(matrix.shape[1] / ncols))
        n_codons = int(len(matrix.index))
        tick_fontsize = 14 if ncols == 1 else 12
        label_width = 3 * 0.42 * tick_fontsize / 72.0 + 0.01
        panel_width = max(5.5, 2.0 + n_codons * label_width)
        panel_height = 3.8
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(panel_width * ncols, panel_height * nrows),
            squeeze=False,
        )

        for ax, sample in zip(axes.ravel(), matrix.columns):
            data = long_df.loc[long_df["Sample"] == sample, :].dropna()
            data = data.sort_values("CDT", ascending=True).reset_index(drop=True)
            x = np.arange(len(data))
            ax.scatter(
                x,
                data["CDT"],
                s=15,
                alpha=0.8,
                edgecolors="none",
                label="_nolegend_",
            )
            top = data.tail(min(6, len(data)))
            ax.scatter(
                top.index.to_numpy(),
                top["CDT"],
                s=24,
                alpha=0.95,
                color="#c0392b",
                edgecolors="none",
                zorder=3,
                label="Top 6 decoded codons",
            )
            ax.set_xticks(x)
            ax.set_xticklabels(
                data["codon"] + "[" + data["Abbr"].astype(str) + "]",
                rotation=90,
                ha="center",
                fontsize=tick_fontsize,
            )
            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlim(-0.5, len(data) - 0.5)
            ax.set_title(sample, fontsize=14)
            ax.set_xlabel("Codon [amino acid]", fontsize=13)
            ax.set_ylabel("RNA-normalized CDT", fontsize=13)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(axis="y", linewidth=0.4, alpha=0.25)

        for ax in axes.ravel()[matrix.shape[1] :]:
            ax.set_axis_off()

        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                frameon=False,
                loc="upper center",
                ncol=len(handles),
                fontsize=12,
            )

        out_pdf = self.output + "_cdt_rankplot.pdf"
        out_png = self.output + "_cdt_rankplot.png"
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["rankplot_pdf"] = out_pdf
        self.output_files["rankplot_png"] = out_png

    # ------------------------------------------------------------------
    # Summary and workflow
    # ------------------------------------------------------------------

    def write_summary(self) -> None:
        """Write a machine-readable CDT summary JSON."""
        path = self.output + "_cdt.summary.json"
        self.output_files["summary_json"] = path
        summary = OrderedDict(
            [
                ("tool", "rpf_CDT"),
                ("version", "0.2.8-dev.001"),
                ("rpf_file", os.path.abspath(self.rpf_file)),
                ("rna_file", os.path.abspath(self.rna_file)),
                ("rpf_format", self.rpf_format),
                ("rna_format", self.rna_format),
                ("sample_pairs", [
                    {"sample": out, "rpf_sample": rpf, "rna_sample": rna}
                    for out, rpf, rna in self.sample_pairs
                ]),
                ("parameters", OrderedDict([
                    ("site", self.site), ("frame", self.frame),
                    ("min_rpf", self.rpf_num), ("min_rna", self.rna_num),
                    ("tis", self.tis), ("tts", self.tts),
                    ("scale", self.scale), ("plot_transform", self.plot_transform),
                    ("rankplot_ncol", self.rankplot_ncol),
                    ("thread", self.thread), ("remove_outlier", self.remove_outlier),
                    ("outlier_iqr", self.outlier_iqr),
                    ("outlier_window", self.outlier_window),
                    ("outlier_local_fold", self.outlier_local_fold),
                ])),
                ("sample_summary", self.sample_summary.to_dict(orient="records")),
                ("output_files", self.output_files),
            ]
        )
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

    def run(self) -> None:
        """Run the complete CDT workflow."""
        self.import_density()
        self.calculate_cdt()
        self.output_tables()
        self.draw_cdt_corr()
        self.draw_cdt_heat()
        self.draw_cdt_rank()
        self.write_summary()
