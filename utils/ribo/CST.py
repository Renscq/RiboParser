#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-09-10
# Version: dev002
# Function: Calculate sample-specific codon selection time from paired RPF and RNA density data.
# Input: RPF and RNA density files in JSONL or TXT format and an optional transcript filter.
# Output: Codon selection time tables, outlier records, summary JSON, correlation, heatmap, rank, and convergence plots.

"""Core functions for RiboParser codon selection time analysis.

RPF and RNA density files are imported through ``RPFs.RPFData``. RPF/RNA
samples are paired by name when possible and otherwise by column order. Genes
are filtered independently for every sample pair.

For codon ``c`` and sample pair ``s``, the initial codon selection time is::

    CST_0[c, s] = RPF_proportion[c, s] / RNA_proportion[c, s]

For each iteration, gene-specific elongation and initiation rates are estimated
from the current CST vector, and RNA codon proportions are reweighted by the
estimated initiation rates before the CST vector is updated. The implementation
uses gene-by-codon matrices and NumPy matrix multiplication instead of legacy
per-gene/per-codon Python loops.
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

# Use Arial (with common fallbacks when it is unavailable on the system) for
# every figure produced by this module.
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
    }
)


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
STOP_CODONS = {"TAA", "TAG", "TGA"}
HEATMAP_ANNOTATION_MAX_SAMPLES = 12
PLOT_CMAP = "RdBu_r"


class CodonSelectiveTime(object):
    """Calculate sample-specific codon selection time.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_CST``.
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
        self.times = args.times
        self.tolerance = args.tolerance
        self.scale = args.scale
        self.plot_transform = args.plot_transform
        self.thread = max(1, int(args.thread))
        self.output_all = args.all
        self.rankplot_ncol = max(1, int(getattr(args, "rankplot_ncol", 1)))

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

        self.cst = None
        self.iterative_cst = None
        self.gene_metrics = None
        self.outliers = pd.DataFrame()
        self.sample_summary = pd.DataFrame()
        self.output_files = OrderedDict()

        _, codon_table = RPFs.codon_table()
        self.codon_annotation = codon_table.loc[codon_table["Abbr"] != "*", :].copy()
        self.codon_order = self.codon_annotation.sort_values(["Abbr", "codon"]).index.tolist()

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
        # RNA density is not shifted by a ribosomal E/P/A-site offset.
        self.rna_table = self.rna_data.get_frame(frame=self.frame)

        self._validate_density_table(self.rpf_table, self.rpf_samples, "RPF")
        self._validate_density_table(self.rna_table, self.rna_samples, "RNA")

        self.rpf_table = self.rpf_table.loc[
            (self.rpf_table["region"] == "cds")
            & (~self.rpf_table["codon"].isin(STOP_CODONS)),
            BASE_COLUMNS + self.rpf_samples,
        ].copy()
        self.rna_table = self.rna_table.loc[
            (self.rna_table["region"] == "cds")
            & (~self.rna_table["codon"].isin(STOP_CODONS)),
            BASE_COLUMNS + self.rna_samples,
        ].copy()

        for sample in self.rpf_samples:
            self.rpf_table[sample] = pd.to_numeric(
                self.rpf_table[sample], errors="coerce"
            ).fillna(0.0).astype(float)
        for sample in self.rna_samples:
            self.rna_table[sample] = pd.to_numeric(
                self.rna_table[sample], errors="coerce"
            ).fillna(0.0).astype(float)

        if self.rpf_table.empty or self.rna_table.empty:
            raise ValueError("No CDS codons remain after RPF/RNA filtering.")

        print(
            "Imported RPF format={rpf_fmt}, RNA format={rna_fmt}, sample pairs={pairs}.".format(
                rpf_fmt=self.rpf_format,
                rna_fmt=self.rna_format,
                pairs=len(self.sample_pairs),
            ),
            flush=True,
        )
        for output_sample, rpf_sample, rna_sample in self.sample_pairs:
            print(f"Pair {output_sample}: RPF={rpf_sample}, RNA={rna_sample}.", flush=True)

    @staticmethod
    def _validate_density_table(table: pd.DataFrame, samples: list[str], label: str) -> None:
        """Validate an imported codon-level density table."""
        missing = [column for column in BASE_COLUMNS + samples if column not in table.columns]
        if missing:
            raise ValueError(f"{label} density table is missing column(s): {', '.join(missing)}")

    @staticmethod
    def _pair_samples(rpf_samples: list[str], rna_samples: list[str]) -> list[tuple[str, str, str]]:
        """Pair RPF and RNA samples by name or, when necessary, by order."""
        if not rpf_samples or not rna_samples:
            raise ValueError("RPF and RNA density files must both contain samples.")

        if set(rpf_samples) == set(rna_samples):
            return [(sample, sample, sample) for sample in rpf_samples]

        if len(rpf_samples) != len(rna_samples):
            raise ValueError(
                "RPF and RNA sample names differ and sample counts are unequal. "
                "Use matching names or equal ordered sample sets."
            )

        print(
            "Warning: RPF and RNA sample names differ; samples are paired by column order.",
            flush=True,
        )
        return [
            (rpf_sample, rpf_sample, rna_sample)
            for rpf_sample, rna_sample in zip(rpf_samples, rna_samples)
        ]

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
        """Calculate a symmetric local mean excluding the focal codon."""
        values = np.asarray(values, dtype=float)
        result = np.full(values.size, np.nan, dtype=float)
        if values.size == 0:
            return result

        finite = np.isfinite(values)
        clean = np.where(finite, values, 0.0)
        prefix_sum = np.concatenate(([0.0], np.cumsum(clean)))
        prefix_n = np.concatenate(([0], np.cumsum(finite.astype(int))))
        idx = np.arange(values.size)
        left = np.maximum(0, idx - int(window))
        right = np.minimum(values.size, idx + int(window) + 1)
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
        if not self.remove_outlier:
            return work, pd.DataFrame()

        cutoff = self._robust_upper_cutoff(work[sample], self.outlier_iqr)
        if not np.isfinite(cutoff):
            return work, pd.DataFrame()

        candidate_idx = work.index[work[sample] > cutoff]
        if len(candidate_idx) == 0:
            return work, pd.DataFrame()

        outlier_records = []
        candidate_genes = work.loc[candidate_idx, "name"].drop_duplicates().tolist()
        for gene in candidate_genes:
            gene_idx = work.index[work["name"] == gene]
            gene_values = work.loc[gene_idx, sample].to_numpy(dtype=float)
            background = self._local_background(gene_values, self.outlier_window)
            index_to_pos = {idx: pos for pos, idx in enumerate(gene_idx)}
            gene_candidates = [idx for idx in candidate_idx if idx in index_to_pos]

            for idx in gene_candidates:
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
    # CST calculation
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
        """Divide arrays and return NaN for invalid denominators."""
        numerator = np.asarray(numerator, dtype=float)
        denominator = np.asarray(denominator, dtype=float)
        result = np.full(numerator.shape, np.nan, dtype=float)
        valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0)
        np.divide(numerator, denominator, out=result, where=valid)
        return result

    @staticmethod
    def _scale_values(values: pd.Series, method: str) -> pd.Series:
        """Scale one sample's final CST values."""
        values = values.astype(float)
        if method == "none":
            return values
        if method == "minmax":
            maximum = values.max(skipna=True)
            return values / maximum if np.isfinite(maximum) and maximum > 0 else values * np.nan
        if method == "zscore":
            mean = values.mean(skipna=True)
            std = values.std(skipna=True, ddof=0)
            return (values - mean) / std if np.isfinite(std) and std > 0 else values * np.nan
        raise ValueError(f"Unsupported scale method: {method}")

    def _calculate_sample_pair(
        self,
        output_sample: str,
        rpf_sample: str,
        rna_sample: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
        """Calculate CST for one paired RPF/RNA sample."""
        rpf = self.rpf_table.loc[:, BASE_COLUMNS + [rpf_sample]].copy()
        rna = self.rna_table.loc[:, BASE_COLUMNS + [rna_sample]].copy()
        rpf, outliers = self._remove_sample_outliers(rpf, rpf_sample, output_sample)

        rpf_gene_sum = rpf.groupby("name", sort=False)[rpf_sample].sum(min_count=1)
        rna_gene_sum = rna.groupby("name", sort=False)[rna_sample].sum(min_count=1)

        high_rpf = set(rpf_gene_sum.index[rpf_gene_sum >= float(self.rpf_num)])
        high_rna = set(rna_gene_sum.index[rna_gene_sum >= float(self.rna_num)])
        overlap = [gene for gene in rpf_gene_sum.index if gene in high_rpf and gene in high_rna]

        if not overlap:
            empty = pd.DataFrame()
            summary = {
                "Sample": output_sample,
                "RPF_Sample": rpf_sample,
                "RNA_Sample": rna_sample,
                "HighGeneCount": 0,
                "OutlierCount": int(len(outliers)),
                "IterationsCompleted": 0,
                "Converged": False,
                "FinalDelta": None,
            }
            return empty, empty, empty, outliers, summary

        rpf = rpf.loc[rpf["name"].isin(overlap), :].copy()
        rna = rna.loc[rna["name"].isin(overlap), :].copy()

        # Use RPF codon coordinates as the sequence definition. RNA density is
        # aggregated independently by the same transcript/codon labels.
        gene_codon_count = (
            rpf.groupby(["name", "codon"], sort=False)
            .size()
            .unstack(fill_value=0)
            .reindex(index=overlap, columns=self.codon_order, fill_value=0)
            .astype(float)
        )
        gene_length = gene_codon_count.sum(axis=1).to_numpy(dtype=float)

        rpf_gene_codon = (
            rpf.groupby(["name", "codon"], sort=False)[rpf_sample]
            .sum(min_count=1)
            .unstack(fill_value=0)
            .reindex(index=overlap, columns=self.codon_order, fill_value=0)
            .astype(float)
        )
        rna_gene_codon = (
            rna.groupby(["name", "codon"], sort=False)[rna_sample]
            .sum(min_count=1)
            .unstack(fill_value=0)
            .reindex(index=overlap, columns=self.codon_order, fill_value=0)
            .astype(float)
        )

        rpf_gene_values = rpf_gene_codon.sum(axis=1).to_numpy(dtype=float)
        rna_gene_values = rna_gene_codon.sum(axis=1).to_numpy(dtype=float)
        translation_efficiency = self._safe_ratio(rpf_gene_values, rna_gene_values)

        rpf_codon_sum = rpf_gene_codon.sum(axis=0).to_numpy(dtype=float)
        rna_codon_sum = rna_gene_codon.sum(axis=0).to_numpy(dtype=float)
        total_rpf = float(np.nansum(rpf_codon_sum))
        total_rna = float(np.nansum(rna_codon_sum))
        rpf_proportion = rpf_codon_sum / total_rpf if total_rpf > 0 else np.full(len(self.codon_order), np.nan)
        rna_proportion = rna_codon_sum / total_rna if total_rna > 0 else np.full(len(self.codon_order), np.nan)
        cst = self._safe_ratio(rpf_proportion, rna_proportion)

        iteration_records = []
        convergence_records = []

        def append_iteration(iteration: int, values: np.ndarray, rna_prop: np.ndarray, delta: float | None) -> None:
            for codon, value, rp, mp in zip(self.codon_order, values, rpf_proportion, rna_prop):
                iteration_records.append(
                    {
                        "Codon": codon,
                        "Sample": output_sample,
                        "Iteration": int(iteration),
                        "RPFProportion": float(rp) if np.isfinite(rp) else np.nan,
                        "AdjustedRNAProportion": float(mp) if np.isfinite(mp) else np.nan,
                        "AbsoluteCST": float(value) if np.isfinite(value) else np.nan,
                    }
                )
            convergence_records.append(
                {
                    "Sample": output_sample,
                    "Iteration": int(iteration),
                    "MaxAbsoluteDelta": delta,
                }
            )

        append_iteration(0, cst, rna_proportion, None)
        converged = False
        final_delta = np.nan
        iterations_completed = 0
        count_matrix = gene_codon_count.to_numpy(dtype=float)
        rna_matrix = rna_gene_codon.to_numpy(dtype=float)

        for iteration in range(1, int(self.times) + 1):
            cst_for_calc = np.where(np.isfinite(cst), cst, 0.0)
            weighted_time = count_matrix @ cst_for_calc
            elongation_rate = self._safe_ratio(gene_length, weighted_time)
            initiation_rate = elongation_rate * translation_efficiency
            initiation_rate = np.where(np.isfinite(initiation_rate), initiation_rate, 0.0)

            adjusted_rna_codon = initiation_rate @ rna_matrix
            adjusted_total = float(np.nansum(adjusted_rna_codon))
            adjusted_rna_prop = (
                adjusted_rna_codon / adjusted_total
                if adjusted_total > 0
                else np.full(len(self.codon_order), np.nan)
            )
            new_cst = self._safe_ratio(rpf_proportion, adjusted_rna_prop)

            finite = np.isfinite(new_cst) & np.isfinite(cst)
            final_delta = float(np.max(np.abs(new_cst[finite] - cst[finite]))) if finite.any() else np.nan
            cst = new_cst
            iterations_completed = iteration
            append_iteration(iteration, cst, adjusted_rna_prop, final_delta)

            if np.isfinite(final_delta) and final_delta <= float(self.tolerance):
                converged = True
                break

        final_df = pd.DataFrame(
            {
                "Codon": self.codon_order,
                "Sample": output_sample,
                "CodonCount": gene_codon_count.sum(axis=0).to_numpy(dtype=float),
                "ValidCodonCount": (rpf_gene_codon > 0).sum(axis=0).to_numpy(dtype=float),
                "RPFCount": rpf_codon_sum,
                "RNACount": rna_codon_sum,
                "RPFProportion": rpf_proportion,
                "InitialRNAProportion": rna_proportion,
                "AbsoluteCST": cst,
                "IterationsCompleted": iterations_completed,
                "Converged": converged,
            }
        )
        final_df["RelativeCST"] = self._scale_values(final_df["AbsoluteCST"], self.scale)

        gene_metrics = pd.DataFrame(
            {
                "name": overlap,
                "Sample": output_sample,
                "RPFCount": rpf_gene_values,
                "RNACount": rna_gene_values,
                "TranslationEfficiency": translation_efficiency,
            }
        )
        if iterations_completed > 0:
            cst_for_calc = np.where(np.isfinite(cst), cst, 0.0)
            weighted_time = count_matrix @ cst_for_calc
            elongation_rate = self._safe_ratio(gene_length, weighted_time)
            initiation_rate = elongation_rate * translation_efficiency
            gene_metrics["ElongationRate"] = elongation_rate
            gene_metrics["InitiationRate"] = initiation_rate
        else:
            gene_metrics["ElongationRate"] = np.nan
            gene_metrics["InitiationRate"] = np.nan

        summary = {
            "Sample": output_sample,
            "RPF_Sample": rpf_sample,
            "RNA_Sample": rna_sample,
            "HighGeneCount": int(len(overlap)),
            "OutlierCount": int(len(outliers)),
            "IterationsCompleted": int(iterations_completed),
            "Converged": bool(converged),
            "FinalDelta": float(final_delta) if np.isfinite(final_delta) else None,
        }

        return (
            final_df,
            pd.DataFrame.from_records(iteration_records),
            gene_metrics,
            outliers,
            {**summary, "Convergence": pd.DataFrame.from_records(convergence_records)},
        )

    def calculate_cst(self) -> None:
        """Calculate CST independently for all sample pairs."""
        if self.rpf_table is None or self.rna_table is None:
            raise ValueError("Density data have not been imported.")

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

        final_tables = []
        iterative_tables = []
        gene_tables = []
        outlier_tables = []
        summaries = []
        convergence_tables = []

        for output_sample, _, _ in self.sample_pairs:
            final_df, iterative_df, gene_df, outlier_df, summary = results[output_sample]
            if not final_df.empty:
                final_tables.append(final_df)
            if not iterative_df.empty:
                iterative_tables.append(iterative_df)
            if not gene_df.empty:
                gene_tables.append(gene_df)
            if not outlier_df.empty:
                outlier_tables.append(outlier_df)
            convergence = summary.pop("Convergence")
            if not convergence.empty:
                convergence_tables.append(convergence)
            summaries.append(summary)
            print(
                "Sample {sample}: genes={genes:,}, iterations={iterations}, "
                "converged={converged}, outliers={outliers:,}.".format(
                    sample=output_sample,
                    genes=summary["HighGeneCount"],
                    iterations=summary["IterationsCompleted"],
                    converged=summary["Converged"],
                    outliers=summary["OutlierCount"],
                ),
                flush=True,
            )

        if not final_tables:
            raise ValueError("No CST result was generated. Check expression thresholds and sample pairing.")

        self.cst = pd.concat(final_tables, ignore_index=True)
        annotation = self.codon_annotation.reset_index().rename(columns={"codon": "Codon"})
        self.cst = self.cst.merge(annotation[["Codon", "AA", "Abbr"]], on="Codon", how="left")
        self.cst.sort_values(["Abbr", "Codon", "Sample"], inplace=True, ignore_index=True)
        self.iterative_cst = pd.concat(iterative_tables, ignore_index=True)
        self.gene_metrics = pd.concat(gene_tables, ignore_index=True) if gene_tables else pd.DataFrame()
        self.outliers = pd.concat(outlier_tables, ignore_index=True) if outlier_tables else pd.DataFrame()
        self.sample_summary = pd.DataFrame.from_records(summaries)
        self.convergence = pd.concat(convergence_tables, ignore_index=True) if convergence_tables else pd.DataFrame()

    # Backward-compatible method names.
    calc_cst = calculate_cst

    def format_cst_results(self) -> pd.DataFrame:
        """Convert the final long-format CST table to codon-wide format."""
        if self.cst is None:
            raise ValueError("CST has not been calculated.")

        value_columns = [
            "CodonCount",
            "ValidCodonCount",
            "RPFCount",
            "RNACount",
            "RPFProportion",
            "InitialRNAProportion",
            "AbsoluteCST",
            "IterationsCompleted",
            "Converged",
            "RelativeCST",
        ]
        sample_order = [pair[0] for pair in self.sample_pairs]
        annotation = self.codon_annotation.reset_index().rename(columns={"codon": "Codon"})
        annotation = annotation[["Codon", "AA", "Abbr"]].copy()
        annotation["_codon_order"] = pd.Categorical(
            annotation["Codon"],
            categories=self.codon_order,
            ordered=True,
        )

        wide_parts = [annotation]
        for sample in sample_order:
            sample_table = self.cst.loc[
                self.cst["Sample"] == sample,
                ["Codon"] + value_columns,
            ].copy()
            rename_map = {
                metric: "{}_{}".format(sample, metric)
                for metric in value_columns
            }
            sample_table.rename(columns=rename_map, inplace=True)
            wide_parts.append(sample_table)

        wide = wide_parts[0]
        for sample_table in wide_parts[1:]:
            wide = wide.merge(sample_table, on="Codon", how="left", validate="one_to_one")

        wide.sort_values("_codon_order", inplace=True)
        wide.drop(columns="_codon_order", inplace=True)
        wide.reset_index(drop=True, inplace=True)
        return wide

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def output_tables(self) -> None:
        """Write CST result tables."""
        if self.cst is None:
            raise ValueError("CST has not been calculated.")

        final_file = self.output + "_codon_selection_time.txt"
        self.format_cst_results().to_csv(final_file, sep="\t", index=False)
        self.output_files["cst_table"] = final_file

        iterative_file = self.output + "_iterative_codon_selection_time.txt"
        self.iterative_cst.to_csv(iterative_file, sep="\t", index=False)
        self.output_files["iterative_cst_table"] = iterative_file

        if self.output_all and self.gene_metrics is not None:
            gene_file = self.output + "_cst_gene_metrics.txt"
            self.gene_metrics.to_csv(gene_file, sep="\t", index=False)
            self.output_files["gene_metrics"] = gene_file

        if self.remove_outlier:
            outlier_file = self.output + "_cst.outliers.txt"
            self.outliers.to_csv(outlier_file, sep="\t", index=False)
            self.output_files["outlier_table"] = outlier_file

    output_cst = output_tables

    def write_summary(self) -> None:
        """Write a machine-readable CST summary JSON."""
        summary_file = self.output + "_cst.summary.json"
        self.output_files["summary_json"] = summary_file
        summary = OrderedDict(
            [
                ("tool", "rpf_CST"),
                ("version", "dev002"),
                ("input_rpf", os.path.abspath(self.rpf_file)),
                ("input_rna", os.path.abspath(self.rna_file)),
                ("rpf_format", self.rpf_format),
                ("rna_format", self.rna_format),
                ("sample_pairs", [
                    {"sample": out, "rpf_sample": rpf, "rna_sample": rna}
                    for out, rpf, rna in self.sample_pairs
                ]),
                ("parameters", OrderedDict(
                    [
                        ("min_rpf", self.rpf_num),
                        ("min_rna", self.rna_num),
                        ("site", self.site),
                        ("frame", self.frame),
                        ("tis", self.tis),
                        ("tts", self.tts),
                        ("times", self.times),
                        ("tolerance", self.tolerance),
                        ("scale", self.scale),
                        ("thread", self.thread),
                        ("rankplot_ncol", self.rankplot_ncol),
                        ("remove_outlier", self.remove_outlier),
                        ("outlier_iqr", self.outlier_iqr),
                        ("outlier_window", self.outlier_window),
                        ("outlier_local_fold", self.outlier_local_fold),
                    ]
                )),
                ("sample_summary", self.sample_summary.to_dict(orient="records")),
                ("output_files", self.output_files),
            ]
        )
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_name(name: str) -> str:
        """Return a file-safe sample name."""
        return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_") or "sample"

    def _transform_values(self, values: pd.DataFrame | pd.Series):
        """Transform positive CST values for plotting only."""
        transformed = values.astype(float).replace([np.inf, -np.inf], np.nan)
        if self.plot_transform == "none":
            return transformed
        clipped = transformed.clip(lower=0)
        if self.plot_transform == "sqrt":
            return np.sqrt(clipped)
        if self.plot_transform in {"log", "log1p"}:
            return np.log1p(clipped)
        if self.plot_transform == "log2":
            return np.log2(clipped + 1.0)
        if self.plot_transform == "log10":
            return np.log10(clipped + 1.0)
        raise ValueError(f"Unsupported plot transform: {self.plot_transform}")

    def _metric_matrix(self, metric: str) -> pd.DataFrame:
        """Return a codon-by-sample matrix for one CST metric."""
        matrix = self.cst.pivot_table(
            index="Codon", columns="Sample", values=metric, aggfunc="first"
        )
        sample_order = [pair[0] for pair in self.sample_pairs]
        return matrix.reindex(index=self.codon_order, columns=sample_order)

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

    def draw_cst_corr(self) -> None:
        """Draw sample correlation based on absolute CST values."""
        matrix = self._metric_matrix("AbsoluteCST")
        corr = matrix.corr(method="pearson", min_periods=3)
        corr_file = self.output + "_cst_corr.txt"
        corr.to_csv(corr_file, sep="\t")
        self.output_files["correlation_table"] = corr_file

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
        ax.set_title("Codon selection time correlation")
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

        out_pdf = self.output + "_cst_corrplot.pdf"
        out_png = self.output + "_cst_corrplot.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["correlation_pdf"] = out_pdf
        self.output_files["correlation_png"] = out_png

    def draw_cst_heat(self) -> None:
        """Draw absolute and relative CST heatmaps."""
        absolute = self._metric_matrix("AbsoluteCST")
        relative = self._metric_matrix("RelativeCST")
        absolute_plot = self._transform_values(absolute)

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
            (axes[0], absolute_plot, "Absolute CST"),
            (axes[1], relative, f"Relative CST ({self.scale})"),
        ):
            values = matrix.to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            if self.scale == "zscore" and title.startswith("Relative"):
                bound = float(np.quantile(np.abs(finite), 0.98)) if finite.size else 1.0
                bound = bound if bound > 0 else 1.0
                vmin, vmax = -bound, bound
            else:
                vmax = float(np.quantile(finite, 0.98)) if finite.size else 1.0
                vmax = vmax if vmax > 0 else 1.0
                vmin = 0.0
            image = ax.imshow(
                values,
                aspect="auto",
                interpolation="nearest",
                cmap=PLOT_CMAP,
                vmin=vmin,
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

            # Add numeric CST values only when the sample count is small enough
            # to keep the heatmap readable.
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
                            fontsize=5.5,
                            color=text_color,
                        )

            cbar = fig.colorbar(image, ax=ax, shrink=0.55, pad=0.02)
            cbar.set_label(title)

        axes[0].set_ylabel("Codon [amino acid]")
        axes[1].set_ylabel("")

        out_pdf = self.output + "_cst_heatmap.pdf"
        out_png = self.output + "_cst_heatmap.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["heatmap_pdf"] = out_pdf
        self.output_files["heatmap_png"] = out_png

    def draw_cst_rank(self) -> None:
        """Draw sample-wise CST rank profiles with codons on the x-axis."""
        matrix = self._metric_matrix("AbsoluteCST")
        long_df = matrix.reset_index().melt(
            id_vars="Codon", var_name="Sample", value_name="CST"
        )
        long_df = long_df.merge(
            self.codon_annotation.reset_index().rename(columns={"codon": "Codon"}),
            on="Codon",
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
            data = data.sort_values("CST", ascending=True).reset_index(drop=True)
            x = np.arange(len(data))
            ax.scatter(
                x,
                data["CST"],
                s=15,
                alpha=0.8,
                edgecolors="none",
                label="_nolegend_",
            )
            top = data.tail(min(6, len(data)))
            ax.scatter(
                top.index.to_numpy(),
                top["CST"],
                s=24,
                alpha=0.95,
                color="#c0392b",
                edgecolors="none",
                zorder=3,
                label="Top 6 selected codons",
            )
            ax.set_xticks(x)
            ax.set_xticklabels(
                data["Codon"] + "[" + data["Abbr"].astype(str) + "]",
                rotation=90,
                ha="center",
                fontsize=tick_fontsize,
            )
            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlim(-0.5, len(data) - 0.5)
            ax.set_title(sample, fontsize=14)
            ax.set_xlabel("Codon [amino acid]", fontsize=13)
            ax.set_ylabel("Absolute CST", fontsize=13)
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

        out_pdf = self.output + "_cst_rankplot.pdf"
        out_png = self.output + "_cst_rankplot.png"
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["rank_pdf"] = out_pdf
        self.output_files["rank_png"] = out_png

    def draw_cst_convergence(self) -> None:
        """Draw convergence trajectories for iterative CST estimation."""
        if self.convergence is None or self.convergence.empty:
            return

        fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=300)
        for sample, group in self.convergence.groupby("Sample", sort=False):
            group = group.loc[group["Iteration"] > 0].copy()
            if group.empty:
                continue
            ax.plot(group["Iteration"], group["MaxAbsoluteDelta"], marker="o", linewidth=1.2, markersize=3, label=sample)
        ax.axhline(self.tolerance, linestyle="--", linewidth=0.9, color="#333333", label="Tolerance")
        ax.set_yscale("log")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Maximum absolute CST change")
        ax.set_title("CST iteration convergence")
        ax.grid(axis="y", linewidth=0.4, alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=8, ncol=min(4, max(1, len(self.sample_pairs))))
        fig.tight_layout()

        out_pdf = self.output + "_cst_convergence.pdf"
        out_png = self.output + "_cst_convergence.png"
        fig.savefig(out_pdf, bbox_inches="tight")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        self.output_files["convergence_pdf"] = out_pdf
        self.output_files["convergence_png"] = out_png
