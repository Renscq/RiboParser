#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-09
# Version: dev001
# Function: Summarize and write pairwise SeRP peak-overlap analysis results.
# Input: PeakOverlapAnalyzer results from two SeRP conditions.
# Output: Summary, length, relationship, shared-cluster, shared-peak, and specific-peak tables.

"""Output and summary helpers for pairwise SeRP peak-overlap analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from utils.serp.Overlap import PeakOverlapAnalyzer


class SeRPOverlapOutput:
    """Summarize and write one completed SeRP peak-overlap comparison.

    Parameters
    ----------
    analysis : PeakOverlapAnalyzer
        Completed overlap analysis containing normalized peaks, interval
        relationships, shared graph clusters, and condition-specific peaks.
    """

    def __init__(self, analysis: PeakOverlapAnalyzer) -> None:
        """Initialize output tables from one overlap analysis."""
        self.analysis = analysis
        self.summary = pd.DataFrame()
        self.length_summary = pd.DataFrame()
        self.length_compare = pd.DataFrame()

    @staticmethod
    def _public_columns(table: pd.DataFrame) -> pd.DataFrame:
        """Remove private workflow columns before writing one result table."""
        private_columns = [column for column in table.columns if column.startswith("_")]
        return table.drop(columns=private_columns, errors="ignore")

    @staticmethod
    def _distribution_row(
        sample: str,
        category: str,
        length_type: str,
        values: pd.Series,
    ) -> dict[str, object]:
        """Summarize one peak-length distribution."""
        numeric = pd.to_numeric(values, errors="coerce").dropna().astype(float)
        if numeric.empty:
            return {
                "sample": sample,
                "category": category,
                "length_type": length_type,
                "n": 0,
                "mean": np.nan,
                "median": np.nan,
                "std": np.nan,
                "min": np.nan,
                "q1": np.nan,
                "q3": np.nan,
                "max": np.nan,
            }

        return {
            "sample": sample,
            "category": category,
            "length_type": length_type,
            "n": int(len(numeric)),
            "mean": float(numeric.mean()),
            "median": float(numeric.median()),
            "std": float(numeric.std(ddof=1)) if len(numeric) > 1 else 0.0,
            "min": float(numeric.min()),
            "q1": float(numeric.quantile(0.25)),
            "q3": float(numeric.quantile(0.75)),
            "max": float(numeric.max()),
        }

    @staticmethod
    def _unique_count(table: pd.DataFrame, column: str) -> int:
        """Count unique non-missing values from one column."""
        if column not in table.columns or table.empty:
            return 0
        values = table[column].replace("-", pd.NA).dropna().astype(str)
        return int(values.nunique())

    def build_length_summary(self) -> pd.DataFrame:
        """Summarize core and matched interval lengths by overlap class."""
        rows: list[dict[str, object]] = []
        for sample_name, table in (
            (self.analysis.name_a, self.analysis.annotated_a),
            (self.analysis.name_b, self.analysis.annotated_b),
        ):
            categories = {
                "all": table,
                "shared": table.loc[table["overlap_status"] == "shared"],
                "specific": table.loc[table["overlap_status"] == "specific"],
            }
            for category, subset in categories.items():
                rows.append(
                    self._distribution_row(
                        sample_name,
                        category,
                        "core_length",
                        subset["core_length_used"],
                    )
                )
                rows.append(
                    self._distribution_row(
                        sample_name,
                        category,
                        "match_span_length",
                        subset["match_span_length"],
                    )
                )
        self.length_summary = pd.DataFrame(rows)
        return self.length_summary

    def build_summary(self) -> pd.DataFrame:
        """Build comparison-level peak counts and input-filtering QC metrics."""
        annotated_a = self.analysis.annotated_a
        annotated_b = self.analysis.annotated_b
        shared_a = annotated_a["overlap_status"].eq("shared")
        shared_b = annotated_b["overlap_status"].eq("shared")
        specific_a = ~shared_a
        specific_b = ~shared_b

        if self.analysis.relationships.empty:
            qualifying_relationships = 0
            physical_relationships = 0
        else:
            qualifying_relationships = int(
                self.analysis.relationships["qualifying_overlap"].astype(bool).sum()
            )
            physical_relationships = int(len(self.analysis.relationships))

        qc_a = self.analysis.qc_a
        qc_b = self.analysis.qc_b
        metrics = [
            ("input_rows", qc_a["input_rows"], qc_b["input_rows"], "-"),
            (
                "non_peak_rows_removed",
                qc_a["non_peak_rows_removed"],
                qc_b["non_peak_rows_removed"],
                "-",
            ),
            (
                "fdr_filtered_rows",
                qc_a["fdr_filtered_rows"],
                qc_b["fdr_filtered_rows"],
                "-",
            ),
            (
                "fdr_missing_rows",
                qc_a["fdr_missing_rows"],
                qc_b["fdr_missing_rows"],
                "-",
            ),
            ("called_peaks", len(annotated_a), len(annotated_b), "-"),
            ("shared_peaks", int(shared_a.sum()), int(shared_b.sum()), "-"),
            ("specific_peaks", int(specific_a.sum()), int(specific_b.sum()), "-"),
            (
                "shared_peak_fraction",
                float(shared_a.mean()) if len(shared_a) else 0.0,
                float(shared_b.mean()) if len(shared_b) else 0.0,
                "-",
            ),
            (
                "transcripts_with_peaks",
                self._unique_count(annotated_a, "transcripts"),
                self._unique_count(annotated_b, "transcripts"),
                "-",
            ),
            (
                "genes_with_peaks",
                self._unique_count(annotated_a, "gene_name"),
                self._unique_count(annotated_b, "gene_name"),
                "-",
            ),
            ("physical_overlap_relationships", "-", "-", physical_relationships),
            (
                "qualifying_overlap_relationships",
                "-",
                "-",
                qualifying_relationships,
            ),
            (
                "shared_clusters",
                "-",
                "-",
                int(len(self.analysis.shared_clusters)),
            ),
            ("min_overlap_codons", "-", "-", self.analysis.min_overlap),
            (
                "min_reciprocal_overlap",
                "-",
                "-",
                self.analysis.min_reciprocal,
            ),
            ("interval_mode", "-", "-", self.analysis.interval_mode),
            (
                "max_fdr",
                "-",
                "-",
                self.analysis.max_fdr
                if self.analysis.max_fdr is not None
                else "None",
            ),
        ]
        self.summary = pd.DataFrame(
            metrics,
            columns=[
                "metric",
                self.analysis.name_a,
                self.analysis.name_b,
                "comparison",
            ],
        )
        return self.summary

    def build_length_compare(self) -> pd.DataFrame:
        """Build pairwise length metrics for qualifying shared relationships."""
        columns = [
            "cluster_id",
            "transcripts",
            "peak_uid_a",
            "peak_uid_b",
            "core_length_a",
            "core_length_b",
            "core_length_difference_b_minus_a",
            "match_span_length_a",
            "match_span_length_b",
            "span_length_difference_b_minus_a",
            "overlap_length",
            "overlap_fraction_a",
            "overlap_fraction_b",
            "reciprocal_overlap",
            "jaccard",
            "center_shift_b_minus_a",
        ]
        if self.analysis.relationships.empty:
            self.length_compare = pd.DataFrame(columns=columns)
            return self.length_compare

        self.length_compare = self.analysis.relationships.loc[
            self.analysis.relationships["qualifying_overlap"].astype(bool),
            columns,
        ].copy()
        return self.length_compare

    def prepare_tables(self) -> None:
        """Build all derived summary tables before writing output files."""
        self.build_summary()
        self.build_length_summary()
        self.build_length_compare()

    def write_results(self) -> dict[str, str]:
        """Write all overlap, shared, specific, and length-comparison tables.

        Returns
        -------
        dict[str, str]
            Mapping from logical result names to output paths.
        """
        self.prepare_tables()

        output_path = Path(self.analysis.output_prefix)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        files = {
            "summary": self.analysis.output_prefix + ".summary.txt",
            "length_summary": self.analysis.output_prefix + ".length_summary.txt",
            "relationships": self.analysis.output_prefix + ".relationships.txt",
            "shared_clusters": self.analysis.output_prefix + ".shared.clusters.txt",
            "shared_peaks": self.analysis.output_prefix + ".shared.peaks.txt",
            "specific_a": self.analysis.output_prefix
            + "."
            + self.analysis.file_label_a
            + ".specific.peaks.txt",
            "specific_b": self.analysis.output_prefix
            + "."
            + self.analysis.file_label_b
            + ".specific.peaks.txt",
            "length_compare": self.analysis.output_prefix + ".length_compare.txt",
        }

        self.summary.to_csv(files["summary"], sep="\t", index=False)
        self.length_summary.to_csv(files["length_summary"], sep="\t", index=False)
        self._public_columns(self.analysis.relationships).to_csv(
            files["relationships"], sep="\t", index=False
        )
        self.analysis.shared_clusters.to_csv(
            files["shared_clusters"], sep="\t", index=False
        )
        self._public_columns(self.analysis.shared_peaks).to_csv(
            files["shared_peaks"], sep="\t", index=False
        )
        self._public_columns(self.analysis.specific_a).to_csv(
            files["specific_a"], sep="\t", index=False
        )
        self._public_columns(self.analysis.specific_b).to_csv(
            files["specific_b"], sep="\t", index=False
        )
        self.length_compare.to_csv(files["length_compare"], sep="\t", index=False)
        return files