#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Rensc
date: 2026-05-23

@Project : smORF
@Script  : smorf_riboseq_integrate.py

Object-oriented integrator for long-format smORF Ribo-seq evidence tables.
"""

import gzip
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

import pandas as pd


@dataclass
class SmorfRiboSeqIntegrationConfig:
    """Store integration parameters."""
    input_file: str
    output_matrix: Optional[str] = None
    output_integrated: Optional[str] = None
    capture_labels: List[str] = field(
        default_factory=lambda: ["LowConfidence", "MediumConfidence", "HighConfidence"]
    )
    pass_labels: List[str] = field(
        default_factory=lambda: ["MediumConfidence", "HighConfidence"]
    )
    excellent_min_samples: int = 2


class SmorfRiboSeqEvidenceIntegrator:
    """Integrate sample-level smORF Ribo-seq evidence into downstream analysis tables."""

    STATIC_COLUMNS = [
        "orf_id",
        "gene_id",
        "transcript_id",
        "chrom",
        "strand",
        "category",
        "genomic_start",
        "genomic_end",
        "nt_length",
        "coding_nt_length",
        "coding_codon_count",
    ]

    SUM_COLUMNS = [
        "rpf_sum",
        "covered_nt",
        "covered_codon",
    ]

    MAX_COLUMNS = [
        "rpf_mean",
        "coverage_ratio",
        "covered_codon_ratio",
        "max_density",
        "frame0_density",
        "frame1_density",
        "frame2_density",
        "frame0_ratio",
        "frame1_ratio",
        "frame2_ratio",
        "start_pause_mean",
        "body_mean",
        "pre_stop_mean",
        "start_pause_ratio",
        "stop_pause_ratio",
        "post_stop_mean",
        "release_ratio",
        "release_drop_score",
    ]

    MIN_COLUMNS = [
        "codon_gini",
        "max_to_mean_ratio",
        "top10_fraction",
    ]

    LABEL_RANKS = {
        "translation_evidence": {
            "NoEvidence": 0,
            "LowConfidence": 1,
            "MediumConfidence": 2,
            "HighConfidence": 3,
        },
        "periodicity_label": {
            "Weak": 0,
            "Moderate": 1,
            "Strong": 2,
        },
        "pausing_label": {
            "Weak": 0,
            "Moderate": 1,
            "Strong": 2,
        },
        "release_label": {
            "Weak": 0,
            "Moderate": 1,
            "Strong": 2,
        },
        "coverage_shape": {
            "Skewed": 0,
            "Disperse": 1,
            "Uniform": 2,
        },
    }

    MATRIX_METRICS = [
        "rpf_sum",
        "frame0_density",
        "frame1_density",
        "frame2_density",
        "frame0_ratio",
        "frame1_ratio",
        "frame2_ratio",
        "coverage_ratio",
        "translation_evidence",
    ]

    def __init__(self, config: SmorfRiboSeqIntegrationConfig):
        """Initialize the evidence integrator."""
        self.config = config
        self.table: Optional[pd.DataFrame] = None
        self.matrix_table: Optional[pd.DataFrame] = None
        self.integrated_table: Optional[pd.DataFrame] = None

    @staticmethod
    def eprint(message: str) -> None:
        """Print progress information to stderr."""
        print(message, file=sys.stderr, flush=True)

    @staticmethod
    def smart_open(path: str, mode: str = "rt"):
        """Open plain or gzip-compressed files."""
        if path.endswith(".gz"):
            return gzip.open(path, mode)
        return open(path, mode)

    @staticmethod
    def split_items(value: Optional[str]) -> List[str]:
        """Split comma-separated items."""
        if value is None:
            return []

        text = str(value).strip()

        if text == "" or text.lower() in {"none", "na", "nan", "."}:
            return []

        return [item.strip() for item in text.split(",") if item.strip()]

    @staticmethod
    def first_non_missing(values: Iterable) -> str:
        """Return the first non-empty value."""
        for value in values:
            if pd.isna(value):
                continue

            text = str(value).strip()

            if text == "" or text.lower() in {"nan", "none", "na"}:
                continue

            return text

        return "."

    @staticmethod
    def unique_join(values: Iterable, sep: str = ",") -> str:
        """Join unique non-empty values while preserving order."""
        items = []
        seen = set()

        for value in values:
            if pd.isna(value):
                continue

            text = str(value).strip()

            if text == "" or text.lower() in {"nan", "none", "na", "."}:
                continue

            if text not in seen:
                seen.add(text)
                items.append(text)

        return sep.join(items) if items else "."

    @staticmethod
    def best_label(values: Iterable, rank_map: Dict[str, int]) -> str:
        """Select the highest-ranked label."""
        best = "NA"
        best_rank = -1

        for value in values:
            label = str(value).strip() if not pd.isna(value) else "NA"
            rank = rank_map.get(label, -1)

            if rank > best_rank:
                best = label
                best_rank = rank

        return best

    def read(self) -> "SmorfRiboSeqEvidenceIntegrator":
        """Read the input evidence table."""
        self.eprint(f"[Info] Reading evidence table: {self.config.input_file}")
        self.table = pd.read_csv(
            self.config.input_file,
            sep="\t",
            dtype=str,
            low_memory=False,
        )
        self._ensure_required_columns()
        self._coerce_numeric_columns()
        self._normalize_labels()

        return self

    def run(self) -> "SmorfRiboSeqEvidenceIntegrator":
        """Run all requested integration steps."""
        if self.table is None:
            self.read()

        if self.config.output_matrix:
            self.eprint(f"[Info] Building frame density matrix: {self.config.output_matrix}")
            self.matrix_table = self.build_frame_density_matrix()
            self.write_table(self.matrix_table, self.config.output_matrix)

        if self.config.output_integrated:
            self.eprint(f"[Info] Building integrated ORF evidence table: {self.config.output_integrated}")
            self.integrated_table = self.build_integrated_evidence_table()
            self.write_table(self.integrated_table, self.config.output_integrated)

        self.eprint("[Info] Done.")
        return self

    def _ensure_required_columns(self) -> None:
        """Check required columns."""
        required = {"sample", "orf_id"}
        missing = required - set(self.table.columns)

        if missing:
            raise ValueError(f"Evidence table is missing required columns: {missing}")

    def _coerce_numeric_columns(self) -> None:
        """Convert known metric columns to numeric values."""
        numeric_columns = set(self.SUM_COLUMNS + self.MAX_COLUMNS + self.MIN_COLUMNS)

        for column in numeric_columns:
            if column in self.table.columns:
                self.table[column] = pd.to_numeric(self.table[column], errors="coerce")

        for column in ["genomic_start", "genomic_end", "nt_length", "coding_nt_length", "coding_codon_count"]:
            if column in self.table.columns:
                self.table[column] = pd.to_numeric(self.table[column], errors="coerce")

    def _normalize_labels(self) -> None:
        """Normalize missing label columns."""
        for column in self.LABEL_RANKS:
            if column not in self.table.columns:
                self.table[column] = "NA"

            self.table[column] = self.table[column].fillna("NA").astype(str)

    def _get_static_table(self) -> pd.DataFrame:
        """Extract one static annotation row per ORF."""
        columns = [column for column in self.STATIC_COLUMNS if column in self.table.columns]

        if "orf_id" not in columns:
            columns = ["orf_id"] + columns

        static = self.table[columns].copy()
        static = static.drop_duplicates(subset=["orf_id"], keep="first")

        return static

    def build_frame_density_matrix(self) -> pd.DataFrame:
        """Build an ORF-by-sample matrix for RPF and frame-specific density."""
        static = self._get_static_table()
        result = static

        metrics = [column for column in self.MATRIX_METRICS if column in self.table.columns]

        for metric in metrics:
            wide = self._pivot_metric(metric)
            result = result.merge(wide, on="orf_id", how="left")

        result = result.fillna(0)

        return result

    def _pivot_metric(self, metric: str) -> pd.DataFrame:
        """Pivot one sample-specific metric into wide format."""
        sub = self.table[["orf_id", "sample", metric]].copy()

        if metric != "translation_evidence":
            sub[metric] = pd.to_numeric(sub[metric], errors="coerce").fillna(0)

        wide = sub.pivot_table(
            index="orf_id",
            columns="sample",
            values=metric,
            aggfunc="max" if metric != "translation_evidence" else "first",
            fill_value=0 if metric != "translation_evidence" else "",
        )

        wide.columns = [f"{sample}__{metric}" for sample in wide.columns]
        wide = wide.reset_index()

        return wide

    def build_integrated_evidence_table(self) -> pd.DataFrame:
        """Build one ORF-level integrated evidence table."""
        work = self.table.copy()

        if "rpf_sum" not in work.columns:
            work["rpf_sum"] = 0

        work["__is_captured"] = (
            work["translation_evidence"].isin(set(self.config.capture_labels))
            & (work["rpf_sum"].fillna(0) > 0)
        )
        work["__is_pass"] = work["translation_evidence"].isin(set(self.config.pass_labels))

        rows = []
        static_columns = [
            column for column in self.STATIC_COLUMNS
            if column in work.columns and column != "orf_id"
        ]

        for orf_id, group in work.groupby("orf_id", sort=False):
            row = self._summarize_one_orf(orf_id, group, static_columns)
            rows.append(row)

        result = pd.DataFrame(rows)

        best_sample = self._get_best_sample_table(work)
        result = result.merge(best_sample, on="orf_id", how="left")

        return result

    def _summarize_one_orf(self, orf_id: str, group: pd.DataFrame, static_columns: List[str]) -> dict:
        """Summarize all sample rows for one ORF."""
        row = {"orf_id": orf_id}

        for column in static_columns:
            row[column] = self.first_non_missing(group[column])

        captured = group[group["__is_captured"]].copy()
        passed = group[group["__is_pass"]].copy()

        row["samples"] = self.unique_join(captured["sample"]) if not captured.empty else "."
        row["sample_count"] = int(captured["sample"].nunique()) if not captured.empty else 0
        row["pass_samples"] = self.unique_join(passed["sample"]) if not passed.empty else "."
        row["pass_sample_count"] = int(passed["sample"].nunique()) if not passed.empty else 0

        self._add_sum_metrics(row, group)
        self._add_max_metrics(row, group)
        self._add_min_metrics(row, group)
        self._add_best_labels(row, group)
        self._add_multi_sample_status(row)

        return row

    def _add_sum_metrics(self, row: dict, group: pd.DataFrame) -> None:
        """Add summed metrics."""
        for column in self.SUM_COLUMNS:
            if column in group.columns:
                row[f"{column}_sum"] = float(group[column].fillna(0).sum())

    def _add_max_metrics(self, row: dict, group: pd.DataFrame) -> None:
        """Add max metrics."""
        for column in self.MAX_COLUMNS:
            if column in group.columns:
                row[f"{column}_max"] = float(group[column].fillna(0).max())

    def _add_min_metrics(self, row: dict, group: pd.DataFrame) -> None:
        """Add min metrics."""
        for column in self.MIN_COLUMNS:
            if column in group.columns:
                valid = group[column].dropna()
                row[f"{column}_min"] = float(valid.min()) if len(valid) > 0 else 0.0

    def _add_best_labels(self, row: dict, group: pd.DataFrame) -> None:
        """Add best labels."""
        for label_column, rank_map in self.LABEL_RANKS.items():
            if label_column in group.columns:
                row[f"{label_column}_best"] = self.best_label(group[label_column], rank_map)

    def _add_multi_sample_status(self, row: dict) -> None:
        """Add multi-sample evidence status."""
        if row["pass_sample_count"] >= self.config.excellent_min_samples:
            row["multi_sample_status"] = "Excellent"
        elif row["pass_sample_count"] == 1:
            row["multi_sample_status"] = "SingleSample"
        elif row["sample_count"] > 0:
            row["multi_sample_status"] = "CapturedButWeak"
        else:
            row["multi_sample_status"] = "NoEvidence"

    def _get_best_sample_table(self, table: pd.DataFrame) -> pd.DataFrame:
        """Get the best sample row per ORF by evidence rank and RPF sum."""
        work = table.copy()
        work["__evidence_rank"] = self._label_to_rank(work["translation_evidence"], "translation_evidence")
        work["__rpf_sum"] = pd.to_numeric(work["rpf_sum"], errors="coerce").fillna(0)

        work = work.sort_values(
            by=["orf_id", "__evidence_rank", "__rpf_sum"],
            ascending=[True, False, False],
        )

        best = work.drop_duplicates(subset=["orf_id"], keep="first").copy()

        keep_columns = ["orf_id", "sample", "translation_evidence"]
        rename_map = {
            "sample": "best_sample",
            "translation_evidence": "best_translation_evidence",
        }

        for column in ["rpf_sum", "coverage_ratio", "frame0_ratio", "release_ratio", "coverage_shape"]:
            if column in best.columns:
                keep_columns.append(column)
                rename_map[column] = f"best_{column}"

        return best[keep_columns].rename(columns=rename_map)

    def _label_to_rank(self, series: pd.Series, label_name: str) -> pd.Series:
        """Convert label series to rank values."""
        rank_map = self.LABEL_RANKS.get(label_name, {})
        return series.fillna("NA").astype(str).map(lambda value: rank_map.get(value, -1))

    @staticmethod
    def write_table(table: pd.DataFrame, path: str) -> None:
        """Write a TSV table."""
        if path.endswith(".gz"):
            table.to_csv(path, sep="\t", index=False, compression="gzip")
        else:
            table.to_csv(path, sep="\t", index=False)


def run_integration(args) -> None:
    """Run integration workflow from CLI arguments."""
    config = SmorfRiboSeqIntegrationConfig(
        input_file=args.input,
        output_matrix=args.output_matrix,
        output_integrated=args.output_integrated,
        capture_labels=SmorfRiboSeqEvidenceIntegrator.split_items(args.capture_labels)
        or ["LowConfidence", "MediumConfidence", "HighConfidence"],
        pass_labels=SmorfRiboSeqEvidenceIntegrator.split_items(args.pass_labels)
        or ["MediumConfidence", "HighConfidence"],
        excellent_min_samples=args.excellent_min_samples,
    )

    SmorfRiboSeqEvidenceIntegrator(config).run()
