#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev001
# Function: Build shared/specific SeRP peak groups and extract max-site-centered local sequences.
# Input: Validated CDS records and serp_summary peak-classification table.
# Output: Local sequence QC, peak-local properties, and significant-specific peak properties.

"""Peak-group parsing and local sequence extraction for SeRP properties."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from utils.serp.Properties import (
    AA_ORDER,
    AROMATIC_AA,
    NEGATIVE_AA,
    POSITIVE_AA,
    SeRPProperties,
    calculate_gc_properties,
)


class SeRPPeakProperties:
    """Parse SeRP peak groups and extract fixed local sequence contexts."""

    REQUIRED_PEAK_COLUMNS = {
        "transcripts",
        "peak_category",
        "comparison_sample",
        "peak_start",
        "peak_end",
    }

    def __init__(self, args, sequence_workflow: SeRPProperties) -> None:
        """Initialize peak grouping and local-window settings."""
        self.peak_classification_file = Path(args.peak_classification)
        self.output_prefix = str(args.output)
        self.local_up = int(args.local_up)
        self.local_down = int(args.local_down)
        self.require_complete_local = bool(args.require_complete_local)
        self.specific_set = str(args.specific_set)
        self.sequence_workflow = sequence_workflow
        self.sequence_lookup = {
            sequence.sequence_id: sequence
            for sequence in sequence_workflow.valid_sequences
        }

        self.peak_table = pd.DataFrame()
        self.shared_representatives = pd.DataFrame()
        self.base_peak_table = pd.DataFrame()
        self.local_sequence_qc = pd.DataFrame()
        self.local_peak_properties = pd.DataFrame()
        self.significant_specific_properties = pd.DataFrame()

        self.sample_order: list[str] = []
        self.base_groups: list[str] = []
        self.significant_groups: list[str] = []
        self.plot_groups: list[str] = []

        self.local_qc_file = self.output_prefix + ".PeakSequenceQC.txt"
        self.local_properties_file = self.output_prefix + ".PeakProperties.txt"
        self.significant_specific_file = (
            self.output_prefix + ".SignificantSpecificPeakProperties.txt"
        )

    @staticmethod
    def _safe_numeric(series: pd.Series) -> pd.Series:
        """Convert one column to numeric values while accepting historical dash fields."""
        return pd.to_numeric(series.replace("-", pd.NA), errors="coerce")

    @staticmethod
    def _safe_bool(series: pd.Series) -> pd.Series:
        """Convert common text and numeric Boolean representations to bool."""
        return (
            series.fillna(False)
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(["true", "1", "yes", "y"])
        )

    @staticmethod
    def _require_columns(
        table: pd.DataFrame,
        required: set[str],
        label: str,
    ) -> None:
        """Validate required columns in one input table."""
        missing = sorted(required.difference(table.columns))
        if missing:
            raise ValueError(
                "{0} is missing required column(s): {1}".format(
                    label,
                    ", ".join(missing),
                )
            )

    @staticmethod
    def _peak_midpoint(row: pd.Series) -> float:
        """Return one peak midpoint as a fallback local-sequence anchor."""
        return (float(row["peak_start"]) + float(row["peak_end"])) / 2.0

    def read_peak_classification(self) -> None:
        """Read serp_summary peak classification and build shared/specific base groups."""
        table = pd.read_csv(
            self.peak_classification_file,
            sep="\t",
            low_memory=False,
        )
        self._require_columns(
            table,
            self.REQUIRED_PEAK_COLUMNS,
            str(self.peak_classification_file),
        )
        if table.empty:
            raise ValueError(
                "Peak-classification table contains no called peaks: {0}".format(
                    self.peak_classification_file
                )
            )

        for column in ["peak_start", "peak_end", "peak_num", "max_site"]:
            if column in table.columns:
                table[column] = self._safe_numeric(table[column])
        invalid_coordinate = (
            table["peak_start"].isna()
            | table["peak_end"].isna()
            | (table["peak_start"] > table["peak_end"])
        )
        if invalid_coordinate.any():
            examples = table.loc[
                invalid_coordinate,
                ["transcripts", "peak_start", "peak_end"],
            ].head(10)
            raise ValueError(
                "Invalid peak coordinates were found in peak classification. Example:\n{0}".format(
                    examples.to_string(index=False)
                )
            )

        table["peak_category"] = (
            table["peak_category"].fillna("").astype(str).str.strip().str.lower()
        )
        invalid_category = ~table["peak_category"].isin(["shared", "specific"])
        if invalid_category.any():
            values = sorted(table.loc[invalid_category, "peak_category"].unique())
            raise ValueError(
                "Unsupported peak_category value(s): {0}".format(", ".join(values))
            )

        table["comparison_sample"] = (
            table["comparison_sample"].fillna("").astype(str).str.strip()
        )
        self.sample_order = [
            value
            for value in table["comparison_sample"].drop_duplicates().tolist()
            if value
        ]
        if len(self.sample_order) != 2:
            raise ValueError(
                "Expected exactly two comparison samples in peak classification; found: {0}".format(
                    ", ".join(self.sample_order)
                )
            )

        if "is_significant_specific" in table.columns:
            table["is_significant_specific"] = self._safe_bool(
                table["is_significant_specific"]
            )
        else:
            table["is_significant_specific"] = False
        if "significance_status" not in table.columns:
            table["significance_status"] = "untested"
        else:
            table["significance_status"] = (
                table["significance_status"].fillna("untested").astype(str)
            )

        table["transcripts"] = table["transcripts"].astype(str)
        self.peak_table = table
        self.shared_representatives = self._build_shared_representatives(table)
        specific = self._build_specific_peaks(table)
        self.base_peak_table = pd.concat(
            [self.shared_representatives, specific],
            ignore_index=True,
            sort=False,
        )

        self.base_groups = ["shared"] + [
            "{0}_specific".format(sample) for sample in self.sample_order
        ]
        self.significant_groups = [
            "{0}_significant_specific".format(sample)
            for sample in self.sample_order
        ]
        self.plot_groups = (
            ["shared"] + self.significant_groups
            if self.specific_set == "significant"
            else list(self.base_groups)
        )

    def _build_shared_representatives(self, table: pd.DataFrame) -> pd.DataFrame:
        """Collapse A/B members of each shared cluster to one sequence-analysis anchor."""
        shared = table.loc[table["peak_category"].eq("shared")].copy()
        columns = [
            "Peak_ID",
            "transcripts",
            "gene_name",
            "Analysis_Group",
            "comparison_sample",
            "peak_category",
            "peak_start",
            "peak_end",
            "max_site",
            "Anchor_Source",
            "significance_status",
            "is_significant_specific",
            "Shared_Member_Number",
        ]
        if shared.empty:
            return pd.DataFrame(columns=columns)
        if "shared_cluster_id" not in shared.columns:
            raise ValueError(
                "Shared peaks require shared_cluster_id in the current serp_summary output."
            )
        cluster_id = shared["shared_cluster_id"].fillna("").astype(str).str.strip()
        if cluster_id.isin(["", "-", "nan", "None"]).any():
            raise ValueError("Some shared peaks are missing shared_cluster_id.")
        shared["shared_cluster_id"] = cluster_id

        rows = []
        for cluster, one in shared.groupby("shared_cluster_id", sort=False):
            transcripts = one["transcripts"].dropna().astype(str).unique().tolist()
            if len(transcripts) != 1:
                raise ValueError(
                    "Shared cluster {0} contains multiple transcripts: {1}".format(
                        cluster,
                        ", ".join(transcripts),
                    )
                )
            anchors = (
                self._safe_numeric(one["max_site"]).dropna()
                if "max_site" in one.columns
                else pd.Series(dtype=float)
            )
            if anchors.empty:
                anchors = one.apply(self._peak_midpoint, axis=1)
                anchor_source = "peak_midpoint_median"
            else:
                anchor_source = "max_site_median"
            anchor = int(np.rint(float(np.median(anchors))))

            statuses = one["significance_status"].astype(str).tolist()
            if statuses and all(value == "significant" for value in statuses):
                cluster_status = "all_significant"
            elif any(value == "significant" for value in statuses):
                cluster_status = "partly_significant"
            elif statuses and all(value == "untested" for value in statuses):
                cluster_status = "untested"
            else:
                cluster_status = "no_significant_peak"

            gene_name = "-"
            if "gene_name" in one.columns:
                valid_gene = one["gene_name"].fillna("-").astype(str).replace("", "-")
                valid_gene = valid_gene.loc[~valid_gene.isin(["-", "nan", "None"])]
                if not valid_gene.empty:
                    gene_name = valid_gene.iloc[0]

            rows.append(
                {
                    "Peak_ID": str(cluster),
                    "transcripts": transcripts[0],
                    "gene_name": gene_name,
                    "Analysis_Group": "shared",
                    "comparison_sample": "shared",
                    "peak_category": "shared",
                    "peak_start": int(one["peak_start"].min()),
                    "peak_end": int(one["peak_end"].max()),
                    "max_site": anchor,
                    "Anchor_Source": anchor_source,
                    "significance_status": cluster_status,
                    "is_significant_specific": False,
                    "Shared_Member_Number": int(len(one)),
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def _build_specific_peaks(self, table: pd.DataFrame) -> pd.DataFrame:
        """Build one local-analysis row per condition-specific called peak."""
        specific = table.loc[table["peak_category"].eq("specific")].copy()
        if specific.empty:
            return pd.DataFrame(columns=self.shared_representatives.columns)

        rows = []
        for index, row in specific.iterrows():
            sample = str(row["comparison_sample"])
            if not sample:
                raise ValueError("A specific peak is missing comparison_sample.")
            max_site = row.get("max_site", np.nan)
            if pd.isna(max_site):
                max_site = self._peak_midpoint(row)
                anchor_source = "peak_midpoint"
            else:
                anchor_source = "max_site"
            peak_id = row.get("peak_uid", "")
            if pd.isna(peak_id) or str(peak_id).strip() in ["", "-", "nan"]:
                peak_id = "specific_{0}".format(index + 1)
            rows.append(
                {
                    "Peak_ID": str(peak_id),
                    "transcripts": str(row["transcripts"]),
                    "gene_name": str(row.get("gene_name", "-")),
                    "Analysis_Group": "{0}_specific".format(sample),
                    "comparison_sample": sample,
                    "peak_category": "specific",
                    "peak_start": int(row["peak_start"]),
                    "peak_end": int(row["peak_end"]),
                    "max_site": int(np.rint(float(max_site))),
                    "Anchor_Source": anchor_source,
                    "significance_status": str(row.get("significance_status", "untested")),
                    "is_significant_specific": bool(row["is_significant_specific"]),
                    "Shared_Member_Number": 1,
                }
            )
        return pd.DataFrame(rows)

    @staticmethod
    def _calculate_local_properties(
        amino_acid_sequence: str,
        nucleotide_sequence: str,
    ) -> dict[str, float]:
        """Calculate properties appropriate for an internal local protein segment."""
        length = len(amino_acid_sequence)
        counts = {aa: amino_acid_sequence.count(aa) for aa in AA_ORDER}
        positive = sum(counts[aa] for aa in POSITIVE_AA)
        negative = sum(counts[aa] for aa in NEGATIVE_AA)
        aromatic = sum(counts[aa] for aa in AROMATIC_AA)
        gc = calculate_gc_properties(nucleotide_sequence)
        analysis = ProteinAnalysis(amino_acid_sequence)
        return {
            "Local_Length_aa": length,
            "GC": gc["GC"],
            "GC3": gc["GC3"],
            "Gravy": float(analysis.gravy()),
            "Signed_Charge_Fraction": (positive - negative) / float(length),
            "Positive_Fraction": positive / float(length),
            "Negative_Fraction": negative / float(length),
            "Aromatic_Fraction": aromatic / float(length),
            "Proline_Fraction": counts["P"] / float(length),
            "Glycine_Fraction": counts["G"] / float(length),
        }

    def extract_local_sequences(self) -> None:
        """Extract max-site-centered CDS/protein contexts and calculate local properties."""
        qc_rows = []
        property_rows = []
        expected_length = self.local_up + self.local_down + 1

        for _, peak in self.base_peak_table.iterrows():
            transcript = str(peak["transcripts"])
            sequence = self.sequence_lookup.get(transcript)
            anchor = int(peak["max_site"])
            qc = {
                "Peak_ID": peak["Peak_ID"],
                "transcripts": transcript,
                "gene_name": peak.get("gene_name", "-"),
                "Analysis_Group": peak["Analysis_Group"],
                "max_site": anchor,
                "Anchor_Source": peak["Anchor_Source"],
                "Status": "PASS",
                "Issue": "-",
                "Requested_Start": anchor - self.local_up,
                "Requested_End": anchor + self.local_down,
                "Extracted_Start": np.nan,
                "Extracted_End": np.nan,
                "Protein_Length_aa": np.nan,
            }
            if sequence is None:
                qc["Status"] = "SKIP"
                qc["Issue"] = "transcript_absent_from_valid_cds"
                qc_rows.append(qc)
                continue

            protein_length = len(sequence.protein_sequence)
            qc["Protein_Length_aa"] = protein_length
            if anchor < 0 or anchor >= protein_length:
                qc["Status"] = "SKIP"
                qc["Issue"] = "max_site_outside_cds"
                qc_rows.append(qc)
                continue

            requested_start = anchor - self.local_up
            requested_end = anchor + self.local_down
            if self.require_complete_local and (
                requested_start < 0 or requested_end >= protein_length
            ):
                qc["Status"] = "SKIP"
                qc["Issue"] = "incomplete_local_window"
                qc_rows.append(qc)
                continue

            start = max(0, requested_start)
            end = min(protein_length - 1, requested_end)
            aa_sequence = sequence.protein_sequence[start:end + 1]
            nt_sequence = sequence.coding_sequence[start * 3:(end + 1) * 3]
            qc["Extracted_Start"] = start
            qc["Extracted_End"] = end
            if start != requested_start or end != requested_end:
                qc["Status"] = "WARN"
                qc["Issue"] = "partial_local_window"
            qc_rows.append(qc)

            property_rows.append(
                {
                    **peak.to_dict(),
                    "Context_Start": start,
                    "Context_End": end,
                    "Relative_Start": start - anchor,
                    "Relative_End": end - anchor,
                    "Context_AA": aa_sequence,
                    "Context_NT": nt_sequence,
                    **self._calculate_local_properties(
                        amino_acid_sequence=aa_sequence,
                        nucleotide_sequence=nt_sequence,
                    ),
                }
            )

        self.local_sequence_qc = pd.DataFrame(qc_rows)
        self.local_peak_properties = pd.DataFrame(property_rows)
        if self.local_peak_properties.empty:
            raise ValueError(
                "No peak-local sequences remain after CDS matching and local-window QC."
            )
        self.significant_specific_properties = self.local_peak_properties.loc[
            self.local_peak_properties["is_significant_specific"].astype(bool)
        ].copy()

        if self.require_complete_local:
            lengths = self.local_peak_properties["Local_Length_aa"].astype(int)
            if not lengths.eq(expected_length).all():
                raise RuntimeError("Complete local windows have inconsistent lengths.")

    def local_table_with_significant_groups(self) -> pd.DataFrame:
        """Return base local peaks plus significant-specific subset group labels."""
        base = self.local_peak_properties.copy()
        base["Plot_Group"] = base["Analysis_Group"]
        extra_rows = []
        for sample in self.sample_order:
            base_group = "{0}_specific".format(sample)
            significant_group = "{0}_significant_specific".format(sample)
            one = base.loc[
                base["Analysis_Group"].eq(base_group)
                & base["is_significant_specific"].astype(bool)
            ].copy()
            one["Plot_Group"] = significant_group
            extra_rows.append(one)
        if extra_rows:
            return pd.concat([base] + extra_rows, ignore_index=True, sort=False)
        return base

    def write_tables(self) -> None:
        """Write peak-local QC and property tables."""
        for table, path in [
            (self.local_sequence_qc, self.local_qc_file),
            (self.local_peak_properties, self.local_properties_file),
            (self.significant_specific_properties, self.significant_specific_file),
        ]:
            result = table.copy()
            float_columns = result.select_dtypes(include=["float", "float32", "float64"]).columns
            if len(float_columns):
                result[float_columns] = result[float_columns].round(6)
            result.to_csv(path, sep="\t", index=False, na_rep="NA")

    def result_counts(self) -> dict[str, object]:
        """Return concise peak parsing and local-sequence counts."""
        skipped = (
            int(self.local_sequence_qc["Status"].eq("SKIP").sum())
            if not self.local_sequence_qc.empty
            else 0
        )
        return {
            "peak_classification_rows": int(len(self.peak_table)),
            "shared_cluster_representatives": int(len(self.shared_representatives)),
            "specific_peaks": int(
                self.base_peak_table["peak_category"].eq("specific").sum()
            ),
            "local_contexts": int(len(self.local_peak_properties)),
            "local_contexts_skipped": skipped,
        }