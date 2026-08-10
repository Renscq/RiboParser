#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev001
# Function: Aggregate and compare SeRP sequence properties across shared and specific peak groups.
# Input: Validated full-sequence properties and max-site-centered peak-local properties.
# Output: Group codon/AA tables, position profiles, transcript-aware statistics, summaries, and figures.

"""Group aggregation and statistics for SeRP sequence-property analysis."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from Bio.SeqUtils import ProtParamData
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from utils.serp.PeakProperties import SeRPPeakProperties
from utils.serp.Properties import (
    AA_ORDER,
    AROMATIC_AA,
    NEGATIVE_AA,
    POSITIVE_AA,
    SeRPProperties,
    _three_letter_aa,
    calculate_rscu_matrix,
)


LOCAL_AA_PSEUDOCOUNT = 0.5


class SeRPPropertyGroups:
    """Aggregate and statistically compare SeRP shared/specific sequence properties."""

    FULL_PROPERTY_FEATURES = [
        "Length",
        "GC3",
        "CAI",
        "Gravy",
        "Aromaticity",
        "Isoelectric_Point",
        "Charge_Density",
        "Positive_Fraction",
        "Negative_Fraction",
        "Proline_Fraction",
        "Glycine_Fraction",
    ]
    LOCAL_PROPERTY_FEATURES = [
        "Local_Length_aa",
        "GC3",
        "Gravy",
        "Signed_Charge_Fraction",
        "Positive_Fraction",
        "Negative_Fraction",
        "Aromatic_Fraction",
        "Proline_Fraction",
        "Glycine_Fraction",
    ]

    def __init__(
        self,
        args,
        sequence_workflow: SeRPProperties,
        peak_workflow: SeRPPeakProperties,
    ) -> None:
        """Initialize group analysis from validated sequence and peak-local workflows."""
        self.output_prefix = str(args.output)
        self.specific_set = str(args.specific_set)
        self.min_group_size = int(args.min_group_size)
        self.sequence_workflow = sequence_workflow
        self.peak_workflow = peak_workflow
        self.sample_order = list(peak_workflow.sample_order)
        self.base_groups = list(peak_workflow.base_groups)
        self.significant_groups = list(peak_workflow.significant_groups)
        self.plot_groups = list(peak_workflow.plot_groups)
        self.sequence_lookup = {
            sequence.sequence_id: sequence
            for sequence in sequence_workflow.valid_sequences
        }
        self.sequence_index = {
            sequence.sequence_id: index
            for index, sequence in enumerate(sequence_workflow.valid_sequences)
        }

        self.group_membership = pd.DataFrame()
        self.group_full_properties = pd.DataFrame()
        self.group_codon_usage = pd.DataFrame()
        self.group_aa_composition = pd.DataFrame()
        self.local_aa_composition = pd.DataFrame()
        self.position_aa_profile = pd.DataFrame()
        self.position_property_profile = pd.DataFrame()
        self.property_statistics = pd.DataFrame()
        self.group_summary = pd.DataFrame()

        self.group_membership_file = self.output_prefix + "_group_membership.txt"
        self.group_full_properties_file = self.output_prefix + "_group_protein_properties.txt"
        self.group_codon_usage_file = self.output_prefix + "_group_codon_usage.txt"
        self.group_aa_composition_file = self.output_prefix + "_group_aa_composition.txt"
        self.local_aa_composition_file = self.output_prefix + "_local_aa_composition.txt"
        self.position_aa_profile_file = self.output_prefix + "_position_aa_profile.txt"
        self.position_property_profile_file = self.output_prefix + "_position_property_profile.txt"
        self.statistics_file = self.output_prefix + "_property_statistics.txt"
        self.group_summary_file = self.output_prefix + "_group_summary.txt"

    def _build_group_membership(self) -> dict[str, set[str]]:
        """Build unique transcript sets for shared, specific, and significant-specific groups."""
        peaks = self.peak_workflow.base_peak_table
        membership: dict[str, set[str]] = {}
        membership["shared"] = set(
            peaks.loc[peaks["Analysis_Group"].eq("shared"), "transcripts"].astype(str)
        )
        for sample in self.sample_order:
            base_group = "{0}_specific".format(sample)
            significant_group = "{0}_significant_specific".format(sample)
            one = peaks.loc[peaks["Analysis_Group"].eq(base_group)]
            membership[base_group] = set(one["transcripts"].astype(str))
            membership[significant_group] = set(
                one.loc[one["is_significant_specific"].astype(bool), "transcripts"].astype(str)
            )
        return membership

    def calculate_group_properties(self) -> None:
        """Calculate full-protein, pooled codon, and pooled amino-acid properties by group."""
        membership = self._build_group_membership()
        membership_rows = []
        full_rows = []
        codon_rows = []
        aa_rows = []
        property_index = self.sequence_workflow.protein_properties.set_index("ID")

        for group, transcripts in membership.items():
            available = sorted(
                transcript for transcript in transcripts if transcript in self.sequence_index
            )
            missing = sorted(transcripts.difference(available))
            for transcript in available:
                membership_rows.append({"Group": group, "ID": transcript, "Status": "included"})
                full_rows.append(
                    {
                        "Group": group,
                        "ID": transcript,
                        **property_index.loc[transcript].to_dict(),
                    }
                )
            for transcript in missing:
                membership_rows.append(
                    {"Group": group, "ID": transcript, "Status": "absent_from_valid_cds"}
                )

            counts = np.zeros(len(self.sequence_workflow.sense_codons), dtype=int)
            if available:
                indices = [self.sequence_index[transcript] for transcript in available]
                counts = self.sequence_workflow.codon_counts[indices, :].sum(axis=0).astype(int)
            total = int(counts.sum())
            frequency = (
                counts.astype(float) / float(total) * 1000.0
                if total > 0
                else np.zeros(len(counts), dtype=float)
            )
            rscu = calculate_rscu_matrix(
                codon_counts=counts,
                sense_codons=self.sequence_workflow.sense_codons,
                codon_to_aa=self.sequence_workflow.codon_to_aa,
            )
            for codon, count, freq, rscu_value in zip(
                self.sequence_workflow.sense_codons,
                counts,
                frequency,
                rscu,
            ):
                aa = self.sequence_workflow.codon_to_aa[codon]
                codon_rows.append(
                    {
                        "Group": group,
                        "Codon": codon,
                        "AA": _three_letter_aa(aa),
                        "Abbr.": aa,
                        "Count": int(count),
                        "Frequency": float(freq),
                        "RSCU": float(rscu_value),
                    }
                )

            aa_counts = {aa: 0 for aa in AA_ORDER}
            for transcript in available:
                protein = self.sequence_lookup[transcript].protein_sequence
                for aa in AA_ORDER:
                    aa_counts[aa] += protein.count(aa)
            aa_total = sum(aa_counts.values())
            for aa in AA_ORDER:
                aa_rows.append(
                    {
                        "Group": group,
                        "AA": aa,
                        "AA3": _three_letter_aa(aa),
                        "Count": int(aa_counts[aa]),
                        "Fraction": aa_counts[aa] / float(aa_total) if aa_total > 0 else np.nan,
                    }
                )

        self.group_membership = pd.DataFrame(membership_rows)
        self.group_full_properties = pd.DataFrame(full_rows)
        self.group_codon_usage = pd.DataFrame(codon_rows)
        self.group_aa_composition = pd.DataFrame(aa_rows)

    def calculate_local_composition(self) -> None:
        """Calculate local AA composition and position-resolved property profiles."""
        local = self.peak_workflow.local_table_with_significant_groups()
        all_groups = self.base_groups + self.significant_groups
        local_rows = []
        position_count_rows = []
        background_count_rows = []
        property_position_rows = []

        for group in all_groups:
            one = local.loc[local["Plot_Group"].eq(group)]
            aa_counts = {aa: 0 for aa in AA_ORDER}
            for sequence in one["Context_AA"].astype(str):
                for aa in AA_ORDER:
                    aa_counts[aa] += sequence.count(aa)
            total = sum(aa_counts.values())
            for aa in AA_ORDER:
                local_rows.append(
                    {
                        "Group": group,
                        "AA": aa,
                        "AA3": _three_letter_aa(aa),
                        "Count": aa_counts[aa],
                        "Fraction": aa_counts[aa] / float(total) if total > 0 else np.nan,
                    }
                )

        for group in all_groups:
            one = local.loc[local["Plot_Group"].eq(group)]
            for _, row in one.iterrows():
                sequence = str(row["Context_AA"])
                relative_start = int(row["Relative_Start"])
                for offset, aa in enumerate(sequence):
                    relative_position = relative_start + offset
                    position_count_rows.append(
                        {
                            "Group": group,
                            "transcripts": row["transcripts"],
                            "Peak_ID": row["Peak_ID"],
                            "Relative_Position": relative_position,
                            "AA": aa,
                        }
                    )
                    property_position_rows.append(
                        {
                            "Group": group,
                            "transcripts": row["transcripts"],
                            "Peak_ID": row["Peak_ID"],
                            "Relative_Position": relative_position,
                            "Hydropathy": float(ProtParamData.kd[aa]),
                            "Signed_Charge_Class": (
                                1.0 if aa in POSITIVE_AA else (-1.0 if aa in NEGATIVE_AA else 0.0)
                            ),
                            "Positive": float(aa in POSITIVE_AA),
                            "Negative": float(aa in NEGATIVE_AA),
                            "Aromatic": float(aa in AROMATIC_AA),
                        }
                    )

        # Build the positional background from unique base peaks only. Significant-
        # specific groups are subsets and must not be counted twice in the background.
        for _, row in self.peak_workflow.local_peak_properties.iterrows():
            relative_start = int(row["Relative_Start"])
            for offset, aa in enumerate(str(row["Context_AA"])):
                background_count_rows.append(
                    {"Relative_Position": relative_start + offset, "AA": aa}
                )

        self.local_aa_composition = pd.DataFrame(local_rows)
        self._calculate_position_aa_profile(
            pd.DataFrame(position_count_rows),
            pd.DataFrame(background_count_rows),
            all_groups,
        )
        position_properties = pd.DataFrame(property_position_rows)
        if position_properties.empty:
            self.position_property_profile = pd.DataFrame()
        else:
            self.position_property_profile = (
                position_properties.groupby(["Group", "Relative_Position"], sort=False)
                .agg(
                    Sequence_Number=("Peak_ID", "count"),
                    Mean_Hydropathy=("Hydropathy", "mean"),
                    Mean_Signed_Charge_Class=("Signed_Charge_Class", "mean"),
                    Positive_Fraction=("Positive", "mean"),
                    Negative_Fraction=("Negative", "mean"),
                    Aromatic_Fraction=("Aromatic", "mean"),
                )
                .reset_index()
            )

    def _calculate_position_aa_profile(
        self,
        position_counts: pd.DataFrame,
        background_counts: pd.DataFrame,
        all_groups: list[str],
    ) -> None:
        """Calculate complete position-by-amino-acid frequency and enrichment grids."""
        if position_counts.empty or background_counts.empty:
            self.position_aa_profile = pd.DataFrame()
            return
        if self.peak_workflow.require_complete_local:
            positions = list(
                range(-self.peak_workflow.local_up, self.peak_workflow.local_down + 1)
            )
        else:
            positions = list(
                range(
                    int(position_counts["Relative_Position"].min()),
                    int(position_counts["Relative_Position"].max()) + 1,
                )
            )

        grouped = (
            position_counts.groupby(["Group", "Relative_Position", "AA"], sort=False)
            .size()
            .rename("Count")
        )
        complete_index = pd.MultiIndex.from_product(
            [all_groups, positions, AA_ORDER],
            names=["Group", "Relative_Position", "AA"],
        )
        counts = grouped.reindex(complete_index, fill_value=0).reset_index()
        counts["Position_Total"] = counts.groupby(
            ["Group", "Relative_Position"], sort=False
        )["Count"].transform("sum")
        counts["Fraction"] = np.where(
            counts["Position_Total"] > 0,
            counts["Count"] / counts["Position_Total"],
            np.nan,
        )

        background_grouped = (
            background_counts.groupby(["Relative_Position", "AA"], sort=False)
            .size()
            .rename("Background_Count")
        )
        background_index = pd.MultiIndex.from_product(
            [positions, AA_ORDER],
            names=["Relative_Position", "AA"],
        )
        background = background_grouped.reindex(background_index, fill_value=0).reset_index()
        background["Background_Total"] = background.groupby(
            "Relative_Position", sort=False
        )["Background_Count"].transform("sum")
        counts = counts.merge(
            background,
            on=["Relative_Position", "AA"],
            how="left",
            validate="many_to_one",
        )
        group_fraction = (counts["Count"] + LOCAL_AA_PSEUDOCOUNT) / (
            counts["Position_Total"] + LOCAL_AA_PSEUDOCOUNT * len(AA_ORDER)
        )
        background_fraction = (counts["Background_Count"] + LOCAL_AA_PSEUDOCOUNT) / (
            counts["Background_Total"] + LOCAL_AA_PSEUDOCOUNT * len(AA_ORDER)
        )
        valid = (counts["Position_Total"] > 0) & (counts["Background_Total"] > 0)
        counts["Log2_Enrichment_vs_All"] = np.where(
            valid,
            np.log2(group_fraction / background_fraction),
            np.nan,
        )
        self.position_aa_profile = counts

    def selected_full_properties(self) -> pd.DataFrame:
        """Return full-protein properties for the current statistics/figure groups."""
        if self.group_full_properties.empty:
            return self.group_full_properties.copy()
        return self.group_full_properties.loc[
            self.group_full_properties["Group"].isin(self.plot_groups)
        ].copy()

    def selected_local_properties(self) -> pd.DataFrame:
        """Return local properties for the current statistics/figure groups."""
        local = self.peak_workflow.local_table_with_significant_groups()
        return local.loc[local["Plot_Group"].isin(self.plot_groups)].copy()

    @staticmethod
    def _rank_biserial_from_u(u_statistic: float, n_a: int, n_b: int) -> float:
        """Calculate rank-biserial effect size with positive values meaning A > B."""
        return 2.0 * float(u_statistic) / float(n_a * n_b) - 1.0

    def _calculate_scope_statistics(
        self,
        table: pd.DataFrame,
        group_column: str,
        id_column: str,
        features: list[str],
        scope: str,
        aggregate_by_id: bool,
    ) -> list[dict[str, object]]:
        """Calculate pairwise tests after removing transcripts present in both groups."""
        rows = []
        existing_groups = set(table[group_column].astype(str))
        groups = [group for group in self.plot_groups if group in existing_groups]
        for group_a, group_b in combinations(groups, 2):
            one_a = table.loc[table[group_column].eq(group_a)].copy()
            one_b = table.loc[table[group_column].eq(group_b)].copy()
            if aggregate_by_id:
                numeric_features = [feature for feature in features if feature in table.columns]
                one_a = one_a.groupby(id_column, as_index=False)[numeric_features].median(
                    numeric_only=True
                )
                one_b = one_b.groupby(id_column, as_index=False)[numeric_features].median(
                    numeric_only=True
                )

            overlap_ids = set(one_a[id_column].astype(str)).intersection(
                set(one_b[id_column].astype(str))
            )
            if overlap_ids:
                one_a = one_a.loc[~one_a[id_column].astype(str).isin(overlap_ids)]
                one_b = one_b.loc[~one_b[id_column].astype(str).isin(overlap_ids)]

            for feature in features:
                if feature not in one_a.columns or feature not in one_b.columns:
                    continue
                values_a = pd.to_numeric(one_a[feature], errors="coerce").dropna()
                values_b = pd.to_numeric(one_b[feature], errors="coerce").dropna()
                row = {
                    "Scope": scope,
                    "Feature": feature,
                    "Group_A": group_a,
                    "Group_B": group_b,
                    "Overlap_Transcript_Number_Removed": len(overlap_ids),
                    "N_A": int(len(values_a)),
                    "N_B": int(len(values_b)),
                    "Median_A": float(values_a.median()) if len(values_a) else np.nan,
                    "Median_B": float(values_b.median()) if len(values_b) else np.nan,
                    "Median_Difference_A_minus_B": (
                        float(values_a.median() - values_b.median())
                        if len(values_a) and len(values_b)
                        else np.nan
                    ),
                    "Rank_Biserial": np.nan,
                    "P_Value": np.nan,
                    "BHFDR": np.nan,
                    "Test_Status": "insufficient_sample_size",
                }
                if len(values_a) >= self.min_group_size and len(values_b) >= self.min_group_size:
                    result = mannwhitneyu(
                        values_a.to_numpy(dtype=float),
                        values_b.to_numpy(dtype=float),
                        alternative="two-sided",
                        method="auto",
                    )
                    row["Rank_Biserial"] = self._rank_biserial_from_u(
                        result.statistic,
                        len(values_a),
                        len(values_b),
                    )
                    row["P_Value"] = float(result.pvalue)
                    row["Test_Status"] = "tested"
                rows.append(row)
        return rows

    def calculate_statistics(self) -> None:
        """Calculate transcript-aware full-protein and local-context comparisons."""
        rows = []
        full = self.selected_full_properties()
        if not full.empty:
            rows.extend(
                self._calculate_scope_statistics(
                    full,
                    "Group",
                    "ID",
                    self.FULL_PROPERTY_FEATURES,
                    "full_protein",
                    False,
                )
            )
        local = self.selected_local_properties()
        if not local.empty:
            rows.extend(
                self._calculate_scope_statistics(
                    local,
                    "Plot_Group",
                    "transcripts",
                    self.LOCAL_PROPERTY_FEATURES,
                    "local_context",
                    True,
                )
            )

        statistics = pd.DataFrame(rows)
        if not statistics.empty:
            valid = statistics["P_Value"].notna()
            if valid.any():
                statistics.loc[valid, "BHFDR"] = multipletests(
                    statistics.loc[valid, "P_Value"].to_numpy(dtype=float),
                    method="fdr_bh",
                )[1]
        self.property_statistics = statistics

    def calculate_group_summary(self) -> None:
        """Summarize peak, transcript, and extracted-context counts by group."""
        membership = self._build_group_membership()
        local = self.peak_workflow.local_table_with_significant_groups()
        peaks = self.peak_workflow.base_peak_table
        rows = []
        for group in self.base_groups + self.significant_groups:
            if group == "shared":
                peak_count = int(peaks["Analysis_Group"].eq("shared").sum())
            elif group.endswith("_significant_specific"):
                sample = group[: -len("_significant_specific")]
                peak_count = int(
                    (
                        peaks["Analysis_Group"].eq("{0}_specific".format(sample))
                        & peaks["is_significant_specific"].astype(bool)
                    ).sum()
                )
            else:
                peak_count = int(peaks["Analysis_Group"].eq(group).sum())
            rows.append(
                {
                    "Group": group,
                    "Peak_Number": peak_count,
                    "Transcript_Number": len(membership.get(group, set())),
                    "Local_Context_Number": int(local["Plot_Group"].eq(group).sum()),
                    "Complete_Local_Window_Required": self.peak_workflow.require_complete_local,
                    "Local_Upstream_aa": self.peak_workflow.local_up,
                    "Local_Downstream_aa": self.peak_workflow.local_down,
                }
            )
        self.group_summary = pd.DataFrame(rows)

    def calculate_all(self) -> None:
        """Run group aggregation, positional profiling, statistics, and summaries."""
        self.calculate_group_properties()
        self.calculate_local_composition()
        self.calculate_statistics()
        self.calculate_group_summary()

    @staticmethod
    def _round_numeric(table: pd.DataFrame, digits: int = 6) -> pd.DataFrame:
        """Round floating-point columns without altering integer count columns."""
        result = table.copy()
        float_columns = result.select_dtypes(include=["float", "float32", "float64"]).columns
        if len(float_columns):
            result[float_columns] = result[float_columns].round(digits)
        return result

    def write_tables(self) -> None:
        """Write all group aggregation, positional, statistical, and summary tables."""
        outputs = [
            (self.group_membership, self.group_membership_file),
            (self.group_full_properties, self.group_full_properties_file),
            (self.group_codon_usage, self.group_codon_usage_file),
            (self.group_aa_composition, self.group_aa_composition_file),
            (self.local_aa_composition, self.local_aa_composition_file),
            (self.position_aa_profile, self.position_aa_profile_file),
            (self.position_property_profile, self.position_property_profile_file),
            (self.property_statistics, self.statistics_file),
            (self.group_summary, self.group_summary_file),
        ]
        for table, path in outputs:
            self._round_numeric(table).to_csv(
                path,
                sep="\t",
                index=False,
                na_rep="NA",
            )

    def result_counts(self) -> dict[str, object]:
        """Return concise group-analysis settings and counts."""
        return {
            "specific_set_for_statistics_and_figures": self.specific_set,
            "property_comparison_group_number": len(self.plot_groups),
        }

    def draw_figures(self, args) -> None:
        """Draw group-comparison figures through the dedicated plotting layer."""
        from utils.serp.PropertiesPlot import SeRPPropertiesPlot

        SeRPPropertiesPlot(
            analyzer=self,
            output_prefix=self.output_prefix,
            output_format=args.output_format,
            dpi=args.dpi,
            font_size=args.font_size,
        ).draw_all()