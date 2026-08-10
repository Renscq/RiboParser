#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev001
# Function: Calculate validated CDS, codon-usage, CAI, amino-acid, and protein properties.
# Input: Coding-sequence FASTA and optional CAI reference FASTA.
# Output: Sequence QC, codon-usage, amino-acid composition, CAI, and protein-property tables.

"""Sequence-property algorithms for selective ribosome profiling analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex, GC123
from Bio.SeqUtils.ProtParam import ProteinAnalysis


AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
POSITIVE_AA = {"K", "R"}
NEGATIVE_AA = {"D", "E"}
AROMATIC_AA = {"F", "W", "Y"}


@dataclass(frozen=True)
class ValidatedCDS:
    """Store one validated coding sequence and its translated protein."""

    sequence_id: str
    input_sequence: str
    coding_sequence: str
    protein_sequence: str
    start_codon: str
    terminal_stop: str


@dataclass(frozen=True)
class SequenceValidation:
    """Store sequence QC fields and an optional validated CDS object."""

    qc_row: dict[str, object]
    validated: ValidatedCDS | None


def _three_letter_aa(amino_acid: str) -> str:
    """Return the conventional three-letter amino-acid abbreviation."""
    return IUPACData.protein_letters_1to3.get(amino_acid, amino_acid)


def build_codon_annotation(genetic_code: int) -> pd.DataFrame:
    """Build codon annotation from one NCBI genetic code table.

    Parameters
    ----------
    genetic_code : int
        NCBI genetic code identifier.

    Returns
    -------
    pandas.DataFrame
        Codon annotation with sense and stop codons.
    """
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    rows = []
    for codon in sorted(table.forward_table):
        amino_acid = table.forward_table[codon]
        rows.append(
            {
                "Codon": codon,
                "AA": _three_letter_aa(amino_acid),
                "Abbr.": amino_acid,
                "CodonType": "sense",
            }
        )
    for codon in sorted(table.stop_codons):
        rows.append(
            {
                "Codon": codon,
                "AA": "Stop",
                "Abbr.": "*",
                "CodonType": "stop",
            }
        )
    return pd.DataFrame(rows)


def _translate_coding_sequence(
    coding_sequence: str,
    genetic_code: int,
    start_codon: str,
) -> str:
    """Translate a stop-trimmed CDS and normalize recognized start codons to Met."""
    if not coding_sequence:
        return ""
    protein = str(Seq(coding_sequence).translate(table=genetic_code))
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    if protein and start_codon in table.start_codons:
        protein = "M" + protein[1:]
    return protein


def validate_cds_sequence(
    sequence_id: str,
    sequence: str,
    genetic_code: int,
) -> SequenceValidation:
    """Validate one CDS and separate the terminal stop from sense codons.

    Parameters
    ----------
    sequence_id : str
        FASTA record identifier.
    sequence : str
        Input nucleotide sequence.
    genetic_code : int
        NCBI genetic code identifier.

    Returns
    -------
    SequenceValidation
        QC information plus a validated CDS when no fatal problem is found.
    """
    input_sequence = str(sequence).upper()
    raw_sequence = input_sequence.replace("U", "T")
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    stop_codons = set(table.stop_codons)
    valid_codons = set(table.forward_table).union(stop_codons)

    issues = []
    fatal_issues = []
    if "U" in input_sequence:
        issues.append("u_converted_to_t")
    input_length = len(raw_sequence)
    start_codon = raw_sequence[:3] if input_length >= 3 else ""

    if input_length == 0:
        fatal_issues.append("empty_sequence")
    if input_length % 3 != 0:
        fatal_issues.append("length_not_multiple_of_3")

    codons = []
    if input_length > 0 and input_length % 3 == 0:
        codons = [raw_sequence[index:index + 3] for index in range(0, input_length, 3)]

    ambiguous_codons = [codon for codon in codons if codon not in valid_codons]
    if ambiguous_codons:
        fatal_issues.append("ambiguous_or_invalid_codon")

    terminal_stop = ""
    internal_stop_count = 0
    if codons and not ambiguous_codons:
        if codons[-1] in stop_codons:
            terminal_stop = codons[-1]
            sense_codons = codons[:-1]
        else:
            sense_codons = codons
            issues.append("missing_terminal_stop")

        internal_stop_count = sum(codon in stop_codons for codon in sense_codons)
        if internal_stop_count > 0:
            fatal_issues.append("internal_stop_codon")
        if not sense_codons:
            fatal_issues.append("no_sense_codon")
    else:
        sense_codons = []

    if start_codon:
        if start_codon == "ATG":
            start_status = "canonical"
        elif start_codon in table.start_codons:
            start_status = "alternative"
            issues.append("non_atg_start")
        else:
            start_status = "unrecognized"
            issues.append("non_atg_start")
    else:
        start_status = "missing"

    status = "FAIL" if fatal_issues else ("WARN" if issues else "PASS")
    all_issues = fatal_issues + issues
    coding_sequence = "".join(sense_codons) if not fatal_issues else ""
    protein_sequence = ""
    validated = None
    if not fatal_issues:
        protein_sequence = _translate_coding_sequence(
            coding_sequence=coding_sequence,
            genetic_code=genetic_code,
            start_codon=start_codon,
        )
        validated = ValidatedCDS(
            sequence_id=str(sequence_id),
            input_sequence=raw_sequence,
            coding_sequence=coding_sequence,
            protein_sequence=protein_sequence,
            start_codon=start_codon,
            terminal_stop=terminal_stop,
        )

    qc_row = {
        "ID": str(sequence_id),
        "Status": status,
        "Issues": ";".join(all_issues) if all_issues else "-",
        "Input_Length_nt": input_length,
        "Coding_Length_nt": len(coding_sequence),
        "Sense_Codon_Number": len(sense_codons) if not fatal_issues else 0,
        "Protein_Length_aa": len(protein_sequence),
        "Start_Codon": start_codon if start_codon else "-",
        "Start_Status": start_status,
        "Terminal_Stop": terminal_stop if terminal_stop else "-",
        "Internal_Stop_Number": internal_stop_count,
        "Invalid_Codon_Number": len(ambiguous_codons),
    }
    return SequenceValidation(qc_row=qc_row, validated=validated)


def count_sense_codons(sequence: str, sense_codons: list[str]) -> np.ndarray:
    """Count ordered sense codons in one stop-trimmed CDS."""
    codon_to_index = {codon: index for index, codon in enumerate(sense_codons)}
    counts = np.zeros(len(sense_codons), dtype=np.int32)
    for index in range(0, len(sequence), 3):
        codon = sequence[index:index + 3]
        counts[codon_to_index[codon]] += 1
    return counts


def calculate_rscu_matrix(
    codon_counts: np.ndarray,
    sense_codons: list[str],
    codon_to_aa: dict[str, str],
) -> np.ndarray:
    """Calculate RSCU for each sequence from sense-codon counts.

    RSCU is zero when an amino acid is absent from one sequence. For amino acids
    that are present, synonymous RSCU values sum to the number of synonymous
    codons.
    """
    counts = np.asarray(codon_counts, dtype=float)
    one_dimensional = counts.ndim == 1
    if one_dimensional:
        counts = counts.reshape(1, -1)

    rscu = np.zeros_like(counts, dtype=float)
    amino_acids = sorted(set(codon_to_aa.values()))
    for amino_acid in amino_acids:
        indices = [
            index
            for index, codon in enumerate(sense_codons)
            if codon_to_aa[codon] == amino_acid
        ]
        group_counts = counts[:, indices]
        group_total = group_counts.sum(axis=1)
        synonymous_number = float(len(indices))
        denominator = group_total / synonymous_number
        valid = denominator > 0
        if np.any(valid):
            rscu[np.ix_(valid, indices)] = (
                group_counts[valid, :] / denominator[valid, np.newaxis]
            )

    return rscu[0] if one_dimensional else rscu


def calculate_gc_properties(sequence: str) -> dict[str, float]:
    """Calculate GC, GC1, GC2, and GC3 from a stop-trimmed CDS."""
    if not sequence:
        return {"GC": np.nan, "GC1": np.nan, "GC2": np.nan, "GC3": np.nan}
    gc, gc1, gc2, gc3 = GC123(sequence)
    return {
        "GC": float(gc) / 100.0,
        "GC1": float(gc1) / 100.0,
        "GC2": float(gc2) / 100.0,
        "GC3": float(gc3) / 100.0,
    }


def calculate_protein_properties(
    protein_sequence: str,
    charge_ph: float,
) -> dict[str, float]:
    """Calculate physicochemical properties for one translated protein."""
    length = len(protein_sequence)
    if length == 0:
        return {
            "Gravy": np.nan,
            "Aromaticity": np.nan,
            "Flexibility": np.nan,
            "Instability": np.nan,
            "Isoelectric_Point": np.nan,
            "Net_Charge_pH": np.nan,
            "Charge_Density": np.nan,
            "Positive_Fraction": np.nan,
            "Negative_Fraction": np.nan,
            "Aromatic_Fraction": np.nan,
            "Proline_Fraction": np.nan,
            "Glycine_Fraction": np.nan,
            "Helix_Propensity_Fraction": np.nan,
            "Turn_Propensity_Fraction": np.nan,
            "Sheet_Propensity_Fraction": np.nan,
        }

    analysis = ProteinAnalysis(protein_sequence)
    flexibility_values = analysis.flexibility() if length >= 9 else []
    helix, turn, sheet = analysis.secondary_structure_fraction()
    net_charge = float(analysis.charge_at_pH(charge_ph))

    counts = {amino_acid: protein_sequence.count(amino_acid) for amino_acid in AA_ORDER}
    positive_count = sum(counts[amino_acid] for amino_acid in POSITIVE_AA)
    negative_count = sum(counts[amino_acid] for amino_acid in NEGATIVE_AA)
    aromatic_count = sum(counts[amino_acid] for amino_acid in AROMATIC_AA)

    return {
        "Gravy": float(analysis.gravy()),
        "Aromaticity": float(analysis.aromaticity()),
        "Flexibility": (
            float(np.mean(flexibility_values)) if flexibility_values else np.nan
        ),
        "Instability": float(analysis.instability_index()),
        "Isoelectric_Point": float(analysis.isoelectric_point()),
        "Net_Charge_pH": net_charge,
        "Charge_Density": net_charge / float(length),
        "Positive_Fraction": positive_count / float(length),
        "Negative_Fraction": negative_count / float(length),
        "Aromatic_Fraction": aromatic_count / float(length),
        "Proline_Fraction": counts["P"] / float(length),
        "Glycine_Fraction": counts["G"] / float(length),
        "Helix_Propensity_Fraction": float(helix),
        "Turn_Propensity_Fraction": float(turn),
        "Sheet_Propensity_Fraction": float(sheet),
    }


def calculate_amino_acid_composition(protein_sequence: str) -> dict[str, tuple[int, float]]:
    """Calculate counts and fractions for the 20 standard amino acids."""
    length = len(protein_sequence)
    result = {}
    for amino_acid in AA_ORDER:
        count = protein_sequence.count(amino_acid)
        fraction = count / float(length) if length > 0 else np.nan
        result[amino_acid] = (count, fraction)
    return result


class SeRPProperties:
    """Calculate validated coding-sequence and protein properties.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed ``serp_properties`` command-line arguments.
    """

    def __init__(self, args) -> None:
        """Initialize analysis settings and output paths."""
        self.fasta = Path(args.fasta)
        self.output_prefix = str(args.output)
        self.genetic_code = int(args.genetic_code)
        self.charge_ph = float(args.charge_ph)
        self.cai_reference = Path(args.cai_reference) if args.cai_reference else None

        try:
            self.codon_table = CodonTable.unambiguous_dna_by_id[self.genetic_code]
        except KeyError as error:
            raise ValueError(
                "Unknown NCBI genetic code: {0}".format(self.genetic_code)
            ) from error

        self.sense_codons = sorted(self.codon_table.forward_table)
        self.codon_to_aa = dict(self.codon_table.forward_table)
        self.codon_annotation = build_codon_annotation(self.genetic_code)

        self.sequence_qc_file = self.output_prefix + ".SequenceQC.txt"
        self.sequence_summary_file = self.output_prefix + "_summary.txt"
        self.frequency_file = self.output_prefix + "_frequency.txt"
        self.rscu_file = self.output_prefix + "_rscu.txt"
        self.cai_file = self.output_prefix + "_cai.txt"
        self.whole_codon_file = self.output_prefix + "_whole_codon_usage.txt"
        self.aa_composition_file = self.output_prefix + "_aa_composition.txt"
        self.properties_file = self.output_prefix + ".Properties.txt"
        self.cai_reference_qc_file = self.output_prefix + "_cai_reference_qc.txt"

        self.valid_sequences: list[ValidatedCDS] = []
        self.sequence_qc = pd.DataFrame()
        self.codon_counts = np.zeros((0, len(self.sense_codons)), dtype=np.int32)
        self.codon_frequency = np.zeros((0, len(self.sense_codons)), dtype=float)
        self.rscu = np.zeros((0, len(self.sense_codons)), dtype=float)
        self.cai_index: CodonAdaptationIndex | None = None
        self.cai_values: dict[str, float] = {}
        self.protein_properties = pd.DataFrame()
        self.aa_composition = pd.DataFrame()

    @staticmethod
    def _read_fasta_records(path: Path) -> list:
        """Read FASTA records and reject duplicate identifiers."""
        records = []
        seen = set()
        duplicates = []
        for record in SeqIO.parse(str(path), "fasta"):
            sequence_id = str(record.id)
            if sequence_id in seen:
                duplicates.append(sequence_id)
            else:
                seen.add(sequence_id)
            records.append(record)
        if duplicates:
            examples = ", ".join(sorted(set(duplicates))[:10])
            raise ValueError(
                "Duplicate FASTA identifiers detected: {0}".format(examples)
            )
        if not records:
            raise ValueError("No FASTA records were found in: {0}".format(path))
        return records

    def read_sequences(self) -> None:
        """Read CDS FASTA, validate all records, and build the codon matrix."""
        records = self._read_fasta_records(self.fasta)
        qc_rows = []
        valid_sequences = []

        for record in records:
            validation = validate_cds_sequence(
                sequence_id=str(record.id),
                sequence=str(record.seq),
                genetic_code=self.genetic_code,
            )
            qc_rows.append(validation.qc_row)
            if validation.validated is not None:
                valid_sequences.append(validation.validated)

        self.sequence_qc = pd.DataFrame(qc_rows)
        self.valid_sequences = valid_sequences
        if not self.valid_sequences:
            raise ValueError("No valid CDS sequences remain after sequence QC.")

        self.codon_counts = np.vstack(
            [
                count_sense_codons(sequence.coding_sequence, self.sense_codons)
                for sequence in self.valid_sequences
            ]
        )
        sense_codon_number = self.codon_counts.sum(axis=1).astype(float)
        self.codon_frequency = np.divide(
            self.codon_counts.astype(float) * 1000.0,
            sense_codon_number[:, np.newaxis],
            out=np.zeros_like(self.codon_counts, dtype=float),
            where=sense_codon_number[:, np.newaxis] > 0,
        )
        self.rscu = calculate_rscu_matrix(
            codon_counts=self.codon_counts,
            sense_codons=self.sense_codons,
            codon_to_aa=self.codon_to_aa,
        )

    def build_cai_index(self) -> None:
        """Build the Sharp-Li CAI reference index when a reference FASTA is given."""
        self.cai_index = None
        self.cai_values = {sequence.sequence_id: np.nan for sequence in self.valid_sequences}
        if self.cai_reference is None:
            return

        records = self._read_fasta_records(self.cai_reference)
        qc_rows = []
        reference_sequences = []
        for record in records:
            validation = validate_cds_sequence(
                sequence_id=str(record.id),
                sequence=str(record.seq),
                genetic_code=self.genetic_code,
            )
            qc_rows.append(validation.qc_row)
            if validation.validated is not None:
                reference_sequences.append(validation.validated.coding_sequence)

        reference_qc = pd.DataFrame(qc_rows)
        reference_qc.to_csv(
            self.cai_reference_qc_file,
            sep="\t",
            index=False,
            na_rep="NA",
        )
        if not reference_sequences:
            raise ValueError("No valid CDS sequences remain in the CAI reference FASTA.")

        self.cai_index = CodonAdaptationIndex(
            reference_sequences,
            table=self.codon_table,
        )
        for sequence in self.valid_sequences:
            try:
                value = float(self.cai_index.calculate(sequence.coding_sequence))
            except (TypeError, ValueError, ZeroDivisionError):
                value = np.nan
            self.cai_values[sequence.sequence_id] = value

    def calculate_properties(self) -> None:
        """Calculate nucleotide, protein, and amino-acid composition properties."""
        property_rows = []
        aa_rows = []

        for sequence in self.valid_sequences:
            gc_properties = calculate_gc_properties(sequence.coding_sequence)
            protein_properties = calculate_protein_properties(
                protein_sequence=sequence.protein_sequence,
                charge_ph=self.charge_ph,
            )
            cai_value = self.cai_values.get(sequence.sequence_id, np.nan)

            row = {
                "ID": sequence.sequence_id,
                "Seq": sequence.protein_sequence,
                "Length": len(sequence.protein_sequence),
                "CDS_Length_nt": len(sequence.coding_sequence),
                "Sense_Codon_Number": len(sequence.coding_sequence) // 3,
                "Start_Codon": sequence.start_codon,
                "Terminal_Stop": sequence.terminal_stop if sequence.terminal_stop else "-",
                "GC": gc_properties["GC"],
                "GC1": gc_properties["GC1"],
                "GC2": gc_properties["GC2"],
                "GC3": gc_properties["GC3"],
                "CAI": cai_value,
                "Charge_pH": self.charge_ph,
                **protein_properties,
            }
            property_rows.append(row)

            aa_composition = calculate_amino_acid_composition(sequence.protein_sequence)
            for amino_acid in AA_ORDER:
                count, fraction = aa_composition[amino_acid]
                aa_rows.append(
                    {
                        "ID": sequence.sequence_id,
                        "AA": amino_acid,
                        "AA3": _three_letter_aa(amino_acid),
                        "Count": count,
                        "Fraction": fraction,
                    }
                )

        self.protein_properties = pd.DataFrame(property_rows)
        self.aa_composition = pd.DataFrame(aa_rows)

    def _wide_codon_table(self, matrix: np.ndarray) -> pd.DataFrame:
        """Convert one sequence-by-codon matrix to the historical wide schema."""
        table = pd.DataFrame(matrix, columns=self.sense_codons)
        table.insert(
            0,
            "Length",
            [len(sequence.coding_sequence) // 3 for sequence in self.valid_sequences],
        )
        table.insert(0, "Gene", [sequence.sequence_id for sequence in self.valid_sequences])
        return table

    def write_tables(self) -> None:
        """Write QC, codon-usage, CAI, amino-acid, and protein-property tables."""
        self.sequence_qc.to_csv(
            self.sequence_qc_file,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        frequency_table = self._wide_codon_table(self.codon_frequency).round(6)
        rscu_table = self._wide_codon_table(self.rscu).round(6)
        frequency_table.to_csv(self.frequency_file, sep="\t", index=False)
        rscu_table.to_csv(self.rscu_file, sep="\t", index=False)

        cai_table = pd.DataFrame(
            {
                "Gene": [sequence.sequence_id for sequence in self.valid_sequences],
                "Length": [
                    len(sequence.coding_sequence) // 3 for sequence in self.valid_sequences
                ],
                "CAI": [
                    self.cai_values.get(sequence.sequence_id, np.nan)
                    for sequence in self.valid_sequences
                ],
                "CAI_Reference": (
                    str(self.cai_reference) if self.cai_reference is not None else "-"
                ),
            }
        )
        cai_table["CAI"] = pd.to_numeric(cai_table["CAI"], errors="coerce").round(6)
        cai_table.to_csv(self.cai_file, sep="\t", index=False, na_rep="NA")

        whole_counts = self.codon_counts.sum(axis=0).astype(int)
        whole_total = int(whole_counts.sum())
        whole_frequency = (
            whole_counts.astype(float) / float(whole_total) * 1000.0
            if whole_total > 0
            else np.zeros(len(self.sense_codons), dtype=float)
        )
        whole_rscu = calculate_rscu_matrix(
            codon_counts=whole_counts,
            sense_codons=self.sense_codons,
            codon_to_aa=self.codon_to_aa,
        )
        whole_table = pd.DataFrame(
            {
                "Codon": self.sense_codons,
                "AA": [_three_letter_aa(self.codon_to_aa[codon]) for codon in self.sense_codons],
                "Abbr.": [self.codon_to_aa[codon] for codon in self.sense_codons],
                "Count": whole_counts,
                "Frequency": whole_frequency,
                "RSCU": whole_rscu,
                "RelativeAdaptiveness": [
                    float(self.cai_index[codon]) if self.cai_index is not None else np.nan
                    for codon in self.sense_codons
                ],
            }
        )
        numeric_whole_columns = ["Frequency", "RSCU", "RelativeAdaptiveness"]
        whole_table[numeric_whole_columns] = whole_table[numeric_whole_columns].round(6)
        whole_table.to_csv(
            self.whole_codon_file,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        aa_composition = self.aa_composition.copy()
        aa_composition["Fraction"] = pd.to_numeric(
            aa_composition["Fraction"], errors="coerce"
        ).round(6)
        aa_composition.to_csv(
            self.aa_composition_file,
            sep="\t",
            index=False,
            na_rep="NA",
        )
        protein_properties = self.protein_properties.copy()
        numeric_property_columns = protein_properties.select_dtypes(
            include=[np.number]
        ).columns.tolist()
        protein_properties[numeric_property_columns] = protein_properties[
            numeric_property_columns
        ].round(6)
        protein_properties.to_csv(
            self.properties_file,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        summary = self.summary_counts()
        pd.DataFrame(
            [{"metric": key, "value": value} for key, value in summary.items()]
        ).to_csv(self.sequence_summary_file, sep="\t", index=False)

    def summary_counts(self) -> dict[str, object]:
        """Return concise sequence-QC and analysis counts."""
        status_counts = self.sequence_qc["Status"].value_counts()
        return {
            "input_sequence_number": int(len(self.sequence_qc)),
            "valid_sequence_number": int(len(self.valid_sequences)),
            "pass_sequence_number": int(status_counts.get("PASS", 0)),
            "warning_sequence_number": int(status_counts.get("WARN", 0)),
            "failed_sequence_number": int(status_counts.get("FAIL", 0)),
            "genetic_code": self.genetic_code,
            "charge_ph": self.charge_ph,
            "cai_calculated": self.cai_reference is not None,
            "cai_reference": str(self.cai_reference) if self.cai_reference is not None else "-",
        }