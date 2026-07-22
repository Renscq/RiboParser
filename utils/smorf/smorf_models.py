#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Define data models used by transcript-centric smORF scanning.
# Input: Parsed transcript annotation and predicted ORF attributes.
# Output: Validated Transcript and ORFRecord objects.

"""Data models used by the transcript-centric smORF scanner.

Coordinate Conventions:
    Genomic coordinates use the UCSC genePred convention: zero-based,
    half-open intervals. Transcript coordinates are also zero-based and
    half-open, but always follow the 5-prime-to-3-prime direction of the
    reconstructed transcript sequence.

The models intentionally contain no file I/O. Validation that depends on a
particular input format belongs in the corresponding parser.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class Transcript:
    """Store one transcript parsed from a genePred record.

    Attributes:
        transcript_id: Unique transcript identifier.
        chrom: Chromosome or contig identifier.
        strand: Genomic strand, either ``+`` or ``-``.
        tx_start: Genomic transcript start, zero-based and inclusive.
        tx_end: Genomic transcript end, zero-based and exclusive.
        cds_start: Genomic CDS start, zero-based and inclusive.
        cds_end: Genomic CDS end, zero-based and exclusive.
        exon_starts: Genomic exon starts in ascending genomic order.
        exon_ends: Genomic exon ends in ascending genomic order.
        gene_id: Gene identifier associated with the transcript.
        tx_seq: Reconstructed spliced transcript sequence in transcript
            orientation.
        cds_tx_start: CDS start in transcript coordinates.
        cds_tx_end: CDS end in transcript coordinates.
    """

    transcript_id: str
    chrom: str
    strand: str
    tx_start: int
    tx_end: int
    cds_start: int
    cds_end: int
    exon_starts: list[int]
    exon_ends: list[int]
    gene_id: str = "."
    tx_seq: str = ""
    cds_tx_start: int | None = None
    cds_tx_end: int | None = None

    def is_coding(self) -> bool:
        """Return whether the transcript contains a non-empty CDS.

        Returns:
            ``True`` when ``cds_end`` is greater than ``cds_start``.
        """
        return self.cds_end > self.cds_start

    def exon_count(self) -> int:
        """Return the number of exon blocks.

        Returns:
            Number of exon intervals.
        """
        return len(self.exon_starts)

    def transcript_length(self) -> int:
        """Return the spliced transcript length.

        Returns:
            Sum of all exon lengths.
        """
        return sum(
            exon_end - exon_start
            for exon_start, exon_end in zip(
                self.exon_starts,
                self.exon_ends,
            )
        )


@dataclass(slots=True)
class ORFRecord:
    """Store one predicted ORF and all output-ready annotations.

    Attributes:
        orf_id: Stable ORF identifier assigned by the pipeline.
        transcript_id: Source transcript identifier.
        gene_id: Source gene identifier.
        chrom: Chromosome or contig identifier.
        strand: Genomic strand of the predicted ORF.
        source_strand: Scanning orientation relative to the transcript,
            either ``sense`` or ``antisense``.
        frame: Reading frame in the scanned sequence.
        tx_orf_start: ORF start in transcript coordinates.
        tx_orf_end: ORF end in transcript coordinates.
        genomic_start: Minimum genomic coordinate covered by the ORF.
        genomic_end: Maximum genomic coordinate covered by the ORF.
        exon_starts: ORF genomic block starts in ascending genomic order.
        exon_ends: ORF genomic block ends in ascending genomic order.
        start_codon: Detected start codon.
        stop_codon: Detected stop codon or ``NA`` for a partial ORF.
        nt_length: ORF nucleotide length, including a terminal stop codon
            when present.
        aa_length: Peptide length excluding the terminal stop symbol.
        nt_seq: ORF nucleotide sequence in coding orientation.
        pep_seq: Translated peptide sequence.
        kozak_seq: Fixed-width, ``N``-padded start-codon context.
        category: Positional ORF category assigned by the classifier.
        priority: Priority label used by overlap filtering.
        overlap_type: Relationship to overlapping ORFs.
        completeness: ``complete`` or ``3prime_partial``.
        kozak_start_index: Zero-based start-codon index in ``kozak_seq``.
        ambiguous_codon_count: Number of ORF codons containing non-ACGT
            characters.
        exon_frames: Per-exon frame values used by genePredExt output.
    """

    orf_id: str
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    source_strand: str
    frame: int
    tx_orf_start: int
    tx_orf_end: int
    genomic_start: int
    genomic_end: int
    exon_starts: list[int]
    exon_ends: list[int]
    start_codon: str
    stop_codon: str
    nt_length: int
    aa_length: int
    nt_seq: str
    pep_seq: str
    kozak_seq: str
    category: str = "unknown"
    priority: str = "primary"
    overlap_type: str = "none"
    completeness: str = "complete"
    kozak_start_index: int = 6
    ambiguous_codon_count: int = 0
    exon_frames: list[int] = field(default_factory=list)
