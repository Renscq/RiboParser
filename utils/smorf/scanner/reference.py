#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Load reference sequences and map transcript/genomic coordinates.
# Input: Genome FASTA, genePred annotation, and transcript sequences.
# Output: Parsed transcripts, sequence utilities, and coordinate mappings.

"""Reference parsing, sequence utilities, and coordinate conversion."""

from __future__ import annotations

import gzip

from .models import Transcript

STOP_CODONS = {"TAA", "TAG", "TGA"}

GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


class SeqUtils:
    """
    Basic sequence operation utilities.
    """

    # Translation table for reverse-complement conversion.
    COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """
        Return the reverse-complement sequence.

        Parameters
        ----------
        seq : str
            Input nucleotide sequence.

        Returns
        -------
        str
            Reverse-complement sequence in uppercase.
        """

        return seq.translate(SeqUtils.COMP)[::-1].upper()

    @staticmethod
    def translate(seq: str) -> str:
        """
        Translate nucleotide sequence into peptide sequence.

        Unknown or ambiguous codons are translated as X.

        Parameters
        ----------
        seq : str
            Input nucleotide sequence.

        Returns
        -------
        str
            Translated peptide sequence.
        """

        aa = []
        seq = seq.upper()

        # Translate only complete codons.
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i : i + 3]
            aa.append(GENETIC_CODE.get(codon, "X"))

        return "".join(aa)


class FastaParser:
    """
    Read genome FASTA files into an in-memory dictionary.
    """

    @staticmethod
    def open_file(path: str):
        """
        Open a plain text or gzip-compressed FASTA file.

        Parameters
        ----------
        path : str
            Input FASTA file path.

        Returns
        -------
        file object
            Opened file handle.
        """

        if path.endswith(".gz"):
            return gzip.open(path, "rt")

        return open(path)

    @staticmethod
    def read_fasta(path: str) -> dict[str, str]:
        """
        Read FASTA file into a dictionary.

        The sequence ID is extracted from the first token after '>'.

        Parameters
        ----------
        path : str
            Input FASTA file path.

        Returns
        -------
        dict
            Dictionary of {sequence_id: sequence}.
        """

        seqs = {}
        name = None
        chunks = []

        with FastaParser.open_file(path) as handle:
            for line in handle:
                line = line.strip()

                if not line:
                    continue

                # Start a new FASTA record.
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(chunks).upper()

                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)

        # Save the final FASTA record.
        if name is not None:
            seqs[name] = "".join(chunks).upper()

        return seqs


class GenePredParser:
    """
    Parse genePred annotation into Transcript objects.
    """

    @staticmethod
    def parse_int_list(value: str) -> list[int]:
        """
        Parse comma-separated integer fields from genePred.

        Parameters
        ----------
        value : str
            Comma-separated integer string.

        Returns
        -------
        list
            List of integers.
        """

        value = value.rstrip(",")

        if not value:
            return []

        return [int(x) for x in value.split(",") if x != ""]

    @staticmethod
    def read_genepred(path: str) -> list[Transcript]:
        """
        Read genePred file and return a list of Transcript objects.

        Parameters
        ----------
        path : str
            Input genePred file.

        Returns
        -------
        list
            List of Transcript objects.
        """

        transcripts = []

        with open(path) as handle:
            for line in handle:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")

                if len(fields) < 10:
                    raise ValueError(f"Invalid genePred line: {line}")

                transcript_id = fields[0]
                chrom = fields[1]
                strand = fields[2]
                tx_start = int(fields[3])
                tx_end = int(fields[4])
                cds_start = int(fields[5])
                cds_end = int(fields[6])
                exon_count = int(fields[7])
                exon_starts = GenePredParser.parse_int_list(fields[8])
                exon_ends = GenePredParser.parse_int_list(fields[9])

                # Validate exon structure.
                if len(exon_starts) != exon_count or len(exon_ends) != exon_count:
                    raise ValueError(f"Exon count mismatch in transcript: {transcript_id}")

                # In genePredExt, column 12 is name2, usually gene ID.
                gene_id = fields[11] if len(fields) >= 12 else transcript_id

                transcripts.append(
                    Transcript(
                        transcript_id=transcript_id,
                        chrom=chrom,
                        strand=strand,
                        tx_start=tx_start,
                        tx_end=tx_end,
                        cds_start=cds_start,
                        cds_end=cds_end,
                        exon_starts=exon_starts,
                        exon_ends=exon_ends,
                        gene_id=gene_id,
                    )
                )

        return transcripts


class CoordinateMapper:
    """
    Coordinate mapper between genomic space and transcript space.
    """

    @staticmethod
    def build_transcript_sequence(tx: Transcript, genome: dict[str, str]) -> None:
        """
        Build spliced transcript sequence from genome FASTA.

        For minus-strand transcripts, exons are extracted in reverse order
        and reverse-complemented.

        Parameters
        ----------
        tx : Transcript
            Transcript object.
        genome : dict
            Genome sequence dictionary.
        """

        if tx.chrom not in genome:
            raise KeyError(f"Chromosome not found in genome FASTA: {tx.chrom}")

        chrom_seq = genome[tx.chrom]
        parts = []

        if tx.strand == "+":
            exon_iter = zip(tx.exon_starts, tx.exon_ends)

            for start, end in exon_iter:
                parts.append(chrom_seq[start:end])
        else:
            exon_iter = zip(reversed(tx.exon_starts), reversed(tx.exon_ends))

            for start, end in exon_iter:
                parts.append(SeqUtils.reverse_complement(chrom_seq[start:end]))

        tx.tx_seq = "".join(parts).upper()

        # Convert CDS genomic interval to transcript interval.
        tx.cds_tx_start, tx.cds_tx_end = CoordinateMapper.genomic_interval_to_tx_interval(
            tx,
            tx.cds_start,
            tx.cds_end,
        )

    @staticmethod
    def genomic_interval_to_tx_interval(
        tx: Transcript,
        g_start: int,
        g_end: int,
    ) -> tuple[int | None, int | None]:
        """
        Convert a genomic interval to transcript coordinates.

        Parameters
        ----------
        tx : Transcript
            Transcript object.
        g_start : int
            Genomic interval start.
        g_end : int
            Genomic interval end.

        Returns
        -------
        tuple
            Transcript start and end coordinates.
        """

        if g_end <= g_start:
            return None, None

        tx_ranges = []
        offset = 0

        if tx.strand == "+":
            exon_iter = list(zip(tx.exon_starts, tx.exon_ends))

            for e_start, e_end in exon_iter:
                ov_start = max(e_start, g_start)
                ov_end = min(e_end, g_end)

                if ov_start < ov_end:
                    t_start = offset + (ov_start - e_start)
                    t_end = offset + (ov_end - e_start)
                    tx_ranges.append((t_start, t_end))

                offset += e_end - e_start
        else:
            exon_iter = list(zip(reversed(tx.exon_starts), reversed(tx.exon_ends)))

            for e_start, e_end in exon_iter:
                ov_start = max(e_start, g_start)
                ov_end = min(e_end, g_end)

                if ov_start < ov_end:
                    t_start = offset + (e_end - ov_end)
                    t_end = offset + (e_end - ov_start)
                    tx_ranges.append((t_start, t_end))

                offset += e_end - e_start

        if not tx_ranges:
            return None, None

        return min(x[0] for x in tx_ranges), max(x[1] for x in tx_ranges)

    @staticmethod
    def tx_interval_to_genomic_blocks(
        tx: Transcript,
        t_start: int,
        t_end: int,
    ) -> tuple[list[int], list[int]]:
        """
        Convert a transcript interval to genomic exon blocks.

        This function is used to map ORFs back to genePred-compatible
        genomic exon block coordinates.

        Parameters
        ----------
        tx : Transcript
            Transcript object.
        t_start : int
            Transcript interval start.
        t_end : int
            Transcript interval end.

        Returns
        -------
        tuple
            Genomic exon start list and exon end list.
        """

        blocks = []
        offset = 0

        if tx.strand == "+":
            exon_iter = list(zip(tx.exon_starts, tx.exon_ends))

            for e_start, e_end in exon_iter:
                exon_len = e_end - e_start
                seg_start = offset
                seg_end = offset + exon_len

                ov_start = max(t_start, seg_start)
                ov_end = min(t_end, seg_end)

                if ov_start < ov_end:
                    g_start = e_start + (ov_start - seg_start)
                    g_end = e_start + (ov_end - seg_start)
                    blocks.append((g_start, g_end))

                offset += exon_len
        else:
            exon_iter = list(zip(reversed(tx.exon_starts), reversed(tx.exon_ends)))

            for e_start, e_end in exon_iter:
                exon_len = e_end - e_start
                seg_start = offset
                seg_end = offset + exon_len

                ov_start = max(t_start, seg_start)
                ov_end = min(t_end, seg_end)

                if ov_start < ov_end:
                    g_start = e_end - (ov_end - seg_start)
                    g_end = e_end - (ov_start - seg_start)
                    blocks.append((g_start, g_end))

                offset += exon_len

        # genePred requires genomic blocks sorted by genomic coordinate.
        blocks = sorted(blocks, key=lambda x: x[0])

        return [x[0] for x in blocks], [x[1] for x in blocks]
