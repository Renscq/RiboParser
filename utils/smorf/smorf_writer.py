# Author: Rensc
# date: 2026-05-21

"""
Output writers for smORF scanning.

Supported outputs:
1. genePredExt-like ORF annotation.
2. Full ORF information table.
3. ORF nucleotide FASTA.
4. ORF peptide FASTA.
"""

from typing import List
from .smorf_models import ORFRecord


class GenePredWriter:
    """
    Write ORF annotations in genePredExt-like format.
    """

    @staticmethod
    def exon_frames(rec: ORFRecord) -> List[int]:
        """
        Calculate exon frame values for genePredExt output.

        Parameters
        ----------
        rec : ORFRecord
            ORF record.

        Returns
        -------
        list
            Exon frame list in genomic block order.
        """

        blocks = list(zip(rec.exon_starts, rec.exon_ends))

        # Coding order depends on strand.
        if rec.strand == "+":
            ordered = sorted(blocks, key=lambda x: x[0])
        else:
            ordered = sorted(blocks, key=lambda x: x[0], reverse=True)

        frames_by_block = {}
        offset = 0

        for block in ordered:
            frames_by_block[block] = offset % 3
            offset += block[1] - block[0]

        # genePred stores exon frames in genomic block order.
        return [frames_by_block[b] for b in sorted(blocks, key=lambda x: x[0])]

    @staticmethod
    def write(path: str, records: List[ORFRecord]) -> None:
        """
        Write ORF records to genePredExt-like file.

        Parameters
        ----------
        path : str
            Output file path.
        records : list
            List of ORFRecord objects.
        """

        with open(path, "w") as out:
            for rec in records:
                exon_count = len(rec.exon_starts)
                exon_starts = ",".join(str(x) for x in rec.exon_starts) + ","
                exon_ends = ",".join(str(x) for x in rec.exon_ends) + ","
                exon_frames = GenePredWriter.exon_frames(rec)
                exon_frames_str = ",".join(str(x) for x in exon_frames) + ","

                fields = [
                    rec.orf_id,
                    rec.chrom,
                    rec.strand,
                    str(rec.genomic_start),
                    str(rec.genomic_end),
                    str(rec.genomic_start),
                    str(rec.genomic_end),
                    str(exon_count),
                    exon_starts,
                    exon_ends,
                    "0",
                    rec.gene_id,
                    "cmpl" if rec.completeness == "complete" else "incmpl",
                    "cmpl" if rec.completeness == "complete" else "incmpl",
                    exon_frames_str,
                ]

                out.write("\t".join(fields) + "\n")


class MessageWriter:
    """
    Write full ORF information table.
    """

    HEADER = [
        "orf_id",
        "gene_id",
        "transcript_id",
        "chrom",
        "strand",
        "source_strand",
        "category",
        "priority",
        "overlap_type",
        "frame",
        "tx_orf_start",
        "tx_orf_end",
        "genomic_start",
        "genomic_end",
        "start_codon",
        "stop_codon",
        "nt_length",
        "aa_length",
        "kozak_seq",
        "completeness",
        "exon_count",
        "exon_starts",
        "exon_ends",
    ]

    @staticmethod
    def write(path: str, records: List[ORFRecord]) -> None:
        """
        Write ORF summary information to a tab-delimited file.

        Parameters
        ----------
        path : str
            Output file path.
        records : list
            List of ORFRecord objects.
        """

        with open(path, "w") as out:
            out.write("\t".join(MessageWriter.HEADER) + "\n")

            for rec in records:
                fields = [
                    rec.orf_id,
                    rec.gene_id,
                    rec.transcript_id,
                    rec.chrom,
                    rec.strand,
                    rec.source_strand,
                    rec.category,
                    rec.priority,
                    rec.overlap_type,
                    str(rec.frame),
                    str(rec.tx_orf_start),
                    str(rec.tx_orf_end),
                    str(rec.genomic_start),
                    str(rec.genomic_end),
                    rec.start_codon,
                    rec.stop_codon,
                    str(rec.nt_length),
                    str(rec.aa_length),
                    rec.kozak_seq,
                    rec.completeness,
                    str(len(rec.exon_starts)),
                    ",".join(str(x) for x in rec.exon_starts),
                    ",".join(str(x) for x in rec.exon_ends),
                ]

                out.write("\t".join(fields) + "\n")


class FastaWriter:
    """
    Write ORF nucleotide and peptide sequences in FASTA format.
    """

    @staticmethod
    def wrap(seq: str, width: int = 60) -> str:
        """
        Wrap FASTA sequence.

        Parameters
        ----------
        seq : str
            Input sequence.
        width : int
            Line width.

        Returns
        -------
        str
            Wrapped sequence.
        """

        return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

    @staticmethod
    def write_nt(path: str, records: List[ORFRecord]) -> None:
        """
        Write ORF nucleotide sequences.

        Parameters
        ----------
        path : str
            Output FASTA file.
        records : list
            List of ORFRecord objects.
        """

        with open(path, "w") as out:
            for rec in records:
                header = (
                    f">{rec.orf_id} gene={rec.gene_id} transcript={rec.transcript_id} "
                    f"type={rec.category} strand={rec.strand} length={rec.nt_length}"
                )

                out.write(header + "\n")
                out.write(FastaWriter.wrap(rec.nt_seq) + "\n")

    @staticmethod
    def write_pep(path: str, records: List[ORFRecord]) -> None:
        """
        Write ORF peptide sequences.

        Parameters
        ----------
        path : str
            Output FASTA file.
        records : list
            List of ORFRecord objects.
        """

        with open(path, "w") as out:
            for rec in records:
                header = (
                    f">{rec.orf_id} gene={rec.gene_id} transcript={rec.transcript_id} "
                    f"type={rec.category} strand={rec.strand} aa_length={rec.aa_length}"
                )

                out.write(header + "\n")
                out.write(FastaWriter.wrap(rec.pep_seq) + "\n")
