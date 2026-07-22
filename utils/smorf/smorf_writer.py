#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Write smORF annotation, metadata, nucleotide, and peptide outputs.
# Input: Predicted ORFRecord objects.
# Output: genePredExt-like, tabular, nucleotide FASTA, and peptide FASTA files.

"""Output writers for transcript-centric smORF scanning.

``ORFOutputWriter`` supports streaming output and writes to temporary files.
Files are atomically moved to their final paths only after the complete scan
finishes successfully. This prevents partially written final outputs when a
worker fails and avoids retaining millions of ORF objects in memory.
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import TextIO

from .smorf_models import ORFRecord


class GenePredWriter:
    """Format and write genePredExt-like ORF annotations."""

    @staticmethod
    def exon_frames(record: ORFRecord) -> list[int]:
        """Calculate genePredExt exon-frame values.

        Args:
            record: ORF record containing genomic blocks.

        Returns:
            Exon frame values in ascending genomic block order.
        """
        blocks = list(zip(record.exon_starts, record.exon_ends))
        if record.strand == "+":
            coding_order = sorted(blocks, key=lambda block: block[0])
        else:
            coding_order = sorted(
                blocks,
                key=lambda block: block[0],
                reverse=True,
            )

        frame_by_block: dict[tuple[int, int], int] = {}
        coding_offset = 0
        for block in coding_order:
            frame_by_block[block] = coding_offset % 3
            coding_offset += block[1] - block[0]

        return [
            frame_by_block[block]
            for block in sorted(blocks, key=lambda block: block[0])
        ]

    @staticmethod
    def format_record(record: ORFRecord) -> str:
        """Format one ORF as a genePredExt-like line.

        Args:
            record: ORF record.

        Returns:
            Tab-delimited output line without a trailing newline.
        """
        exon_count = len(record.exon_starts)
        exon_starts = ",".join(
            str(value) for value in record.exon_starts
        ) + ","
        exon_ends = ",".join(
            str(value) for value in record.exon_ends
        ) + ","
        exon_frames = ",".join(
            str(value)
            for value in GenePredWriter.exon_frames(record)
        ) + ","
        completion = (
            "cmpl" if record.completeness == "complete" else "incmpl"
        )

        fields = [
            record.orf_id,
            record.chrom,
            record.strand,
            str(record.genomic_start),
            str(record.genomic_end),
            str(record.genomic_start),
            str(record.genomic_end),
            str(exon_count),
            exon_starts,
            exon_ends,
            "0",
            record.gene_id,
            completion,
            completion,
            exon_frames,
        ]
        return "\t".join(fields)

    @staticmethod
    def write(path: str | Path, records: list[ORFRecord]) -> None:
        """Write all ORFs to a genePredExt-like file.

        Args:
            path: Output path.
            records: ORF records.
        """
        with Path(path).open("w", encoding="utf-8") as handle:
            for record in records:
                handle.write(GenePredWriter.format_record(record) + "\n")


class MessageWriter:
    """Format and write the full ORF metadata table."""

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
        "kozak_start_index",
        "ambiguous_codon_count",
    ]

    @staticmethod
    def format_record(record: ORFRecord) -> str:
        """Format one ORF metadata row.

        Args:
            record: ORF record.

        Returns:
            Tab-delimited output line without a trailing newline.
        """
        fields = [
            record.orf_id,
            record.gene_id,
            record.transcript_id,
            record.chrom,
            record.strand,
            record.source_strand,
            record.category,
            record.priority,
            record.overlap_type,
            str(record.frame),
            str(record.tx_orf_start),
            str(record.tx_orf_end),
            str(record.genomic_start),
            str(record.genomic_end),
            record.start_codon,
            record.stop_codon,
            str(record.nt_length),
            str(record.aa_length),
            record.kozak_seq,
            record.completeness,
            str(len(record.exon_starts)),
            ",".join(str(value) for value in record.exon_starts),
            ",".join(str(value) for value in record.exon_ends),
            str(record.kozak_start_index),
            str(record.ambiguous_codon_count),
        ]
        return "\t".join(fields)

    @staticmethod
    def write(path: str | Path, records: list[ORFRecord]) -> None:
        """Write all ORF metadata rows.

        Args:
            path: Output path.
            records: ORF records.
        """
        with Path(path).open("w", encoding="utf-8") as handle:
            handle.write("\t".join(MessageWriter.HEADER) + "\n")
            for record in records:
                handle.write(MessageWriter.format_record(record) + "\n")


class FastaWriter:
    """Format and write nucleotide or peptide FASTA records."""

    @staticmethod
    def wrap(sequence: str, width: int = 60) -> str:
        """Wrap a FASTA sequence.

        Args:
            sequence: Sequence to wrap.
            width: Maximum line width.

        Returns:
            Wrapped sequence.

        Raises:
            ValueError: If ``width`` is not positive.
        """
        if width < 1:
            raise ValueError("FASTA line width must be >= 1.")
        return "\n".join(
            sequence[index:index + width]
            for index in range(0, len(sequence), width)
        )

    @staticmethod
    def format_nt(record: ORFRecord) -> str:
        """Format one nucleotide FASTA record.

        Args:
            record: ORF record.

        Returns:
            Complete FASTA record ending with a newline.
        """
        header = (
            f">{record.orf_id} gene={record.gene_id} "
            f"transcript={record.transcript_id} type={record.category} "
            f"strand={record.strand} length={record.nt_length}"
        )
        return f"{header}\n{FastaWriter.wrap(record.nt_seq)}\n"

    @staticmethod
    def format_pep(record: ORFRecord) -> str:
        """Format one peptide FASTA record.

        Args:
            record: ORF record.

        Returns:
            Complete FASTA record ending with a newline.
        """
        header = (
            f">{record.orf_id} gene={record.gene_id} "
            f"transcript={record.transcript_id} type={record.category} "
            f"strand={record.strand} aa_length={record.aa_length}"
        )
        return f"{header}\n{FastaWriter.wrap(record.pep_seq)}\n"

    @staticmethod
    def write_nt(path: str | Path, records: list[ORFRecord]) -> None:
        """Write nucleotide FASTA records.

        Args:
            path: Output path.
            records: ORF records.
        """
        with Path(path).open("w", encoding="utf-8") as handle:
            for record in records:
                handle.write(FastaWriter.format_nt(record))

    @staticmethod
    def write_pep(path: str | Path, records: list[ORFRecord]) -> None:
        """Write peptide FASTA records.

        Args:
            path: Output path.
            records: ORF records.
        """
        with Path(path).open("w", encoding="utf-8") as handle:
            for record in records:
                handle.write(FastaWriter.format_pep(record))


class ORFOutputWriter:
    """Atomically stream all four smORF output formats.

    Args:
        output_prefix: Prefix shared by final output files.
    """

    def __init__(self, output_prefix: str | Path) -> None:
        """Initialize final and temporary output paths.

        Args:
            output_prefix: Prefix shared by all output files.
        """
        prefix = Path(output_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)

        self.final_paths = {
            "genepred": Path(f"{prefix}.genePred"),
            "message": Path(f"{prefix}.message.txt"),
            "nt": Path(f"{prefix}.nt.fa"),
            "pep": Path(f"{prefix}.pep.fa"),
        }
        self.temporary_paths = {
            name: path.with_name(path.name + ".tmp")
            for name, path in self.final_paths.items()
        }
        self.handles: dict[str, TextIO] = {}
        self.record_count = 0

    def __enter__(self) -> "ORFOutputWriter":
        """Open temporary output files.

        Returns:
            The active output writer.
        """
        try:
            self.handles = {
                name: path.open("w", encoding="utf-8")
                for name, path in self.temporary_paths.items()
            }
            self.handles["message"].write(
                "\t".join(MessageWriter.HEADER) + "\n"
            )
        except Exception:
            for handle in self.handles.values():
                handle.close()
            self.handles.clear()
            for temporary_path in self.temporary_paths.values():
                temporary_path.unlink(missing_ok=True)
            raise
        return self

    def write_records(self, records: list[ORFRecord]) -> None:
        """Append one transcript's ORF records.

        Args:
            records: ORF records with final stable identifiers.

        Raises:
            RuntimeError: If the writer is not active.
        """
        if not self.handles:
            raise RuntimeError("ORFOutputWriter is not active.")

        for record in records:
            self.handles["genepred"].write(
                GenePredWriter.format_record(record) + "\n"
            )
            self.handles["message"].write(
                MessageWriter.format_record(record) + "\n"
            )
            self.handles["nt"].write(FastaWriter.format_nt(record))
            self.handles["pep"].write(FastaWriter.format_pep(record))
            self.record_count += 1

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        """Close files and commit or discard temporary outputs.

        Args:
            exception_type: Raised exception type, if any.
            exception: Raised exception, if any.
            traceback: Raised exception traceback, if any.

        Returns:
            ``False`` so exceptions propagate to the caller.
        """
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()

        if exception_type is None:
            for name, temporary_path in self.temporary_paths.items():
                temporary_path.replace(self.final_paths[name])
        else:
            for temporary_path in self.temporary_paths.values():
                temporary_path.unlink(missing_ok=True)

        return False
