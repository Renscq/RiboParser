#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-24
# Version: 0.2.8.24-dev.002
# Function: Write smORF outputs with batched buffered I/O.
# Input: Predicted ORFRecord objects.
# Output: genePredExt-like, metadata, nucleotide FASTA, and peptide FASTA files.

"""Output writers for transcript-centric smORF scanning."""

from __future__ import annotations

from pathlib import Path
import shutil
from types import TracebackType
from typing import TextIO

from .smorf_models import ORFRecord


class GenePredWriter:
    """Format and write genePredExt-like ORF annotations."""

    @staticmethod
    def exon_frames(record: ORFRecord) -> list[int]:
        """Calculate genePredExt exon-frame values."""
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
        """Format one ORF as a genePredExt-like line."""
        exon_starts = ",".join(map(str, record.exon_starts)) + ","
        exon_ends = ",".join(map(str, record.exon_ends)) + ","
        exon_frames = ",".join(
            map(str, GenePredWriter.exon_frames(record))
        ) + ","
        completion = (
            "cmpl" if record.completeness == "complete" else "incmpl"
        )
        return "\t".join(
            [
                record.orf_id,
                record.chrom,
                record.strand,
                str(record.genomic_start),
                str(record.genomic_end),
                str(record.genomic_start),
                str(record.genomic_end),
                str(len(record.exon_starts)),
                exon_starts,
                exon_ends,
                "0",
                record.gene_id,
                completion,
                completion,
                exon_frames,
            ]
        )

    @staticmethod
    def write(path: str | Path, records: list[ORFRecord]) -> None:
        """Write all ORFs to a genePredExt-like file."""
        with Path(path).open("w", encoding="utf-8") as handle:
            handle.write(
                "".join(
                    GenePredWriter.format_record(record) + "\n"
                    for record in records
                )
            )


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
        """Format one ORF metadata row."""
        return "\t".join(
            [
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
                ",".join(map(str, record.exon_starts)),
                ",".join(map(str, record.exon_ends)),
                str(record.kozak_start_index),
                str(record.ambiguous_codon_count),
            ]
        )

    @staticmethod
    def write(path: str | Path, records: list[ORFRecord]) -> None:
        """Write all ORF metadata rows."""
        with Path(path).open("w", encoding="utf-8") as handle:
            handle.write("\t".join(MessageWriter.HEADER) + "\n")
            handle.write(
                "".join(
                    MessageWriter.format_record(record) + "\n"
                    for record in records
                )
            )


class FastaWriter:
    """Format and write nucleotide or peptide FASTA records."""

    @staticmethod
    def wrap(sequence: str, width: int = 60) -> str:
        """Wrap a FASTA sequence."""
        if width < 1:
            raise ValueError("FASTA line width must be >= 1.")
        return "\n".join(
            sequence[index:index + width]
            for index in range(0, len(sequence), width)
        )

    @staticmethod
    def format_nt(record: ORFRecord) -> str:
        """Format one nucleotide FASTA record."""
        header = (
            f">{record.orf_id} gene={record.gene_id} "
            f"transcript={record.transcript_id} type={record.category} "
            f"strand={record.strand} length={record.nt_length}"
        )
        return f"{header}\n{FastaWriter.wrap(record.nt_seq)}\n"

    @staticmethod
    def format_pep(record: ORFRecord) -> str:
        """Format one peptide FASTA record."""
        header = (
            f">{record.orf_id} gene={record.gene_id} "
            f"transcript={record.transcript_id} type={record.category} "
            f"strand={record.strand} aa_length={record.aa_length}"
        )
        return f"{header}\n{FastaWriter.wrap(record.pep_seq)}\n"

    @staticmethod
    def write_nt(path: str | Path, records: list[ORFRecord]) -> None:
        """Write nucleotide FASTA records."""
        with Path(path).open("w", encoding="utf-8") as handle:
            handle.write(
                "".join(FastaWriter.format_nt(record) for record in records)
            )

    @staticmethod
    def write_pep(path: str | Path, records: list[ORFRecord]) -> None:
        """Write peptide FASTA records."""
        with Path(path).open("w", encoding="utf-8") as handle:
            handle.write(
                "".join(FastaWriter.format_pep(record) for record in records)
            )


class ORFOutputWriter:
    """Atomically stream all four smORF output formats."""

    def __init__(
        self,
        output_prefix: str | Path,
        write_message_header: bool = True,
    ) -> None:
        """Initialize final and temporary output paths."""
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
        self.write_message_header = bool(write_message_header)

    def __enter__(self) -> "ORFOutputWriter":
        """Open temporary output files."""
        try:
            self.handles = {
                name: path.open("w", encoding="utf-8", buffering=1024 * 1024)
                for name, path in self.temporary_paths.items()
            }
            if self.write_message_header:
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
        """Append one transcript's ORF records using four batched writes."""
        if not self.handles:
            raise RuntimeError("ORFOutputWriter is not active.")
        if not records:
            return

        self.handles["genepred"].write(
            "".join(
                GenePredWriter.format_record(record) + "\n"
                for record in records
            )
        )
        self.handles["message"].write(
            "".join(
                MessageWriter.format_record(record) + "\n"
                for record in records
            )
        )
        self.handles["nt"].write(
            "".join(FastaWriter.format_nt(record) for record in records)
        )
        self.handles["pep"].write(
            "".join(FastaWriter.format_pep(record) for record in records)
        )
        self.record_count += len(records)

    @staticmethod
    def merge_parts(
        output_prefix: str | Path,
        part_prefixes: list[str | Path],
    ) -> None:
        """Merge ordered worker shards into atomic final output files.

        Parameters
        ----------
        output_prefix : str or pathlib.Path
            Final output prefix.
        part_prefixes : list of str or pathlib.Path
            Worker shard prefixes in transcript input order.
        """
        prefix = Path(output_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        suffixes = {
            "genepred": ".genePred",
            "message": ".message.txt",
            "nt": ".nt.fa",
            "pep": ".pep.fa",
        }
        final_paths = {
            name: Path(f"{prefix}{suffix}")
            for name, suffix in suffixes.items()
        }
        temporary_paths = {
            name: path.with_name(path.name + ".tmp")
            for name, path in final_paths.items()
        }

        handles: dict[str, object] = {}
        try:
            handles = {
                name: path.open("wb", buffering=1024 * 1024)
                for name, path in temporary_paths.items()
            }
            handles["message"].write(
                ("\t".join(MessageWriter.HEADER) + "\n").encode("utf-8")
            )

            for part_prefix in part_prefixes:
                part = Path(part_prefix)
                for name, suffix in suffixes.items():
                    part_path = Path(f"{part}{suffix}")
                    with part_path.open("rb", buffering=1024 * 1024) as source:
                        if name == "message":
                            source.readline()
                        shutil.copyfileobj(
                            source,
                            handles[name],
                            length=16 * 1024 * 1024,
                        )

            for handle in handles.values():
                handle.close()
            handles.clear()
            for name, temporary_path in temporary_paths.items():
                temporary_path.replace(final_paths[name])
        except Exception:
            for handle in handles.values():
                handle.close()
            for temporary_path in temporary_paths.values():
                temporary_path.unlink(missing_ok=True)
            raise

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        """Close files and commit or discard temporary outputs."""
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
