#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Index transcript annotations and annotated coding ORFs.
# Input: Source genePred annotation.
# Output: Transcript metadata and annotated-ORF overlap index.

"""Index transcript annotations and annotated coding ORFs."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from .config import _BUFFER_BYTES, ANNOTATION_BIN_SIZE

if TYPE_CHECKING:
    from .models import Candidate


@dataclass(frozen=True, slots=True)
class TranscriptMeta:
    """Store compact annotation metadata required during clustering."""

    transcript_id: bytes
    gene_id: bytes
    chrom: bytes
    strand: bytes
    transcript_length: int
    annotation_order: int
    is_coding: bool
    cds_blocks: tuple[tuple[int, int], ...]


@dataclass(frozen=True, slots=True)
class AnnotatedORF:
    """Store one unique genePred-derived annotated coding structure."""

    annotated_orf_id: bytes
    transcript_id: bytes
    gene_id: bytes
    chrom: bytes
    strand: bytes
    blocks: tuple[tuple[int, int], ...]
    nt_length: int
    block_offsets: tuple[int, ...]

    @property
    def genomic_start(self) -> int:
        """Return the leftmost coding coordinate."""
        return self.blocks[0][0]

    @property
    def genomic_end(self) -> int:
        """Return the rightmost coding coordinate."""
        return self.blocks[-1][1]


@dataclass(frozen=True, slots=True)
class AnnotatedOverlap:
    """Describe one phase-compatible candidate-to-annotation overlap."""

    annotated_orf: AnnotatedORF
    relation: bytes
    overlap_nt: int
    overlap_codon: int


@dataclass(slots=True)
class AnnotationIndex:
    """Index transcript metadata, gene lifetimes, and annotated mORFs."""

    transcripts: dict[bytes, TranscriptMeta]
    gene_first_order: dict[bytes, int]
    gene_last_order: dict[bytes, int]
    annotated_orfs: tuple[AnnotatedORF, ...]
    annotated_bins: dict[tuple[bytes, bytes, int], tuple[int, ...]]
    safe_cut_after: tuple[bool, ...]
    transcript_count: int

    @classmethod
    def from_genepred(cls, path: str | Path) -> AnnotationIndex:
        """Build a compact index from genePred or genePredExt."""
        annotation_path = Path(path)
        transcripts: dict[bytes, TranscriptMeta] = {}
        gene_first_order: dict[bytes, int] = {}
        gene_last_order: dict[bytes, int] = {}
        annotation_order = 0

        with annotation_path.open("rb", buffering=_BUFFER_BYTES) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith(b"#"):
                    continue
                annotation_order += 1
                fields = line.split(b"\t")
                if len(fields) < 10:
                    raise ValueError(
                        f"Invalid genePred row at line {line_number}: expected at least 10 columns."
                    )

                transcript_id = fields[0].strip()
                chrom = fields[1].strip()
                strand = fields[2].strip()
                if strand not in {b"+", b"-"}:
                    raise ValueError(
                        f"Invalid transcript strand at line {line_number}: {strand!r}."
                    )
                if transcript_id in transcripts:
                    raise ValueError(
                        "Duplicate transcript identifier in annotation: "
                        + transcript_id.decode("utf-8", errors="replace")
                    )

                try:
                    cds_start = int(fields[5])
                    cds_end = int(fields[6])
                    exon_count = int(fields[7])
                    exon_starts = cls._parse_int_list(fields[8])
                    exon_ends = cls._parse_int_list(fields[9])
                except ValueError as error:
                    raise ValueError(
                        f"Invalid integer field at genePred line {line_number}."
                    ) from error

                if (
                    len(exon_starts) != exon_count
                    or len(exon_ends) != exon_count
                    or exon_count < 1
                    or any(
                        exon_end <= exon_start
                        for exon_start, exon_end in zip(
                            exon_starts,
                            exon_ends,
                        )
                    )
                ):
                    raise ValueError(
                        "Invalid exon structure for transcript "
                        + transcript_id.decode("utf-8", errors="replace")
                    )

                gene_id = (
                    fields[11].strip()
                    if len(fields) >= 12 and fields[11].strip()
                    else transcript_id
                )
                gene_first_order.setdefault(gene_id, annotation_order)
                gene_last_order[gene_id] = annotation_order
                transcript_length = sum(
                    exon_end - exon_start
                    for exon_start, exon_end in zip(
                        exon_starts,
                        exon_ends,
                    )
                )
                is_coding = cds_end > cds_start
                cds_blocks: tuple[tuple[int, int], ...] = ()
                if is_coding:
                    cds_blocks = tuple(
                        (max(exon_start, cds_start), min(exon_end, cds_end))
                        for exon_start, exon_end in zip(
                            exon_starts,
                            exon_ends,
                        )
                        if min(exon_end, cds_end) > max(exon_start, cds_start)
                    )
                    if not cds_blocks:
                        raise ValueError(
                            "Coding transcript has no CDS blocks: "
                            + transcript_id.decode("utf-8", errors="replace")
                        )
                transcripts[transcript_id] = TranscriptMeta(
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    transcript_length=transcript_length,
                    annotation_order=annotation_order,
                    is_coding=is_coding,
                    cds_blocks=cds_blocks,
                )

        if not transcripts:
            raise ValueError(f"No transcript records found: {annotation_path}")

        difference = [0] * (annotation_order + 2)
        for gene_id, first_order in gene_first_order.items():
            last_order = gene_last_order[gene_id]
            if first_order < last_order:
                difference[first_order] += 1
                difference[last_order] -= 1
        active = 0
        safe_cut_after = [False] * (annotation_order + 1)
        for order in range(1, annotation_order + 1):
            active += difference[order]
            safe_cut_after[order] = active == 0

        annotated_orfs: list[AnnotatedORF] = []
        seen_structures: set[tuple[bytes, bytes, bytes, tuple[tuple[int, int], ...]]] = set()
        for transcript_meta in transcripts.values():
            if not transcript_meta.is_coding:
                continue
            structure_key = (
                transcript_meta.gene_id,
                transcript_meta.chrom,
                transcript_meta.strand,
                transcript_meta.cds_blocks,
            )
            if structure_key in seen_structures:
                continue
            seen_structures.add(structure_key)
            cds_length = sum(end - start for start, end in transcript_meta.cds_blocks)
            annotated_orfs.append(
                AnnotatedORF(
                    annotated_orf_id=(b"ANNOTATED_TRANSCRIPT:" + transcript_meta.transcript_id),
                    transcript_id=transcript_meta.transcript_id,
                    gene_id=transcript_meta.gene_id,
                    chrom=transcript_meta.chrom,
                    strand=transcript_meta.strand,
                    blocks=transcript_meta.cds_blocks,
                    nt_length=cds_length,
                    block_offsets=cls._block_offsets(
                        transcript_meta.cds_blocks,
                        transcript_meta.strand,
                    ),
                )
            )

        mutable_bins: dict[tuple[bytes, bytes, int], list[int]] = defaultdict(list)
        for annotated_index, annotated_orf in enumerate(annotated_orfs):
            first_bin = annotated_orf.genomic_start // ANNOTATION_BIN_SIZE
            last_bin = (annotated_orf.genomic_end - 1) // ANNOTATION_BIN_SIZE
            for bin_number in range(first_bin, last_bin + 1):
                mutable_bins[(annotated_orf.chrom, annotated_orf.strand, bin_number)].append(
                    annotated_index
                )

        return cls(
            transcripts=transcripts,
            gene_first_order=gene_first_order,
            gene_last_order=gene_last_order,
            annotated_orfs=tuple(annotated_orfs),
            annotated_bins={key: tuple(values) for key, values in mutable_bins.items()},
            safe_cut_after=tuple(safe_cut_after),
            transcript_count=annotation_order,
        )

    @staticmethod
    def _parse_int_list(value: bytes) -> tuple[int, ...]:
        """Parse a comma-separated integer byte field."""
        text = value.strip().rstrip(b",")
        if not text:
            return ()
        return tuple(int(item) for item in text.split(b",") if item)

    @staticmethod
    def _block_offsets(
        blocks: Sequence[tuple[int, int]],
        strand: bytes,
    ) -> tuple[int, ...]:
        """Return translation offsets aligned to genomic-order blocks."""
        lengths = [end - start for start, end in blocks]
        offsets = [0] * len(lengths)
        running = 0
        indices = range(len(lengths)) if strand == b"+" else range(len(lengths) - 1, -1, -1)
        for index in indices:
            offsets[index] = running
            running += lengths[index]
        return tuple(offsets)

    def overlapping_annotated_orfs(
        self,
        candidate: Candidate,
    ) -> Iterator[AnnotatedORF]:
        """Yield unique annotated ORFs overlapping a candidate span."""
        span_first_bin = candidate.blocks[0][0] // ANNOTATION_BIN_SIZE
        span_last_bin = (candidate.blocks[-1][1] - 1) // ANNOTATION_BIN_SIZE
        if span_first_bin == span_last_bin:
            indices: Iterable[int] = self.annotated_bins.get(
                (candidate.chrom, candidate.strand, span_first_bin),
                (),
            )
        else:
            unique_indices: set[int] = set()
            for block_start, block_end in candidate.blocks:
                first_bin = block_start // ANNOTATION_BIN_SIZE
                last_bin = (block_end - 1) // ANNOTATION_BIN_SIZE
                for bin_number in range(first_bin, last_bin + 1):
                    unique_indices.update(
                        self.annotated_bins.get(
                            (candidate.chrom, candidate.strand, bin_number),
                            (),
                        )
                    )
            indices = unique_indices

        candidate_start = candidate.blocks[0][0]
        candidate_end = candidate.blocks[-1][1]
        for annotated_index in indices:
            annotated_orf = self.annotated_orfs[annotated_index]
            if (
                annotated_orf.genomic_start < candidate_end
                and annotated_orf.genomic_end > candidate_start
            ):
                yield annotated_orf

    def is_safe_cut(self, annotation_order: int) -> bool:
        """Return whether no gene spans the cut after one transcript order."""
        return (
            0 < annotation_order < len(self.safe_cut_after)
            and self.safe_cut_after[annotation_order]
        )
