#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-11-24
# Version: dev001
# Function: Audit smORFs overlapping annotated coding sequences.
# Input: Removed candidates, annotated CDS structures, P-sites, and optional external evidence.
# Output: Candidate-specific and shared CDS support with evidence provenance.

"""Audit annotated CDS overlaps without changing the independent smORF export."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

import numpy as np

from ..cluster.annotation import AnnotatedORF, AnnotationIndex
from .geometry import _geometry_blocks
from .input import _coding_nt_length, _required_columns, _smart_open
from .reporting import (
    MIN_PRIORITY_FRAME0_RATIO,
    MIN_PRIORITY_RPF,
    CountSummary,
    DensityProfiles,
    coding_positions,
    summarize_counts,
)


CONTEXT_COLUMNS = (
    "annotated_cds_ids", "shared_cds_phase0_codons", "outside_cds_codons",
    "shared_cds_frame0_rpf", "annotated_competing_frame0_rpf",
    "outside_cds_rpf", "outside_cds_frame0_rpf", "outside_cds_frame0_codons",
    "outside_cds_sample_pass", "tis_exact_start", "tis_alternative_start",
    "tis_ltm_reads", "tis_source", "tis_nterminal_evidence",
    "junction_context", "peptide_evidence",
)
CDS_COLUMNS = (
    "orf_id", "alias_orf_ids", "gene_ids", "transcript_ids", "categories", "chrom", "strand",
    "start0", "aa_length", "exon_starts", "exon_ends", "cluster_overlap_relations",
    "sample", "group", "total_rpf", "frame0_rpf", "frame0_ratio", "frame0_codons",
    "frame0_segments", "cds_overlap_status",
) + CONTEXT_COLUMNS


@dataclass(frozen=True, slots=True)
class CDSMasks:
    """Mark bases and codons accountable by any overlapping annotated isoform."""

    parent_ids: tuple[str, ...]
    parent_starts: frozenset[int]
    shared_phase0: np.ndarray
    competing_phase0: np.ndarray
    outside_codons: np.ndarray


def cds_masks(
    positions: np.ndarray, strand: str, parents: Iterable[AnnotatedORF]
) -> CDSMasks:
    """Project CDS phases onto the candidate's physical coding bases."""
    covered = np.zeros(positions.size, dtype=bool)
    annotated_phase0 = np.zeros_like(covered)
    parent_ids: set[str] = set()
    parent_starts: set[int] = set()
    for parent in parents:
        touched = False
        for (start, end), offset in zip(parent.blocks, parent.block_offsets):
            indices = np.flatnonzero((positions >= start) & (positions < end))
            if not indices.size:
                continue
            touched = True
            covered[indices] = True
            bases = positions[indices]
            distance = bases - start if strand == "+" else end - 1 - bases
            annotated_phase0[indices[(distance + offset) % 3 == 0]] = True
        if touched:
            parent_ids.add(parent.annotated_orf_id.decode())
            parent_starts.add(parent.blocks[0][0] if strand == "+" else parent.blocks[-1][1] - 1)

    phase0 = annotated_phase0.reshape(-1, 3)
    competing = phase0.copy()
    competing[:, 0] = False
    return CDSMasks(
        tuple(sorted(parent_ids)), frozenset(parent_starts), phase0[:, 0],
        competing, ~covered.reshape(-1, 3).any(axis=1),
    )


@dataclass(frozen=True, slots=True)
class CDSSampleEvidence:
    """Keep annotated, candidate-exclusive, and body P-site support distinct."""

    body: CountSummary
    shared_frame0_rpf: float
    competing_phase0_rpf: float
    outside: CountSummary
    outside_pass: bool

    @property
    def body_pass(self) -> bool:
        """Check minimum distributed body support for a start-site claim."""
        return (
            self.body.total >= MIN_PRIORITY_RPF
            and self.body.frame0_ratio >= MIN_PRIORITY_FRAME0_RATIO
            and self.body.frame0_codons >= 3
            and self.body.frame0_segments >= 2
        )


def cds_sample_evidence(counts: np.ndarray, masks: CDSMasks) -> CDSSampleEvidence:
    """Summarize stop-excluded candidate and host signals in one library."""
    outside = summarize_counts(counts[masks.outside_codons])
    outside_codons = int(masks.outside_codons.sum())
    return CDSSampleEvidence(
        body=summarize_counts(counts),
        shared_frame0_rpf=float(counts[masks.shared_phase0, 0].sum()),
        competing_phase0_rpf=float(counts[masks.competing_phase0].sum()),
        outside=outside,
        outside_pass=(
            outside_codons >= 2
            and outside.total >= MIN_PRIORITY_RPF
            and outside.frame0_ratio >= MIN_PRIORITY_FRAME0_RATIO
            and outside.frame0_codons >= min(3, outside_codons)
        ),
    )


def overlap_status(
    samples: Sequence[CDSSampleEvidence], required_samples: int, alternative_start: bool
) -> str:
    """Classify independent extent and exact initiation separately from shared density."""
    if alternative_start and sum(sample.body_pass for sample in samples) >= required_samples:
        return "AlternativeStartSupported"
    if sum(sample.outside_pass for sample in samples) >= required_samples:
        return "IndependentExtentSupported"
    if any(sample.shared_frame0_rpf > 0 for sample in samples):
        return "SharedCDSSignal"
    return "NoIndependentSupport"


def _read_rows(path: str | Path, columns: Sequence[str]) -> Iterator[dict[str, str]]:
    """Read a normalized input table, rejecting malformed rows."""
    with _smart_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _required_columns(reader.fieldnames or (), columns, str(path))
        for line, row in enumerate(reader, start=2):
            if None in row or any(row[column] is None for column in columns):
                raise ValueError("Malformed TSV row at {}:{}.".format(path, line))
            yield row


@dataclass(frozen=True, slots=True)
class InitiationEvidence:
    """Record an external initiation observation for a physical genomic start."""

    gene_id: str
    ltm_reads: float
    source: str
    nterminal: str


class EvidenceContext:
    """Look up annotated CDS, external initiation, RNA junction, and peptide evidence."""

    def __init__(
        self,
        annotation_path: str | Path | None,
        tis_path: str | Path | None,
        splice_path: str | Path | None,
        peptide_path: str | Path | None,
    ) -> None:
        self.annotation = AnnotationIndex.from_genepred(annotation_path) if annotation_path else None
        self.has_tis = bool(tis_path)
        self.tis: dict[tuple[str, str, int], list[InitiationEvidence]] = defaultdict(list)
        self.junctions: dict[str, set[str]] = defaultdict(set)
        self.peptides: dict[str, set[str]] = defaultdict(set)
        if tis_path:
            for row in _read_rows(tis_path, ("chrom", "strand", "start0", "ltm_reads")):
                reads = float(row["ltm_reads"])
                delta = float(row.get("ltm_minus_chx") or reads)
                start0 = int(row["start0"])
                if (
                    not math.isfinite(reads) or not math.isfinite(delta) or reads < 0
                    or row["strand"] not in {"+", "-"} or start0 < 0 or not row["chrom"]
                ):
                    raise ValueError("Invalid initiation evidence at {}:{}.".format(row["chrom"], start0))
                if reads > 0 and delta > 0:
                    source = "|".join(filter(None, (
                        Path(tis_path).name, row.get("source"), row.get("source_class"), row.get("tissue"),
                    )))
                    self.tis[(row["chrom"], row["strand"], start0)].append(InitiationEvidence(
                        row.get("gene_id") or "", reads, source, row.get("nterminal_peptides") or ""
                    ))
        if splice_path:
            for row in _read_rows(splice_path, ("transcript_id", "supporting_reads", "source")):
                supporting_reads = int(row["supporting_reads"])
                if supporting_reads < 0 or not row["transcript_id"] or not row["source"]:
                    raise ValueError("Invalid junction evidence for {}.".format(row["transcript_id"]))
                self.junctions[row["transcript_id"]].add(
                    "{}|{}|{}".format(row["transcript_id"], row["source"], supporting_reads)
                )
        if peptide_path:
            for row in _read_rows(peptide_path, ("orf_id", "evidence_id")):
                if not row["orf_id"] or not row["evidence_id"]:
                    raise ValueError("Peptide evidence requires nonempty orf_id and evidence_id.")
                self.peptides[row["orf_id"]].add(row["evidence_id"])

    def describe(
        self,
        chrom: str,
        strand: str,
        blocks: tuple[tuple[int, int], ...],
        positions: np.ndarray,
        genes: set[str],
        transcripts: set[str],
        orfs: set[str],
    ) -> tuple[CDSMasks | None, dict[str, object]]:
        """Match physical coordinates and preserve provenance across aliases."""
        masks = None
        if self.annotation is not None:
            parents = self.annotation.annotated_orfs_at(chrom.encode(), strand.encode(), blocks)
            masks = cds_masks(positions, strand, parents)
        start0 = int(positions[0])
        matches = [
            hit for hit in self.tis.get((chrom, strand, start0), ())
            if not hit.gene_id or hit.gene_id in genes
        ]
        alternative_start = bool(matches and masks and masks.parent_ids and start0 not in masks.parent_starts)
        metadata: dict[str, object] = {
            "tis_exact_start": int(bool(matches)) if self.has_tis else "",
            "tis_alternative_start": int(alternative_start) if self.has_tis and masks is not None else "",
            "tis_ltm_reads": max((hit.ltm_reads for hit in matches), default=""),
            "tis_source": ",".join(sorted({hit.source for hit in matches})),
            "tis_nterminal_evidence": ",".join(sorted({hit.nterminal for hit in matches if hit.nterminal})),
            "junction_context": ",".join(sorted({
                value for transcript in transcripts for value in self.junctions.get(transcript, ())
            })),
            "peptide_evidence": ",".join(sorted({
                value for orf in orfs for value in self.peptides.get(orf, ())
            })),
        }
        if masks is not None:
            metadata.update(
                annotated_cds_ids=",".join(masks.parent_ids),
                shared_cds_phase0_codons=int(masks.shared_phase0.sum()),
                outside_cds_codons=int(masks.outside_codons.sum()),
            )
        return masks, metadata


def cds_sample_columns(evidence: CDSSampleEvidence) -> dict[str, float | int]:
    """Serialize the sample-specific part of annotated-CDS evidence."""
    return {
        "shared_cds_frame0_rpf": evidence.shared_frame0_rpf,
        "annotated_competing_frame0_rpf": evidence.competing_phase0_rpf,
        "outside_cds_rpf": evidence.outside.total,
        "outside_cds_frame0_rpf": evidence.outside.frame0,
        "outside_cds_frame0_codons": evidence.outside.frame0_codons,
        "outside_cds_sample_pass": int(evidence.outside_pass),
    }


@dataclass(slots=True)
class _RemovedCandidate:
    """Combine synonymous aliases of an identical removed coding path."""

    row: dict[str, str]
    orfs: set[str] = field(default_factory=set)
    genes: set[str] = field(default_factory=set)
    transcripts: set[str] = field(default_factory=set)
    categories: set[str] = field(default_factory=set)
    relations: set[str] = field(default_factory=set)


def write_cds_overlap_audit(
    removed_path: str | Path,
    handle: TextIO,
    profiles: DensityProfiles,
    context: EvidenceContext,
    required_samples: int,
) -> Counter[str]:
    """Write one record per unique removed geometry and independent sample."""
    if context.annotation is None:
        raise ValueError("An annotated genePred is required for the CDS overlap audit.")
    columns = (
        "orf_id", "gene_id", "transcript_id", "category", "chrom", "strand", "nt_length", "aa_length",
        "exon_starts", "exon_ends", "stop_codon", "completeness", "annotated_overlap_relation",
    )
    candidates: dict[tuple[str, str, tuple[tuple[int, int], ...]], _RemovedCandidate] = {}
    for row in _read_rows(removed_path, columns):
        relation = row["annotated_overlap_relation"]
        if not relation:
            continue
        blocks = _geometry_blocks(row)
        key = row["chrom"], row["strand"], blocks
        candidate = candidates.get(key)
        if candidate is None:
            candidate = candidates[key] = _RemovedCandidate(row)
        candidate.orfs.add(row["orf_id"])
        candidate.genes.add(row["gene_id"])
        candidate.transcripts.add(row["transcript_id"])
        candidate.categories.add(row["category"])
        candidate.relations.add(relation)

    writer = csv.DictWriter(handle, fieldnames=CDS_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    statuses: Counter[str] = Counter()
    for (chrom, strand, blocks), candidate in sorted(candidates.items()):
        row = candidate.row
        positions = coding_positions(blocks, strand, _coding_nt_length(row))
        if not positions.size:
            raise ValueError("Removed ORF lacks complete coding codons: {}.".format(row["orf_id"]))
        masks, metadata = context.describe(
            chrom, strand, blocks, positions, candidate.genes, candidate.transcripts, candidate.orfs,
        )
        assert masks is not None
        samples = [
            (sample, cds_sample_evidence(counts, masks) if counts is not None else None)
            for sample, counts in profiles.counts(chrom, strand, positions)
        ]
        status = overlap_status(
            [evidence for _, evidence in samples if evidence is not None],
            required_samples, bool(metadata["tis_alternative_start"]),
        ) if masks.parent_ids else "NoAnnotatedCDSOverlap"
        statuses[status] += 1
        common = {
            "orf_id": row["orf_id"], "alias_orf_ids": ",".join(sorted(candidate.orfs)),
            "gene_ids": ",".join(sorted(candidate.genes)), "transcript_ids": ",".join(sorted(candidate.transcripts)),
            "categories": ",".join(sorted(candidate.categories)), "chrom": chrom, "strand": strand,
            "start0": int(positions[0]), "aa_length": row["aa_length"],
            "exon_starts": row["exon_starts"], "exon_ends": row["exon_ends"],
            "cluster_overlap_relations": ",".join(sorted(candidate.relations)), "cds_overlap_status": status,
            **metadata,
        }
        for sample, evidence in samples:
            result = {**common, "sample": sample, "group": profiles.groups[sample]}
            if evidence is not None:
                result.update(
                    total_rpf=evidence.body.total, frame0_rpf=evidence.body.frame0,
                    frame0_ratio=evidence.body.frame0_ratio, frame0_codons=evidence.body.frame0_codons,
                    frame0_segments=evidence.body.frame0_segments, **cds_sample_columns(evidence),
                )
            writer.writerow(result)
    return statuses
