#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Define compact smORF cluster data models.
# Input: Candidate records and Kozak configuration.
# Output: Typed family, summary, and worker records.

"""Define compact smORF cluster data models."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any

from .annotation import TranscriptMeta
from .config import _BASE_BYTES, _DEFAULT_START_INDEX, _KOZAK_CACHE_SIZE
from .kozak import KozakModel


@dataclass(frozen=True, slots=True)
class _KozakSpec:
    """Store the pickle-safe part of a Kozak model."""

    rows: tuple[tuple[float, float, float, float], ...]
    name: bytes
    start_index: int

    @classmethod
    def from_model(cls, model: KozakModel | None) -> _KozakSpec | None:
        """Build a compact specification from a Kozak model."""
        if model is None:
            return None
        return cls(
            rows=tuple(tuple(float(value) for value in row) for row in model.rows),
            name=str(model.name).encode("utf-8"),
            start_index=int(model.start_index),
        )


@dataclass(frozen=True, slots=True)
class _ClusterConfig:
    """Store immutable clustering and filtering settings."""

    keep_start_codons: frozenset[bytes]
    min_aa: int
    max_aa: int
    keep_categories: frozenset[bytes]
    remove_categories: frozenset[bytes]
    require_sense: bool
    max_ambiguous_codons: int
    min_kozak_score: float
    start_rank: tuple[tuple[bytes, int], ...]
    kozak_spec: _KozakSpec | None


class _FastKozakScorer:
    """Score short Kozak contexts with precomputed per-base tables."""

    def __init__(self, spec: _KozakSpec) -> None:
        self.spec = spec
        excluded = {
            spec.start_index,
            spec.start_index + 1,
            spec.start_index + 2,
        }
        informative: list[tuple[int, tuple[float, ...]]] = []
        for position, row in enumerate(spec.rows):
            if position in excluded:
                continue
            row_min = min(row)
            row_max = max(row)
            if math.isclose(row_min, row_max, abs_tol=1e-15):
                continue
            minimum = max(row_min, 1e-300)
            maximum = max(row_max, 1e-300)
            denominator = math.log(maximum) - math.log(minimum)
            table = [-1.0] * 256
            for base_byte, probability in zip(_BASE_BYTES, row):
                if probability <= 0 or denominator <= 0:
                    value = 0.0
                else:
                    value = (math.log(max(probability, 1e-300)) - math.log(minimum)) / denominator
                    value = max(0.0, min(1.0, value))
                table[base_byte] = value
            informative.append((position, tuple(table)))
        self.informative = tuple(informative)
        self.informative_count = len(informative)

    @staticmethod
    def _clean(value: bytes) -> bytes:
        cleaned = value.strip().upper()
        return cleaned.replace(b"U", b"T") if b"U" in cleaned else cleaned

    @staticmethod
    def _parse_start_index(value: bytes) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _resolve_start_index(
        sequence: bytes,
        provided: bytes,
        start_codon: bytes,
    ) -> int | None:
        explicit = _FastKozakScorer._parse_start_index(provided)
        if (
            explicit is not None
            and 0 <= explicit <= len(sequence) - 3
            and (not start_codon or sequence[explicit : explicit + 3] == start_codon)
        ):
            return explicit
        if len(start_codon) == 3:
            candidates: list[int] = []
            start = 0
            while True:
                position = sequence.find(start_codon, start)
                if position < 0:
                    break
                candidates.append(position)
                start = position + 1
            if candidates:
                return min(
                    candidates,
                    key=lambda position: (
                        abs(position - _DEFAULT_START_INDEX),
                        position,
                    ),
                )
        if len(sequence) - 3 >= _DEFAULT_START_INDEX:
            return _DEFAULT_START_INDEX
        return None

    @lru_cache(maxsize=_KOZAK_CACHE_SIZE)
    def score(
        self,
        raw_sequence: bytes,
        raw_start_index: bytes,
        raw_start_codon: bytes,
    ) -> tuple[float | None, bytes, float]:
        """Return normalized score, level, and valid-base ratio."""
        sequence = self._clean(raw_sequence)
        start_codon = self._clean(raw_start_codon)
        anchor = self._resolve_start_index(
            sequence,
            raw_start_index,
            start_codon,
        )
        if anchor is None or not self.informative:
            return None, b"NA", 0.0

        model_offset = self.spec.start_index - anchor
        total = 0.0
        valid = 0
        for model_position, table in self.informative:
            sequence_position = model_position - model_offset
            if not 0 <= sequence_position < len(sequence):
                continue
            value = table[sequence[sequence_position]]
            if value < 0:
                continue
            total += value
            valid += 1

        valid_ratio = valid / self.informative_count if self.informative_count else 0.0
        if valid == 0 or valid_ratio < 0.5:
            return None, b"NA", valid_ratio
        score = total / valid
        if score >= 0.75:
            level = b"Strong"
        elif score >= 0.5:
            level = b"Moderate"
        else:
            level = b"Weak"
        return score, level, valid_ratio


@dataclass(frozen=True, slots=True)
class Candidate:
    """Store compact parsed ORF metadata for one active gene."""

    fields: tuple[bytes, ...]
    input_order: int
    transcript_meta: TranscriptMeta
    blocks: tuple[tuple[int, int], ...]
    oriented_blocks: tuple[tuple[int, int], ...]
    nt_length: int
    aa_length: int
    start_codon: bytes
    stop_codon: bytes
    kozak_score: float | None
    orf_id: bytes
    transcript_id: bytes
    gene_id: bytes
    chrom: bytes
    strand: bytes
    source_strand: bytes
    category: bytes
    stop_boundary: int

    @property
    def transcript_length(self) -> int:
        return self.transcript_meta.transcript_length

    @property
    def exact_key(
        self,
    ) -> tuple[bytes, bytes, bytes, tuple[tuple[int, int], ...]]:
        return (
            self.chrom,
            self.strand,
            self.source_strand,
            self.blocks,
        )

    @property
    def morf_key(
        self,
    ) -> tuple[bytes, bytes, tuple[tuple[int, int], ...]]:
        return self.chrom, self.strand, self.blocks

    def path_key(self, block_index: int = 0) -> tuple[Any, ...]:
        """Return a splice-suffix key from one oriented block."""
        block = self.oriented_blocks[block_index]
        terminal = block[1] if self.strand == b"+" else block[0]
        return (
            terminal,
            self.oriented_blocks[block_index + 1 :],
            self.nt_length % 3,
        )

    def suffix_extent(self, block_index: int) -> int:
        """Return the variable start-side boundary of a suffix block."""
        block = self.oriented_blocks[block_index]
        return block[0] if self.strand == b"+" else block[1]

    def start_extent(self) -> int:
        """Return the variable start-side boundary of the ORF first block."""
        block = self.oriented_blocks[0]
        return block[0] if self.strand == b"+" else block[1]


@dataclass(slots=True)
class ExactRepresentative:
    """Store one exact ORF representative and transcript duplicates."""

    primary: Candidate
    members: list[Candidate]


@dataclass(slots=True)
class Family:
    """Store one alternative-start family."""

    primary: ExactRepresentative
    representatives: list[ExactRepresentative]


@dataclass(slots=True)
class ClusterSummary:
    """Accumulate filtering, clustering, and runtime statistics."""

    input_orfs: int = 0
    basic_passed: int = 0
    basic_removed: int = 0
    annotated_overlap_removed: int = 0
    lnc_morf_removed: int = 0
    exact_duplicates_collapsed: int = 0
    alt_starts_collapsed: int = 0
    family_count: int = 0
    singleton_families: int = 0
    multi_member_families: int = 0
    genes_processed: int = 0
    effective_workers: int = 1
    removal_reasons: Counter[str] = field(default_factory=Counter)
    primary_categories: Counter[str] = field(default_factory=Counter)

    def merge(self, other: ClusterSummary) -> None:
        """Merge one independent shard summary."""
        self.input_orfs += other.input_orfs
        self.basic_passed += other.basic_passed
        self.basic_removed += other.basic_removed
        self.annotated_overlap_removed += other.annotated_overlap_removed
        self.lnc_morf_removed += other.lnc_morf_removed
        self.exact_duplicates_collapsed += other.exact_duplicates_collapsed
        self.alt_starts_collapsed += other.alt_starts_collapsed
        self.family_count += other.family_count
        self.singleton_families += other.singleton_families
        self.multi_member_families += other.multi_member_families
        self.genes_processed += other.genes_processed
        self.removal_reasons.update(other.removal_reasons)
        self.primary_categories.update(other.primary_categories)


class _UnionFind:
    """Compact disjoint-set implementation."""

    __slots__ = ("parent", "rank")

    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = bytearray(size)

    def find(self, item: int) -> int:
        parent = self.parent
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        parent = self.parent
        rank = self.rank
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if rank[left_root] < rank[right_root]:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        if rank[left_root] == rank[right_root]:
            rank[left_root] += 1


@dataclass(slots=True)
class _SuffixState:
    """Index prior longer suffixes for one alternative-start path."""

    anchor: int | None = None
    pending: list[tuple[int, int]] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class _RangeTask:
    """Describe one gene-safe input byte range."""

    shard_id: int
    start: int
    end: int
    temp_directory: str


@dataclass(frozen=True, slots=True)
class _RangeResult:
    """Return one worker shard result."""

    shard_id: int
    paths: tuple[str, str, str]
    summary: ClusterSummary
