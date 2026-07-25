#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.005
# Function: Cluster smORFs with indexed streaming and gene-safe parallel shards.
# Input: smorf_scanner message table and source genePred annotation.
# Output: Primary families, compact member mappings, removed ORFs, and summary.

"""High-performance gene-level smORF family clustering.

Workflow
--------
1. Parse genePred into compact transcript and gene-boundary indexes.
2. Stream the scanner table with compiled integer column indexes.
3. Apply structural filtering without constructing per-row dictionaries.
4. Remove same-gene lncORFs that exactly reproduce an annotated mORF.
5. Collapse exact genomic ORF duplicates across transcript isoforms.
6. Build alternative-start families with a splice-suffix index instead of
   quadratic all-pairs comparisons.
7. Optionally split large files at gene-safe byte boundaries and process
   independent ranges in parallel while preserving deterministic output order.

The public ``SmORFCluster`` interface is retained. The implementation does not
use the scanner ``frame`` column across transcripts because that frame is
transcript-relative. Family compatibility is derived from genomic exon blocks,
translation direction, terminal stop, and codon phase.
"""

from __future__ import annotations

import bisect
import heapq
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from types import TracebackType
from typing import Any, BinaryIO, Final

from utils.ribo.ArgsParser import progress_print

from .smorf_kozak import KozakModel

STOP_CODONS: Final[frozenset[bytes]] = frozenset({b"TAA", b"TAG", b"TGA"})
CLUSTER_REQUIRED_COLUMNS: Final[tuple[bytes, ...]] = (
    b"orf_id",
    b"gene_id",
    b"transcript_id",
    b"chrom",
    b"strand",
    b"source_strand",
    b"category",
    b"priority",
    b"start_codon",
    b"stop_codon",
    b"nt_length",
    b"aa_length",
    b"completeness",
    b"exon_count",
    b"exon_starts",
    b"exon_ends",
)
FILTER_COLUMNS: Final[tuple[bytes, ...]] = (
    b"kozak_pwm_name",
    b"kozak_pwm_score",
    b"kozak_pwm_level",
    b"kozak_valid_ratio",
    b"structure_status",
    b"filter_status",
    b"filter_reason",
)
CLUSTER_COLUMNS: Final[tuple[bytes, ...]] = (
    b"family_id",
    b"family_role",
    b"family_type",
    b"family_size",
    b"family_representative_count",
    b"family_transcript_count",
    b"family_exact_duplicate_count",
    b"family_alt_start_count",
    b"family_collapsed_count",
    b"family_categories",
    b"family_start_codons",
    b"primary_selection_rule",
)
MEMBER_COLUMNS: Final[tuple[bytes, ...]] = (
    b"family_id",
    b"primary_orf_id",
    b"representative_orf_id",
    b"member_orf_id",
    b"member_transcript_id",
    b"member_category",
    b"member_start_codon",
    b"member_aa_length",
    b"member_transcript_length",
    b"family_role",
    b"collapse_reason",
)
PRIMARY_SELECTION_RULE: Final[bytes] = (
    b"annotated>start_codon>aa_length>kozak_score>"
    b"transcript_length>input_order"
)
_BUFFER_BYTES: Final[int] = 8 * 1024 * 1024
_FLUSH_BYTES: Final[int] = 8 * 1024 * 1024
_MIN_PARALLEL_BYTES: Final[int] = 128 * 1024 * 1024
_KOZAK_CACHE_SIZE: Final[int] = 131_072
_DEFAULT_START_INDEX: Final[int] = 6
_BASE_BYTES: Final[tuple[int, ...]] = (65, 67, 71, 84)


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


@dataclass(slots=True)
class AnnotationIndex:
    """Index transcript metadata, gene lifetimes, and annotated mORFs."""

    transcripts: dict[bytes, TranscriptMeta]
    gene_first_order: dict[bytes, int]
    gene_last_order: dict[bytes, int]
    morf_keys_by_gene: dict[
        bytes,
        set[tuple[bytes, bytes, tuple[tuple[int, int], ...]]],
    ]
    safe_cut_after: tuple[bool, ...]
    transcript_count: int

    @classmethod
    def from_genepred(cls, path: str | Path) -> "AnnotationIndex":
        """Build a compact index from genePred or genePredExt."""
        annotation_path = Path(path)
        transcripts: dict[bytes, TranscriptMeta] = {}
        gene_first_order: dict[bytes, int] = {}
        gene_last_order: dict[bytes, int] = {}
        morf_keys_by_gene: dict[
            bytes,
            set[tuple[bytes, bytes, tuple[tuple[int, int], ...]]],
        ] = defaultdict(set)
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
                        f"Invalid genePred row at line {line_number}: "
                        "expected at least 10 columns."
                    )

                transcript_id = fields[0].strip()
                chrom = fields[1].strip()
                strand = fields[2].strip()
                if strand not in {b"+", b"-"}:
                    raise ValueError(
                        f"Invalid transcript strand at line {line_number}: "
                        f"{strand!r}."
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
                    morf_keys_by_gene[gene_id].add(
                        (chrom, strand, cds_blocks)
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

        return cls(
            transcripts=transcripts,
            gene_first_order=gene_first_order,
            gene_last_order=gene_last_order,
            morf_keys_by_gene=dict(morf_keys_by_gene),
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

    def is_safe_cut(self, annotation_order: int) -> bool:
        """Return whether no gene spans the cut after one transcript order."""
        return (
            0 < annotation_order < len(self.safe_cut_after)
            and self.safe_cut_after[annotation_order]
        )


@dataclass(frozen=True, slots=True)
class _KozakSpec:
    """Store the pickle-safe part of a Kozak model."""

    rows: tuple[tuple[float, float, float, float], ...]
    name: bytes
    start_index: int

    @classmethod
    def from_model(cls, model: KozakModel | None) -> "_KozakSpec | None":
        """Build a compact specification from a Kozak model."""
        if model is None:
            return None
        return cls(
            rows=tuple(
                tuple(float(value) for value in row)
                for row in model.rows
            ),
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


@dataclass(frozen=True, slots=True)
class _ColumnPlan:
    """Compile table column names into direct integer positions."""

    input_header: tuple[bytes, ...]
    filter_header: tuple[bytes, ...]
    family_header: tuple[bytes, ...]
    input_count: int
    filter_count: int
    family_count: int
    lookup: dict[bytes, int]
    family_lookup: dict[bytes, int]

    @classmethod
    def from_header(cls, header: Sequence[bytes]) -> "_ColumnPlan":
        """Build and validate one column plan."""
        input_header = tuple(header)
        lookup = {name: index for index, name in enumerate(input_header)}
        missing = [name for name in CLUSTER_REQUIRED_COLUMNS if name not in lookup]
        if missing:
            raise ValueError(
                "ORF message table is missing clustering column(s): "
                + ", ".join(
                    item.decode("utf-8", errors="replace")
                    for item in missing
                )
            )

        filter_header = list(input_header)
        filter_lookup = dict(lookup)
        for name in FILTER_COLUMNS:
            if name not in filter_lookup:
                filter_lookup[name] = len(filter_header)
                filter_header.append(name)

        family_header = list(filter_header)
        family_lookup = dict(filter_lookup)
        for name in CLUSTER_COLUMNS:
            if name not in family_lookup:
                family_lookup[name] = len(family_header)
                family_header.append(name)

        return cls(
            input_header=input_header,
            filter_header=tuple(filter_header),
            family_header=tuple(family_header),
            input_count=len(input_header),
            filter_count=len(filter_header),
            family_count=len(family_header),
            lookup=filter_lookup,
            family_lookup=family_lookup,
        )

    def index(self, name: bytes) -> int:
        """Return one filter-table field position."""
        return self.lookup[name]


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
                    value = (
                        math.log(max(probability, 1e-300))
                        - math.log(minimum)
                    ) / denominator
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
            and (
                not start_codon
                or sequence[explicit:explicit + 3] == start_codon
            )
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
        if _DEFAULT_START_INDEX <= len(sequence) - 3:
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

        valid_ratio = (
            valid / self.informative_count
            if self.informative_count
            else 0.0
        )
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
            self.oriented_blocks[block_index + 1:],
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

    def merge(self, other: "ClusterSummary") -> None:
        """Merge one independent shard summary."""
        self.input_orfs += other.input_orfs
        self.basic_passed += other.basic_passed
        self.basic_removed += other.basic_removed
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


class _ShardWriter:
    """Write one headerless worker shard with large byte buffers."""

    def __init__(self, directory: Path, shard_id: int) -> None:
        self.paths = {
            "family": directory / f"family.{shard_id:04d}.part",
            "members": directory / f"members.{shard_id:04d}.part",
            "removed": directory / f"removed.{shard_id:04d}.part",
        }
        self.handles: dict[str, BinaryIO] = {}
        self.buffers = {
            "family": bytearray(),
            "members": bytearray(),
            "removed": bytearray(),
        }

    def __enter__(self) -> "_ShardWriter":
        self.handles = {
            name: path.open("wb", buffering=_BUFFER_BYTES)
            for name, path in self.paths.items()
        }
        return self

    def _append(self, name: str, line: bytes) -> None:
        buffer = self.buffers[name]
        buffer.extend(line)
        if len(buffer) >= _FLUSH_BYTES:
            self.handles[name].write(buffer)
            buffer.clear()

    def write_family(self, fields: Sequence[bytes]) -> None:
        self._append("family", b"\t".join(fields) + b"\n")

    def write_member(self, fields: Sequence[bytes]) -> None:
        self._append("members", b"\t".join(fields) + b"\n")

    def write_removed(self, fields: Sequence[bytes]) -> None:
        self._append("removed", b"\t".join(fields) + b"\n")

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        for name, buffer in self.buffers.items():
            if buffer:
                self.handles[name].write(buffer)
                buffer.clear()
        for handle in self.handles.values():
            handle.close()
        if exception_type is not None:
            for path in self.paths.values():
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass
        return False


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


_WORKER_ANNOTATION: AnnotationIndex | None = None
_WORKER_PLAN: _ColumnPlan | None = None
_WORKER_CONFIG: _ClusterConfig | None = None
_WORKER_INPUT: str | None = None


def _worker_initialize(
    annotation_path: str,
    input_path: str,
    header: tuple[bytes, ...],
    config: _ClusterConfig,
) -> None:
    """Initialize one long-lived clustering worker."""
    global _WORKER_ANNOTATION, _WORKER_PLAN, _WORKER_CONFIG, _WORKER_INPUT
    _WORKER_ANNOTATION = AnnotationIndex.from_genepred(annotation_path)
    _WORKER_PLAN = _ColumnPlan.from_header(header)
    _WORKER_CONFIG = config
    _WORKER_INPUT = input_path


def _process_range_worker(task: _RangeTask) -> _RangeResult:
    """Process one independent gene-safe byte range."""
    if (
        _WORKER_ANNOTATION is None
        or _WORKER_PLAN is None
        or _WORKER_CONFIG is None
        or _WORKER_INPUT is None
    ):
        raise RuntimeError("smorf_cluster worker was not initialized.")
    writer = _ShardWriter(Path(task.temp_directory), task.shard_id)
    engine = _ClusterEngine(
        annotation=_WORKER_ANNOTATION,
        plan=_WORKER_PLAN,
        config=_WORKER_CONFIG,
    )
    with writer:
        summary = engine.process_range(
            input_path=_WORKER_INPUT,
            start=task.start,
            end=task.end,
            writer=writer,
        )
    return _RangeResult(
        shard_id=task.shard_id,
        paths=(
            str(writer.paths["family"]),
            str(writer.paths["members"]),
            str(writer.paths["removed"]),
        ),
        summary=summary,
    )


class _ClusterEngine:
    """Process one complete-gene input range."""

    def __init__(
        self,
        annotation: AnnotationIndex,
        plan: _ColumnPlan,
        config: _ClusterConfig,
    ) -> None:
        self.annotation = annotation
        self.plan = plan
        self.config = config
        self.start_rank = dict(config.start_rank)
        self.kozak_scorer = (
            _FastKozakScorer(config.kozak_spec)
            if config.kozak_spec is not None
            else None
        )

        lookup = plan.lookup
        self.i_orf_id = lookup[b"orf_id"]
        self.i_gene_id = lookup[b"gene_id"]
        self.i_transcript_id = lookup[b"transcript_id"]
        self.i_chrom = lookup[b"chrom"]
        self.i_strand = lookup[b"strand"]
        self.i_source_strand = lookup[b"source_strand"]
        self.i_category = lookup[b"category"]
        self.i_priority = lookup[b"priority"]
        self.i_start_codon = lookup[b"start_codon"]
        self.i_stop_codon = lookup[b"stop_codon"]
        self.i_nt_length = lookup[b"nt_length"]
        self.i_aa_length = lookup[b"aa_length"]
        self.i_completeness = lookup[b"completeness"]
        self.i_exon_count = lookup[b"exon_count"]
        self.i_exon_starts = lookup[b"exon_starts"]
        self.i_exon_ends = lookup[b"exon_ends"]
        self.i_ambiguous = lookup.get(b"ambiguous_codon_count", -1)
        self.i_kozak_seq = lookup.get(b"kozak_seq", -1)
        self.i_kozak_start = lookup.get(b"kozak_start_index", -1)
        self.i_kozak_name = lookup[b"kozak_pwm_name"]
        self.i_kozak_score = lookup[b"kozak_pwm_score"]
        self.i_kozak_level = lookup[b"kozak_pwm_level"]
        self.i_kozak_valid = lookup[b"kozak_valid_ratio"]
        self.i_structure = lookup[b"structure_status"]
        self.i_filter_status = lookup[b"filter_status"]
        self.i_filter_reason = lookup[b"filter_reason"]

    @staticmethod
    def _parse_int(value: bytes) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _parse_blocks(
        starts_value: bytes,
        ends_value: bytes,
        exon_count: int | None,
    ) -> tuple[tuple[int, int], ...] | None:
        try:
            starts = AnnotationIndex._parse_int_list(starts_value)
            ends = AnnotationIndex._parse_int_list(ends_value)
        except ValueError:
            return None
        if (
            not starts
            or len(starts) != len(ends)
            or exon_count != len(starts)
        ):
            return None
        blocks = tuple(zip(starts, ends))
        if any(end <= start for start, end in blocks):
            return None
        return blocks

    @staticmethod
    def _unique_reasons(reasons: Iterable[bytes]) -> list[bytes]:
        return list(dict.fromkeys(reasons))

    def _annotate_kozak(
        self,
        fields: list[bytes],
        start_codon: bytes,
    ) -> float | None:
        scorer = self.kozak_scorer
        if scorer is None:
            fields[self.i_kozak_name] = b"none"
            fields[self.i_kozak_score] = b"NA"
            fields[self.i_kozak_level] = b"NA"
            fields[self.i_kozak_valid] = b"NA"
            return None
        sequence = fields[self.i_kozak_seq] if self.i_kozak_seq >= 0 else b""
        start_index = (
            fields[self.i_kozak_start]
            if self.i_kozak_start >= 0
            else b""
        )
        score, level, valid_ratio = scorer.score(
            sequence,
            start_index,
            start_codon,
        )
        fields[self.i_kozak_name] = scorer.spec.name
        fields[self.i_kozak_score] = (
            b"NA" if score is None else f"{score:.6f}".encode("ascii")
        )
        fields[self.i_kozak_level] = level
        fields[self.i_kozak_valid] = f"{valid_ratio:.6f}".encode("ascii")
        return score

    def _parse_candidate(
        self,
        raw_fields: list[bytes],
        input_order: int,
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> Candidate | None:
        summary.input_orfs += 1
        if len(raw_fields) != self.plan.input_count:
            fields = raw_fields[:self.plan.input_count]
            fields.extend(
                [b""] * (self.plan.filter_count - len(fields))
            )
            fields[self.i_structure] = b"FAIL"
            fields[self.i_filter_status] = b"FAIL"
            fields[self.i_filter_reason] = b"invalid_column_count"
            summary.basic_removed += 1
            summary.removal_reasons["invalid_column_count"] += 1
            writer.write_removed(fields)
            return None

        fields = raw_fields
        if len(fields) < self.plan.filter_count:
            fields.extend([b""] * (self.plan.filter_count - len(fields)))

        reasons: list[bytes] = []
        structure_reasons: list[bytes] = []
        source_strand = fields[self.i_source_strand].strip().lower()
        priority = fields[self.i_priority].strip()
        completeness = fields[self.i_completeness].strip()
        category = fields[self.i_category].strip()
        start_codon = fields[self.i_start_codon].strip().upper().replace(b"U", b"T")
        stop_codon = fields[self.i_stop_codon].strip().upper().replace(b"U", b"T")
        transcript_id = fields[self.i_transcript_id].strip()
        gene_id = fields[self.i_gene_id].strip()
        chrom = fields[self.i_chrom].strip()
        strand = fields[self.i_strand].strip()

        if source_strand not in {b"sense", b"antisense"}:
            reasons.append(b"invalid_source_strand")
        elif self.config.require_sense and source_strand != b"sense":
            reasons.append(b"non_sense_strand")
        if not priority:
            reasons.append(b"missing_priority")
        if not completeness:
            reasons.append(b"missing_completeness")
        elif completeness != b"complete":
            reasons.append(b"incomplete_orf:" + completeness)
        if not category:
            reasons.append(b"missing_category")
        else:
            if category in self.config.remove_categories:
                reasons.append(b"removed_category:" + category)
            if (
                self.config.keep_categories
                and category not in self.config.keep_categories
            ):
                reasons.append(b"category_not_allowed:" + category)
        if not start_codon:
            reasons.append(b"missing_start_codon")
        elif (
            self.config.keep_start_codons
            and start_codon not in self.config.keep_start_codons
        ):
            reasons.append(b"start_codon_not_allowed:" + start_codon)

        aa_length = self._parse_int(fields[self.i_aa_length])
        nt_length = self._parse_int(fields[self.i_nt_length])
        exon_count = self._parse_int(fields[self.i_exon_count])
        if aa_length is None:
            reasons.append(b"invalid_aa_length")
        else:
            if aa_length < self.config.min_aa:
                reasons.append(b"too_short")
            if aa_length > self.config.max_aa:
                reasons.append(b"too_long")
        if nt_length is None or nt_length <= 0:
            structure_reasons.append(b"invalid_nt_length")
        elif aa_length is not None and nt_length != aa_length * 3 + 3:
            structure_reasons.append(b"inconsistent_orf_length")
        if completeness == b"complete" and stop_codon not in STOP_CODONS:
            structure_reasons.append(b"invalid_stop_codon")

        if self.i_ambiguous >= 0:
            ambiguous = self._parse_int(fields[self.i_ambiguous])
            if ambiguous is None or ambiguous < 0:
                structure_reasons.append(b"invalid_ambiguous_codon_count")
            elif ambiguous > self.config.max_ambiguous_codons:
                structure_reasons.append(b"too_many_ambiguous_codons")

        blocks = self._parse_blocks(
            fields[self.i_exon_starts],
            fields[self.i_exon_ends],
            exon_count,
        )
        if blocks is None:
            structure_reasons.append(b"invalid_exon_blocks")
        elif nt_length is not None and sum(end - start for start, end in blocks) != nt_length:
            structure_reasons.append(b"inconsistent_block_length")

        transcript_meta = self.annotation.transcripts.get(transcript_id)
        if transcript_meta is None:
            reasons.append(b"transcript_not_in_annotation")
        else:
            if gene_id != transcript_meta.gene_id:
                reasons.append(b"annotation_gene_mismatch")
            if chrom != transcript_meta.chrom:
                reasons.append(b"annotation_chrom_mismatch")
            expected_strand = transcript_meta.strand
            if source_strand == b"antisense":
                expected_strand = b"-" if expected_strand == b"+" else b"+"
            if strand != expected_strand:
                reasons.append(b"annotation_strand_mismatch")

        kozak_score = self._annotate_kozak(fields, start_codon)
        if self.config.kozak_spec is not None and self.config.min_kozak_score > 0:
            if kozak_score is None:
                reasons.append(b"insufficient_kozak_context")
            elif kozak_score < self.config.min_kozak_score:
                reasons.append(b"weak_kozak_pwm")

        reasons.extend(structure_reasons)
        unique_reasons = self._unique_reasons(reasons)
        fields[self.i_structure] = b"FAIL" if structure_reasons else b"PASS"
        if unique_reasons:
            fields[self.i_filter_status] = b"FAIL"
            fields[self.i_filter_reason] = b";".join(unique_reasons)
            summary.basic_removed += 1
            summary.removal_reasons.update(
                reason.decode("utf-8", errors="replace")
                for reason in unique_reasons
            )
            writer.write_removed(fields)
            return None

        fields[self.i_filter_status] = b"PASS"
        fields[self.i_filter_reason] = b"PASS"
        summary.basic_passed += 1
        assert transcript_meta is not None
        assert blocks is not None
        assert nt_length is not None
        assert aa_length is not None
        oriented_blocks = blocks if strand == b"+" else tuple(reversed(blocks))
        terminal_block = oriented_blocks[-1]
        stop_boundary = terminal_block[1] if strand == b"+" else terminal_block[0]
        return Candidate(
            fields=tuple(fields),
            input_order=input_order,
            transcript_meta=transcript_meta,
            blocks=blocks,
            oriented_blocks=oriented_blocks,
            nt_length=nt_length,
            aa_length=aa_length,
            start_codon=start_codon,
            stop_codon=stop_codon,
            kozak_score=kozak_score,
            orf_id=fields[self.i_orf_id],
            transcript_id=transcript_id,
            gene_id=gene_id,
            chrom=chrom,
            strand=strand,
            source_strand=source_strand,
            category=category,
            stop_boundary=stop_boundary,
        )

    def _iter_transcript_groups(
        self,
        input_path: str | Path,
        start: int,
        end: int,
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> Iterator[tuple[TranscriptMeta, list[Candidate]]]:
        """Stream and filter contiguous transcript blocks in one byte range."""
        current_transcript_id: bytes | None = None
        current_meta: TranscriptMeta | None = None
        current_candidates: list[Candidate] = []
        last_annotation_order = 0
        input_order = 0

        with Path(input_path).open("rb", buffering=_BUFFER_BYTES) as handle:
            handle.seek(start)
            while handle.tell() < end:
                line = handle.readline()
                if not line:
                    break
                input_order += 1
                raw_fields = line.rstrip(b"\r\n").split(b"\t")
                transcript_id = (
                    raw_fields[self.i_transcript_id].strip()
                    if len(raw_fields) > self.i_transcript_id
                    else b""
                )

                if (
                    current_transcript_id is not None
                    and transcript_id != current_transcript_id
                ):
                    if current_meta is not None:
                        if current_meta.annotation_order <= last_annotation_order:
                            raise ValueError(
                                "ORF message table is not ordered like the source "
                                "genePred annotation at transcript "
                                + current_transcript_id.decode(
                                    "utf-8",
                                    errors="replace",
                                )
                            )
                        last_annotation_order = current_meta.annotation_order
                        yield current_meta, current_candidates
                    current_candidates = []
                    current_transcript_id = None
                    current_meta = None

                if current_transcript_id is None and transcript_id:
                    current_transcript_id = transcript_id
                    current_meta = self.annotation.transcripts.get(transcript_id)

                candidate = self._parse_candidate(
                    raw_fields=raw_fields,
                    input_order=input_order,
                    summary=summary,
                    writer=writer,
                )
                if candidate is not None:
                    if current_meta is None:
                        current_meta = candidate.transcript_meta
                        current_transcript_id = candidate.transcript_id
                    current_candidates.append(candidate)

            if current_transcript_id is not None and current_meta is not None:
                if current_meta.annotation_order <= last_annotation_order:
                    raise ValueError(
                        "ORF message table is not transcript-contiguous or "
                        "annotation-ordered at transcript "
                        + current_transcript_id.decode(
                            "utf-8",
                            errors="replace",
                        )
                    )
                yield current_meta, current_candidates

    def _iter_gene_groups(
        self,
        transcript_groups: Iterator[
            tuple[TranscriptMeta, list[Candidate]]
        ],
    ) -> Iterator[tuple[bytes, list[Candidate]]]:
        """Aggregate active genes with a last-order min-heap."""
        pending: dict[bytes, list[Candidate]] = {}
        close_heap: list[tuple[int, int, bytes]] = []

        def flush_before(order: int) -> Iterator[tuple[bytes, list[Candidate]]]:
            while close_heap and close_heap[0][0] < order:
                _last, _first, gene_id = heapq.heappop(close_heap)
                candidates = pending.pop(gene_id, None)
                if candidates is not None:
                    yield gene_id, candidates

        for transcript_meta, candidates in transcript_groups:
            yield from flush_before(transcript_meta.annotation_order)
            gene_id = transcript_meta.gene_id
            if gene_id not in pending:
                pending[gene_id] = []
                heapq.heappush(
                    close_heap,
                    (
                        self.annotation.gene_last_order[gene_id],
                        self.annotation.gene_first_order[gene_id],
                        gene_id,
                    ),
                )
            pending[gene_id].extend(candidates)
            if transcript_meta.annotation_order == self.annotation.gene_last_order[gene_id]:
                pending_candidates = pending.pop(gene_id)
                yield gene_id, pending_candidates

        while close_heap:
            _last, _first, gene_id = heapq.heappop(close_heap)
            candidates = pending.pop(gene_id, None)
            if candidates is not None:
                yield gene_id, candidates

    @staticmethod
    def _category_rank(category: bytes) -> int:
        if category in {b"annotated_ORF", b"annotated_mORF"}:
            return 0
        if category in {
            b"uORF",
            b"dORF",
            b"iORF",
            b"same_frame_iORF",
            b"emORF",
            b"overlap_uORF",
            b"overlap_dORF",
        }:
            return 1
        if category == b"lncORF":
            return 2
        if category == b"other_ORF":
            return 3
        return 4

    def _exact_rank(self, candidate: Candidate) -> tuple[Any, ...]:
        return (
            self._category_rank(candidate.category),
            -candidate.transcript_length,
            candidate.input_order,
            candidate.orf_id,
        )

    def _primary_rank(
        self,
        representative: ExactRepresentative,
    ) -> tuple[Any, ...]:
        candidate = representative.primary
        protected = (
            0
            if candidate.category in {b"annotated_ORF", b"annotated_mORF"}
            else 1
        )
        kozak_rank = (
            -candidate.kozak_score
            if candidate.kozak_score is not None
            else float("inf")
        )
        return (
            protected,
            self.start_rank.get(candidate.start_codon, len(self.start_rank)),
            -candidate.aa_length,
            kozak_rank,
            -candidate.transcript_length,
            candidate.input_order,
            candidate.orf_id,
        )

    def _remove_lnc_morf_matches(
        self,
        gene_id: bytes,
        candidates: list[Candidate],
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> list[Candidate]:
        morf_keys = self.annotation.morf_keys_by_gene.get(gene_id)
        if not morf_keys:
            return candidates
        retained: list[Candidate] = []
        for candidate in candidates:
            if (
                candidate.source_strand == b"sense"
                and candidate.category == b"lncORF"
                and candidate.morf_key in morf_keys
            ):
                fields = list(candidate.fields)
                fields[self.i_filter_status] = b"FAIL"
                fields[self.i_filter_reason] = b"lncORF_matches_annotated_mORF"
                writer.write_removed(fields)
                summary.lnc_morf_removed += 1
                summary.removal_reasons[
                    "lncORF_matches_annotated_mORF"
                ] += 1
            else:
                retained.append(candidate)
        return retained

    def _deduplicate_exact(
        self,
        candidates: list[Candidate],
        summary: ClusterSummary,
    ) -> list[ExactRepresentative]:
        grouped: dict[
            tuple[bytes, bytes, bytes, tuple[tuple[int, int], ...]],
            ExactRepresentative,
        ] = {}
        for candidate in candidates:
            key = candidate.exact_key
            representative = grouped.get(key)
            if representative is None:
                grouped[key] = ExactRepresentative(
                    primary=candidate,
                    members=[candidate],
                )
                continue
            representative.members.append(candidate)
            if self._exact_rank(candidate) < self._exact_rank(representative.primary):
                representative.primary = candidate
            summary.exact_duplicates_collapsed += 1
        representatives = list(grouped.values())
        representatives.sort(key=lambda item: item.primary.input_order)
        return representatives

    @staticmethod
    def _suffix_compatible(
        candidate: Candidate,
        registered_extent: int,
    ) -> bool:
        current_extent = candidate.start_extent()
        if candidate.strand == b"+":
            return registered_extent <= current_extent
        return registered_extent >= current_extent

    def _cluster_bucket_indexed(
        self,
        bucket_indices: list[int],
        representatives: list[ExactRepresentative],
        union_find: _UnionFind,
    ) -> None:
        """Union one stop bucket with an indexed splice-suffix algorithm."""
        ordered = sorted(
            bucket_indices,
            key=lambda index: (
                -representatives[index].primary.nt_length,
                representatives[index].primary.input_order,
            ),
        )
        states: dict[tuple[Any, ...], _SuffixState] = {}

        for index in ordered:
            candidate = representatives[index].primary
            own_key = candidate.path_key(0)
            state = states.get(own_key)
            if state is None:
                state = _SuffixState()
                states[own_key] = state
            else:
                if state.anchor is not None:
                    union_find.union(index, state.anchor)
                if state.pending:
                    retained_pending: list[tuple[int, int]] = []
                    for pending_index, pending_extent in state.pending:
                        if self._suffix_compatible(candidate, pending_extent):
                            union_find.union(index, pending_index)
                        else:
                            retained_pending.append(
                                (pending_index, pending_extent)
                            )
                    state.pending = retained_pending
            if state.anchor is None:
                state.anchor = index

            for block_index in range(1, len(candidate.oriented_blocks)):
                suffix_key = candidate.path_key(block_index)
                suffix_state = states.get(suffix_key)
                if suffix_state is None:
                    suffix_state = _SuffixState()
                    states[suffix_key] = suffix_state
                suffix_state.pending.append(
                    (index, candidate.suffix_extent(block_index))
                )

    def _cluster_representatives(
        self,
        representatives: list[ExactRepresentative],
    ) -> list[Family]:
        if not representatives:
            return []
        union_find = _UnionFind(len(representatives))
        buckets: dict[
            tuple[bytes, bytes, bytes, int, bytes],
            list[int],
        ] = defaultdict(list)
        for index, representative in enumerate(representatives):
            candidate = representative.primary
            buckets[
                (
                    candidate.chrom,
                    candidate.strand,
                    candidate.source_strand,
                    candidate.stop_boundary,
                    candidate.stop_codon,
                )
            ].append(index)
        for bucket_indices in buckets.values():
            if len(bucket_indices) > 1:
                self._cluster_bucket_indexed(
                    bucket_indices,
                    representatives,
                    union_find,
                )

        components: dict[int, list[ExactRepresentative]] = defaultdict(list)
        for index, representative in enumerate(representatives):
            components[union_find.find(index)].append(representative)

        families: list[Family] = []
        for component in components.values():
            component.sort(key=lambda item: item.primary.input_order)
            primary = min(component, key=self._primary_rank)
            families.append(Family(primary=primary, representatives=component))
        families.sort(key=lambda item: item.primary.primary.input_order)
        return families

    def _ordered_start_codons(
        self,
        members: Sequence[Candidate],
    ) -> bytes:
        codons = {member.start_codon for member in members}
        return b",".join(
            sorted(
                codons,
                key=lambda codon: (
                    self.start_rank.get(codon, len(self.start_rank)),
                    codon,
                ),
            )
        )

    @staticmethod
    def _family_type(
        representative_count: int,
        exact_duplicate_count: int,
    ) -> bytes:
        has_alt = representative_count > 1
        has_exact = exact_duplicate_count > 0
        if has_alt and has_exact:
            return b"MIXED"
        if has_alt:
            return b"ALTERNATIVE_START"
        if has_exact:
            return b"TRANSCRIPT_DUPLICATE"
        return b"SINGLETON"

    @staticmethod
    def _all_members(family: Family) -> list[Candidate]:
        members = [
            member
            for representative in family.representatives
            for member in representative.members
        ]
        members.sort(key=lambda item: item.input_order)
        return members

    def _write_family(
        self,
        family: Family,
        local_family_number: int,
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> None:
        primary_candidate = family.primary.primary
        all_members = self._all_members(family)
        representative_count = len(family.representatives)
        exact_duplicate_count = sum(
            len(representative.members) - 1
            for representative in family.representatives
        )
        alt_start_count = representative_count - 1
        collapsed_count = len(all_members) - 1

        fields = list(primary_candidate.fields)
        fields.extend([b""] * (self.plan.family_count - len(fields)))
        family_lookup = self.plan.family_lookup
        fields[family_lookup[b"family_id"]] = str(local_family_number).encode("ascii")
        fields[family_lookup[b"family_role"]] = b"PRIMARY"
        fields[family_lookup[b"family_type"]] = self._family_type(
            representative_count,
            exact_duplicate_count,
        )
        fields[family_lookup[b"family_size"]] = str(len(all_members)).encode("ascii")
        fields[family_lookup[b"family_representative_count"]] = str(
            representative_count
        ).encode("ascii")
        fields[family_lookup[b"family_transcript_count"]] = str(
            len({member.transcript_id for member in all_members})
        ).encode("ascii")
        fields[family_lookup[b"family_exact_duplicate_count"]] = str(
            exact_duplicate_count
        ).encode("ascii")
        fields[family_lookup[b"family_alt_start_count"]] = str(
            alt_start_count
        ).encode("ascii")
        fields[family_lookup[b"family_collapsed_count"]] = str(
            collapsed_count
        ).encode("ascii")
        fields[family_lookup[b"family_categories"]] = b",".join(
            sorted({member.category for member in all_members})
        )
        fields[family_lookup[b"family_start_codons"]] = self._ordered_start_codons(
            all_members
        )
        fields[family_lookup[b"primary_selection_rule"]] = PRIMARY_SELECTION_RULE
        writer.write_family(fields)

        primary_id = family.primary.primary.orf_id
        ordered_representatives = [family.primary] + [
            representative
            for representative in family.representatives
            if representative is not family.primary
        ]
        family_number_bytes = str(local_family_number).encode("ascii")
        for representative in ordered_representatives:
            representative_id = representative.primary.orf_id
            ordered_members = [representative.primary] + [
                member
                for member in representative.members
                if member is not representative.primary
            ]
            for member in ordered_members:
                if member.orf_id == primary_id:
                    family_role = b"PRIMARY"
                    collapse_reason = b"PRIMARY"
                elif member.orf_id == representative_id:
                    family_role = b"COLLAPSED"
                    collapse_reason = b"ALTERNATIVE_START"
                else:
                    family_role = b"COLLAPSED"
                    collapse_reason = b"EXACT_TRANSCRIPT_DUPLICATE"
                writer.write_member(
                    (
                        family_number_bytes,
                        primary_id,
                        representative_id,
                        member.orf_id,
                        member.transcript_id,
                        member.category,
                        member.start_codon,
                        str(member.aa_length).encode("ascii"),
                        str(member.transcript_length).encode("ascii"),
                        family_role,
                        collapse_reason,
                    )
                )

        summary.family_count += 1
        summary.primary_categories[
            primary_candidate.category.decode("utf-8", errors="replace")
        ] += 1
        summary.alt_starts_collapsed += alt_start_count
        if len(all_members) == 1:
            summary.singleton_families += 1
        else:
            summary.multi_member_families += 1

    def process_range(
        self,
        input_path: str | Path,
        start: int,
        end: int,
        writer: _ShardWriter,
    ) -> ClusterSummary:
        """Filter and cluster one gene-safe input range."""
        summary = ClusterSummary()
        transcript_groups = self._iter_transcript_groups(
            input_path=input_path,
            start=start,
            end=end,
            summary=summary,
            writer=writer,
        )
        local_family_number = 0
        for gene_id, candidates in self._iter_gene_groups(transcript_groups):
            summary.genes_processed += 1
            candidates = self._remove_lnc_morf_matches(
                gene_id,
                candidates,
                summary,
                writer,
            )
            representatives = self._deduplicate_exact(candidates, summary)
            families = self._cluster_representatives(representatives)
            for family in families:
                local_family_number += 1
                self._write_family(
                    family,
                    local_family_number,
                    summary,
                    writer,
                )
        return summary


class SmORFCluster:
    """Filter and cluster smORFs with deterministic high-performance I/O."""

    def __init__(
        self,
        annotation: str | Path,
        keep_start_codons: str = "",
        min_aa: int = 8,
        max_aa: int = 10000,
        keep_categories: str = (
            "uORF,dORF,lncORF,iORF,same_frame_iORF,emORF,"
            "overlap_uORF,overlap_dORF,other_ORF,"
            "annotated_ORF,annotated_mORF"
        ),
        remove_categories: str = "antisense_ORF",
        require_sense: bool = True,
        max_ambiguous_codons: int = 0,
        kozak_model: KozakModel | None = None,
        min_kozak_score: float = 0.0,
        start_priority: str = "ATG",
        threads: int = 0,
    ) -> None:
        self.annotation_path = str(Path(annotation).resolve())
        self.annotation = AnnotationIndex.from_genepred(annotation)
        self.threads = int(threads)
        if self.threads < 0:
            raise ValueError("threads must be >= 0.")
        if min_aa < 1:
            raise ValueError("min_aa must be >= 1.")
        if max_aa < min_aa:
            raise ValueError("max_aa must be >= min_aa.")
        if max_ambiguous_codons < 0:
            raise ValueError("max_ambiguous_codons must be >= 0.")
        if not 0 <= min_kozak_score <= 1:
            raise ValueError("min_kozak_score must be in [0, 1].")
        if min_kozak_score > 0 and kozak_model is None:
            raise ValueError(
                "A Kozak model is required when min_kozak_score > 0."
            )

        keep_start_set = self._parse_bytes_set(keep_start_codons, normalize_codon=True)
        keep_category_set = self._parse_bytes_set(keep_categories)
        remove_category_set = self._parse_bytes_set(remove_categories)
        conflict = keep_category_set.intersection(remove_category_set)
        if conflict:
            raise ValueError(
                "Category cannot be both retained and removed: "
                + ", ".join(
                    item.decode("utf-8", errors="replace")
                    for item in sorted(conflict)
                )
            )
        start_order = self._parse_start_order(start_priority)
        self.config = _ClusterConfig(
            keep_start_codons=frozenset(keep_start_set),
            min_aa=int(min_aa),
            max_aa=int(max_aa),
            keep_categories=frozenset(keep_category_set),
            remove_categories=frozenset(remove_category_set),
            require_sense=bool(require_sense),
            max_ambiguous_codons=int(max_ambiguous_codons),
            min_kozak_score=float(min_kozak_score),
            start_rank=tuple((codon, rank) for rank, codon in enumerate(start_order)),
            kozak_spec=_KozakSpec.from_model(kozak_model),
        )

    @staticmethod
    def _parse_bytes_set(
        value: str,
        normalize_codon: bool = False,
    ) -> set[bytes]:
        output: set[bytes] = set()
        for item in str(value).split(","):
            text = item.strip()
            if not text:
                continue
            if normalize_codon:
                text = text.upper().replace("U", "T")
            output.add(text.encode("utf-8"))
        return output

    @staticmethod
    def _parse_start_order(value: str) -> list[bytes]:
        codons: list[bytes] = []
        seen: set[bytes] = set()
        for item in str(value).split(","):
            codon_text = item.strip().upper().replace("U", "T")
            if not codon_text:
                continue
            if len(codon_text) != 3 or set(codon_text).difference("ACGT"):
                raise ValueError(f"Invalid start-priority codon: {codon_text}")
            codon = codon_text.encode("ascii")
            if codon in STOP_CODONS:
                raise ValueError(
                    f"Stop codon cannot be a start-priority codon: {codon_text}"
                )
            if codon not in seen:
                seen.add(codon)
                codons.append(codon)
        if not codons:
            raise ValueError("Start-codon priority cannot be empty.")
        return codons

    @staticmethod
    def read_header(path: str | Path) -> tuple[bytes, ...]:
        """Read the first non-empty table header line."""
        with Path(path).open("rb", buffering=_BUFFER_BYTES) as handle:
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError(f"Empty ORF message table: {path}")
                stripped = line.strip()
                if stripped:
                    return tuple(stripped.split(b"\t"))

    @staticmethod
    def _data_start(path: str | Path) -> int:
        with Path(path).open("rb", buffering=_BUFFER_BYTES) as handle:
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError(f"Empty ORF message table: {path}")
                if line.strip():
                    return handle.tell()

    @staticmethod
    def _extract_field(line: bytes, index: int) -> bytes:
        """Extract one tab field without splitting the complete row."""
        start = 0
        for current in range(index):
            position = line.find(b"\t", start)
            if position < 0:
                return b""
            start = position + 1
        end = line.find(b"\t", start)
        if end < 0:
            end = len(line)
        return line[start:end].strip()

    def _effective_workers(self, file_size: int) -> int:
        if self.threads > 0:
            requested = self.threads
        else:
            slurm = os.environ.get("SLURM_CPUS_PER_TASK")
            requested = int(slurm) if slurm and slurm.isdigit() else (os.cpu_count() or 1)
        requested = max(1, requested)
        if file_size < _MIN_PARALLEL_BYTES:
            return 1
        return requested

    def _discover_ranges(
        self,
        input_path: str | Path,
        plan: _ColumnPlan,
        workers: int,
    ) -> list[tuple[int, int]]:
        """Find approximately balanced cuts that never split one gene."""
        input_file = Path(input_path)
        file_size = input_file.stat().st_size
        data_start = self._data_start(input_file)
        if workers <= 1 or file_size <= data_start:
            return [(data_start, file_size)]

        transcript_index = plan.lookup[b"transcript_id"]
        safe_offsets = [data_start]
        current_transcript: bytes | None = None
        current_order: int | None = None
        unknown_transcript = False

        with input_file.open("rb", buffering=_BUFFER_BYTES) as handle:
            handle.seek(data_start)
            while True:
                line_start = handle.tell()
                line = handle.readline()
                if not line:
                    break
                transcript_id = self._extract_field(
                    line.rstrip(b"\r\n"),
                    transcript_index,
                )
                if transcript_id != current_transcript:
                    if (
                        current_transcript is not None
                        and current_order is not None
                        and self.annotation.is_safe_cut(current_order)
                    ):
                        safe_offsets.append(line_start)
                    current_transcript = transcript_id
                    meta = self.annotation.transcripts.get(transcript_id)
                    if meta is None:
                        unknown_transcript = True
                        current_order = None
                    else:
                        current_order = meta.annotation_order

        if unknown_transcript:
            return [(data_start, file_size)]
        if safe_offsets[-1] != file_size:
            safe_offsets.append(file_size)
        safe_offsets = sorted(set(safe_offsets))
        if len(safe_offsets) <= 2:
            return [(data_start, file_size)]

        selected = [data_start]
        span = file_size - data_start
        for part in range(1, workers):
            target = data_start + span * part // workers
            position = bisect.bisect_left(safe_offsets, target)
            candidates = []
            if position < len(safe_offsets):
                candidates.append(safe_offsets[position])
            if position > 0:
                candidates.append(safe_offsets[position - 1])
            boundary = min(candidates, key=lambda value: abs(value - target))
            if boundary > selected[-1] and boundary < file_size:
                selected.append(boundary)
        selected.append(file_size)
        return [
            (left, right)
            for left, right in zip(selected, selected[1:])
            if right > left
        ]

    @staticmethod
    def _copy_file(source: Path, target: BinaryIO) -> None:
        with source.open("rb", buffering=_BUFFER_BYTES) as handle:
            shutil.copyfileobj(handle, target, length=_BUFFER_BYTES)

    @staticmethod
    def _format_family_id(number: int) -> bytes:
        return f"SMORF_FAM{number:09d}".encode("ascii")

    def _merge_family_shard(
        self,
        source: Path,
        target: BinaryIO,
        family_id_index: int,
        offset: int,
    ) -> None:
        with source.open("rb", buffering=_BUFFER_BYTES) as handle:
            buffer = bytearray()
            for line in handle:
                fields = line.rstrip(b"\r\n").split(b"\t")
                local_number = int(fields[family_id_index])
                fields[family_id_index] = self._format_family_id(
                    offset + local_number
                )
                buffer.extend(b"\t".join(fields) + b"\n")
                if len(buffer) >= _FLUSH_BYTES:
                    target.write(buffer)
                    buffer.clear()
            if buffer:
                target.write(buffer)

    def _merge_member_shard(
        self,
        source: Path,
        target: BinaryIO,
        offset: int,
    ) -> None:
        with source.open("rb", buffering=_BUFFER_BYTES) as handle:
            buffer = bytearray()
            for line in handle:
                fields = line.rstrip(b"\r\n").split(b"\t")
                fields[0] = self._format_family_id(offset + int(fields[0]))
                buffer.extend(b"\t".join(fields) + b"\n")
                if len(buffer) >= _FLUSH_BYTES:
                    target.write(buffer)
                    buffer.clear()
            if buffer:
                target.write(buffer)

    @staticmethod
    def _write_summary(path: Path, summary: ClusterSummary) -> None:
        with path.open("wb", buffering=_BUFFER_BYTES) as handle:
            handle.write(b"section\titem\tcount\n")
            overall = {
                "input_orfs": summary.input_orfs,
                "basic_passed": summary.basic_passed,
                "basic_removed": summary.basic_removed,
                "lnc_morf_removed": summary.lnc_morf_removed,
                "post_annotation_retained": (
                    summary.basic_passed - summary.lnc_morf_removed
                ),
                "exact_duplicates_collapsed": summary.exact_duplicates_collapsed,
                "alt_starts_collapsed": summary.alt_starts_collapsed,
                "total_collapsed": (
                    summary.exact_duplicates_collapsed
                    + summary.alt_starts_collapsed
                ),
                "family_count": summary.family_count,
                "singleton_families": summary.singleton_families,
                "multi_member_families": summary.multi_member_families,
                "genes_processed": summary.genes_processed,
                "effective_workers": summary.effective_workers,
            }
            for item, count in overall.items():
                handle.write(f"overall\t{item}\t{count}\n".encode("utf-8"))
            for reason, count in sorted(summary.removal_reasons.items()):
                handle.write(
                    f"removal_reason\t{reason}\t{count}\n".encode("utf-8")
                )
            for category, count in sorted(summary.primary_categories.items()):
                handle.write(
                    f"primary_category\t{category}\t{count}\n".encode("utf-8")
                )

    def cluster_file(
        self,
        input_path: str | Path,
        output_prefix: str | Path,
    ) -> ClusterSummary:
        """Cluster a scanner table with optional gene-safe multiprocessing."""
        input_file = Path(input_path).resolve()
        header = self.read_header(input_file)
        plan = _ColumnPlan.from_header(header)
        file_size = input_file.stat().st_size
        requested_workers = self._effective_workers(file_size)
        ranges = self._discover_ranges(input_file, plan, requested_workers)
        effective_workers = len(ranges)

        output_prefix_path = Path(output_prefix)
        output_prefix_path.parent.mkdir(parents=True, exist_ok=True)
        final_paths = {
            "family": Path(f"{output_prefix_path}.family.message.txt"),
            "members": Path(f"{output_prefix_path}.family.members.txt"),
            "removed": Path(f"{output_prefix_path}.removed.message.txt"),
            "summary": Path(f"{output_prefix_path}.cluster_summary.txt"),
        }
        temp_final = {
            name: path.with_name(path.name + ".tmp")
            for name, path in final_paths.items()
        }

        with tempfile.TemporaryDirectory(
            prefix=f".{output_prefix_path.name}.cluster.",
            dir=str(output_prefix_path.parent),
        ) as temp_directory:
            tasks = [
                _RangeTask(
                    shard_id=index,
                    start=start,
                    end=end,
                    temp_directory=temp_directory,
                )
                for index, (start, end) in enumerate(ranges)
            ]
            if effective_workers == 1:
                writer = _ShardWriter(Path(temp_directory), 0)
                engine = _ClusterEngine(
                    annotation=self.annotation,
                    plan=plan,
                    config=self.config,
                )
                with writer:
                    shard_summary = engine.process_range(
                        input_path=input_file,
                        start=ranges[0][0],
                        end=ranges[0][1],
                        writer=writer,
                    )
                results = [
                    _RangeResult(
                        shard_id=0,
                        paths=(
                            str(writer.paths["family"]),
                            str(writer.paths["members"]),
                            str(writer.paths["removed"]),
                        ),
                        summary=shard_summary,
                    )
                ]
            else:
                with ProcessPoolExecutor(
                    max_workers=effective_workers,
                    initializer=_worker_initialize,
                    initargs=(
                        self.annotation_path,
                        str(input_file),
                        header,
                        self.config,
                    ),
                ) as executor:
                    results = list(executor.map(_process_range_worker, tasks))

            results.sort(key=lambda item: item.shard_id)
            summary = ClusterSummary(effective_workers=effective_workers)
            family_id_index = plan.family_lookup[b"family_id"]
            family_offset = 0
            try:
                with temp_final["family"].open("wb", buffering=_BUFFER_BYTES) as family_handle, temp_final["members"].open("wb", buffering=_BUFFER_BYTES) as member_handle, temp_final["removed"].open("wb", buffering=_BUFFER_BYTES) as removed_handle:
                    family_handle.write(b"\t".join(plan.family_header) + b"\n")
                    member_handle.write(b"\t".join(MEMBER_COLUMNS) + b"\n")
                    removed_handle.write(b"\t".join(plan.filter_header) + b"\n")
                    for result in results:
                        family_path, member_path, removed_path = map(Path, result.paths)
                        self._merge_family_shard(
                            family_path,
                            family_handle,
                            family_id_index,
                            family_offset,
                        )
                        self._merge_member_shard(
                            member_path,
                            member_handle,
                            family_offset,
                        )
                        self._copy_file(removed_path, removed_handle)
                        family_offset += result.summary.family_count
                        summary.merge(result.summary)
                summary.effective_workers = effective_workers
                self._write_summary(temp_final["summary"], summary)
                for name, final_path in final_paths.items():
                    os.replace(temp_final[name], final_path)
            except Exception:
                for path in temp_final.values():
                    try:
                        path.unlink()
                    except FileNotFoundError:
                        pass
                raise

        return summary
