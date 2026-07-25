#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.004
# Function: Filter candidate smORFs with indexed binary streaming and optional parallel shards.
# Input: smorf_scanner ORF message table and an optional Kozak model.
# Output: Passed, removed, optional all-record, and summary tables.

"""High-performance structural filtering for candidate smORFs.

The file-level implementation is optimized for scanner-generated tab-delimited
records. It avoids per-row dictionaries, reconstructs output rows only once,
uses large binary buffers, and can process independent byte ranges in parallel.
The public in-memory APIs from the previous ``ORFFilter`` implementation remain
available for compatibility with other RiboParser modules.
"""

from __future__ import annotations

import math
import os
import shutil
import tempfile
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, BinaryIO

from utils.ribo.ArgsParser import progress_print

from .smorf_kozak import KozakModel, KozakResult

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
STOP_CODON_BYTES = frozenset({b"TAA", b"TAG", b"TGA"})
REQUIRED_COLUMNS = (
    "orf_id",
    "source_strand",
    "priority",
    "completeness",
    "category",
    "start_codon",
    "aa_length",
)
OUTPUT_COLUMNS = (
    "kozak_pwm_name",
    "kozak_pwm_score",
    "kozak_pwm_level",
    "kozak_valid_ratio",
    "structure_status",
    "filter_status",
    "filter_reason",
)
_BUFFER_BYTES = 8 * 1024 * 1024
_FLUSH_BYTES = 8 * 1024 * 1024
_MIN_PARALLEL_BYTES = 32 * 1024 * 1024
_KOZAK_CACHE_SIZE = 131_072
_DEFAULT_START_INDEX = 6
_DEFAULT_VALID_RATIO = 0.5
_BASE_BYTES = (65, 67, 71, 84)  # A, C, G, T


@dataclass(frozen=True, slots=True)
class FilterDecision:
    """Store one filtering decision."""

    keep: bool
    reasons: tuple[str, ...]
    structure_status: str


@dataclass(slots=True)
class FilterSummary:
    """Accumulate filtering statistics."""

    total: int = 0
    passed: int = 0
    removed: int = 0
    reason_counts: Counter[str] = field(default_factory=Counter)
    category_counts: Counter[str] = field(default_factory=Counter)
    workers: int = 1
    write_all: bool = False
    strict_block_check: bool = False
    kozak_scored: bool = False

    def update(
        self,
        record: Mapping[str, str],
        decision: FilterDecision,
    ) -> None:
        """Update counts for compatibility with in-memory callers."""
        self.total += 1
        category = str(record.get("category", "")).strip() or "NA"
        self.category_counts[category] += 1
        if decision.keep:
            self.passed += 1
        else:
            self.removed += 1
            self.reason_counts.update(decision.reasons)

    def merge(self, other: "_RangeSummary") -> None:
        """Merge one worker summary."""
        self.total += other.total
        self.passed += other.passed
        self.removed += other.removed
        self.reason_counts.update(other.reason_counts)
        self.category_counts.update(
            {
                category.decode("utf-8", errors="replace"): count
                for category, count in other.category_counts.items()
            }
        )


@dataclass(slots=True)
class _RangeSummary:
    """Store byte-range filtering counts."""

    total: int = 0
    passed: int = 0
    removed: int = 0
    reason_counts: Counter[str] = field(default_factory=Counter)
    category_counts: Counter[bytes] = field(default_factory=Counter)


@dataclass(frozen=True, slots=True)
class _KozakSpec:
    """Store the small, pickle-safe part of a Kozak model."""

    rows: tuple[tuple[float, float, float, float], ...]
    name: str
    start_index: int

    @classmethod
    def from_model(cls, model: KozakModel | None) -> "_KozakSpec | None":
        """Build a worker specification from a model."""
        if model is None:
            return None
        return cls(
            rows=tuple(tuple(float(value) for value in row) for row in model.rows),
            name=str(model.name),
            start_index=int(model.start_index),
        )


@dataclass(frozen=True, slots=True)
class _FilterConfig:
    """Store immutable row-filter settings for workers."""

    keep_start_codons: frozenset[bytes]
    min_aa: int
    max_aa: int
    keep_categories: frozenset[bytes]
    remove_categories: frozenset[bytes]
    require_sense: bool
    require_primary: bool
    require_complete: bool
    max_ambiguous_codons: int
    min_kozak_score: float
    validate_exon_blocks: bool
    write_all: bool
    kozak_spec: _KozakSpec | None


@dataclass(frozen=True, slots=True)
class _ColumnPlan:
    """Store integer positions of frequently accessed columns."""

    header: tuple[str, ...]
    output_header: tuple[str, ...]
    field_count: int
    output_field_count: int
    source_strand: int
    priority: int
    completeness: int
    category: int
    start_codon: int
    aa_length: int
    nt_length: int
    stop_codon: int
    ambiguous_codon_count: int
    exon_starts: int
    exon_ends: int
    exon_count: int
    kozak_seq: int
    kozak_start_index: int
    output_positions: tuple[int, ...]
    append_outputs: bool

    @classmethod
    def build(
        cls,
        header: Sequence[str],
        require_kozak: bool,
    ) -> "_ColumnPlan":
        """Compile a header into direct integer indexes."""
        lookup = {name: index for index, name in enumerate(header)}
        missing = [name for name in REQUIRED_COLUMNS if name not in lookup]
        if require_kozak and "kozak_seq" not in lookup:
            missing.append("kozak_seq")
        if missing:
            raise ValueError(
                "ORF message table is missing column(s): "
                + ", ".join(missing)
            )

        output_header = list(header)
        for name in OUTPUT_COLUMNS:
            if name not in lookup:
                lookup[name] = len(output_header)
                output_header.append(name)

        output_positions = tuple(lookup[name] for name in OUTPUT_COLUMNS)
        append_outputs = all(
            name not in header for name in OUTPUT_COLUMNS
        )
        return cls(
            header=tuple(header),
            output_header=tuple(output_header),
            field_count=len(header),
            output_field_count=len(output_header),
            source_strand=lookup["source_strand"],
            priority=lookup["priority"],
            completeness=lookup["completeness"],
            category=lookup["category"],
            start_codon=lookup["start_codon"],
            aa_length=lookup["aa_length"],
            nt_length=lookup.get("nt_length", -1),
            stop_codon=lookup.get("stop_codon", -1),
            ambiguous_codon_count=lookup.get(
                "ambiguous_codon_count",
                -1,
            ),
            exon_starts=lookup.get("exon_starts", -1),
            exon_ends=lookup.get("exon_ends", -1),
            exon_count=lookup.get("exon_count", -1),
            kozak_seq=lookup.get("kozak_seq", -1),
            kozak_start_index=lookup.get("kozak_start_index", -1),
            output_positions=output_positions,
            append_outputs=append_outputs,
        )


@dataclass(frozen=True, slots=True)
class _RangeTask:
    """Describe one independent file byte range."""

    input_path: str
    start: int
    end: int
    header: tuple[str, ...]
    config: _FilterConfig
    passed_path: str
    removed_path: str
    all_path: str | None
    write_header: bool
    report_progress: bool


class _FastKozakScorer:
    """Score short Kozak contexts without per-record alignment objects."""

    def __init__(self, spec: _KozakSpec) -> None:
        self.spec = spec
        self.name_bytes = spec.name.encode("utf-8")
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
                if probability <= 0:
                    value = 0.0
                elif denominator <= 0:
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
        """Return uppercase DNA bytes."""
        cleaned = value.strip().upper()
        return cleaned.replace(b"U", b"T") if b"U" in cleaned else cleaned

    @staticmethod
    def _parse_explicit(value: bytes) -> int | None:
        """Parse an explicit start index."""
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
        """Resolve a start anchor using scanner metadata first."""
        explicit = _FastKozakScorer._parse_explicit(provided)
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
    def score_raw(
        self,
        sequence_value: bytes,
        start_index_value: bytes,
        start_codon_value: bytes,
    ) -> tuple[float | None, bytes, float]:
        """Return score, level bytes, and valid ratio."""
        sequence = self._clean(sequence_value)
        start_codon = self._clean(start_codon_value)
        resolved = self._resolve_start_index(
            sequence=sequence,
            provided=start_index_value,
            start_codon=start_codon,
        )
        if resolved is None or self.informative_count == 0:
            return None, b"NA", 0.0

        offset = self.spec.start_index - resolved
        score_sum = 0.0
        valid_count = 0
        for target_position, score_table in self.informative:
            source_position = target_position - offset
            if not 0 <= source_position < len(sequence):
                continue
            value = score_table[sequence[source_position]]
            if value < 0:
                continue
            score_sum += value
            valid_count += 1

        valid_ratio = valid_count / self.informative_count
        if valid_count == 0 or valid_ratio < _DEFAULT_VALID_RATIO:
            return None, b"NA", valid_ratio

        score = score_sum / valid_count
        if score >= 0.75:
            level = b"strong"
        elif score >= 0.50:
            level = b"moderate"
        else:
            level = b"weak"
        return score, level, valid_ratio

    def score_result(
        self,
        sequence: str,
        start_index: Any,
        start_codon: str,
    ) -> KozakResult:
        """Return the public KozakResult representation."""
        score, level, valid_ratio = self.score_raw(
            str(sequence).encode("utf-8"),
            str(start_index if start_index is not None else "").encode(
                "utf-8"
            ),
            str(start_codon).encode("utf-8"),
        )
        return KozakResult(
            score=score,
            level=level.decode("ascii"),
            valid_ratio=valid_ratio,
        )


class ORFTable:
    """Read and write strict tab-delimited ORF tables."""

    @staticmethod
    def read_header(path: str | Path) -> list[str]:
        """Read and validate the input header."""
        input_path = Path(path)
        with input_path.open(
            "r",
            encoding="utf-8",
            buffering=_BUFFER_BYTES,
        ) as handle:
            line = handle.readline().rstrip("\n\r")
        if not line:
            raise ValueError(f"Empty ORF message table: {input_path}")
        header = line.split("\t")
        duplicates = [
            field
            for field, count in Counter(header).items()
            if count > 1
        ]
        if duplicates:
            raise ValueError(
                "Duplicate ORF column(s): "
                + ", ".join(sorted(duplicates))
            )
        return header

    @staticmethod
    def validate_header(
        header: Sequence[str],
        require_kozak: bool = False,
    ) -> None:
        """Validate fields required by the filtering engine."""
        _ColumnPlan.build(header, require_kozak=require_kozak)

    @staticmethod
    def iter_table(
        path: str | Path,
        header: Sequence[str] | None = None,
    ) -> Iterator[dict[str, str]]:
        """Yield dictionaries for backward-compatible callers."""
        input_path = Path(path)
        expected = list(header) if header is not None else None
        with input_path.open(
            "r",
            encoding="utf-8",
            buffering=_BUFFER_BYTES,
        ) as handle:
            actual = handle.readline().rstrip("\n\r").split("\t")
            if expected is None:
                expected = actual
            elif actual != expected:
                raise ValueError("ORF header changed between read passes.")
            expected_count = len(expected)
            for line_number, raw_line in enumerate(handle, start=2):
                line = raw_line.rstrip("\n\r")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) != expected_count:
                    raise ValueError(
                        "ORF field-count mismatch at line "
                        f"{line_number}: expected {expected_count}, "
                        f"observed {len(fields)}."
                    )
                yield dict(zip(expected, fields))

    @staticmethod
    def iter_kozak_training_records(
        path: str | Path,
        header: Sequence[str] | None = None,
    ) -> Iterator[dict[str, str]]:
        """Yield only columns needed to train an annotated Kozak model."""
        input_path = Path(path)
        expected = list(header) if header is not None else ORFTable.read_header(
            input_path
        )
        lookup = {name: index for index, name in enumerate(expected)}
        required = (
            "category",
            "source_strand",
            "priority",
            "completeness",
            "kozak_seq",
            "start_codon",
        )
        missing = [name for name in required if name not in lookup]
        if missing:
            raise ValueError(
                "ORF message table is missing Kozak column(s): "
                + ", ".join(missing)
            )
        selected_names = (
            "category",
            "source_strand",
            "priority",
            "completeness",
            "kozak_seq",
            "kozak_start_index",
            "start_codon",
            "gene_id",
        )
        selected_indices = tuple(lookup.get(name, -1) for name in selected_names)
        expected_count = len(expected)

        with input_path.open(
            "r",
            encoding="utf-8",
            buffering=_BUFFER_BYTES,
        ) as handle:
            actual = handle.readline().rstrip("\n\r").split("\t")
            if actual != expected:
                raise ValueError("ORF header changed between read passes.")
            for line_number, raw_line in enumerate(handle, start=2):
                line = raw_line.rstrip("\n\r")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) != expected_count:
                    raise ValueError(
                        "ORF field-count mismatch at line "
                        f"{line_number}: expected {expected_count}, "
                        f"observed {len(fields)}."
                    )
                yield {
                    name: fields[index] if index >= 0 else ""
                    for name, index in zip(selected_names, selected_indices)
                }

    @staticmethod
    def read_table(
        path: str | Path,
    ) -> tuple[list[str], list[dict[str, str]]]:
        """Read a complete table into memory."""
        header = ORFTable.read_header(path)
        ORFTable.validate_header(header)
        return header, list(ORFTable.iter_table(path, header))

    @staticmethod
    def write_table(
        path: str | Path,
        header: Sequence[str],
        records: Iterable[Mapping[str, Any]],
    ) -> None:
        """Write a complete ORF table."""
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open(
            "w",
            encoding="utf-8",
            buffering=_BUFFER_BYTES,
        ) as handle:
            handle.write("\t".join(header) + "\n")
            buffer: list[str] = []
            buffer_size = 0
            for record in records:
                line = (
                    "\t".join(str(record.get(field, "")) for field in header)
                    + "\n"
                )
                buffer.append(line)
                buffer_size += len(line)
                if buffer_size >= _FLUSH_BYTES:
                    handle.write("".join(buffer))
                    buffer.clear()
                    buffer_size = 0
            if buffer:
                handle.write("".join(buffer))


class ORFFilter:
    """Filter ORFs using metadata, structural QC, and optional Kozak score."""

    def __init__(
        self,
        keep_start_codons: str = "ATG,CTG,GTG,TTG",
        min_aa: int = 8,
        max_aa: int = 10000,
        keep_categories: str = (
            "uORF,dORF,lncORF,iORF,emORF,overlap_uORF,"
            "overlap_dORF,other_ORF,annotated_ORF,annotated_mORF"
        ),
        remove_categories: str = "same_frame_iORF,antisense_ORF",
        require_sense: bool = True,
        require_primary: bool = True,
        require_complete: bool = True,
        max_ambiguous_codons: int = 0,
        kozak_model: KozakModel | None = None,
        min_kozak_score: float = 0.0,
        validate_exon_blocks: bool = False,
        write_all: bool = False,
        threads: int = 0,
    ) -> None:
        self.keep_start_codons = {
            value.upper().replace("U", "T")
            for value in self._parse_set(keep_start_codons)
        }
        self.min_aa = int(min_aa)
        self.max_aa = int(max_aa)
        self.keep_categories = self._parse_set(keep_categories)
        self.remove_categories = self._parse_set(remove_categories)
        self.require_sense = bool(require_sense)
        self.require_primary = bool(require_primary)
        self.require_complete = bool(require_complete)
        self.max_ambiguous_codons = int(max_ambiguous_codons)
        self.kozak_model = kozak_model
        self.min_kozak_score = float(min_kozak_score)
        self.validate_exon_blocks = bool(validate_exon_blocks)
        self.write_all = bool(write_all)
        self.threads = int(threads)
        self._kozak_spec = _KozakSpec.from_model(kozak_model)
        self._fast_kozak = (
            _FastKozakScorer(self._kozak_spec)
            if self._kozak_spec is not None
            else None
        )

        if self.min_aa < 1:
            raise ValueError("min_aa must be >= 1.")
        if self.max_aa < self.min_aa:
            raise ValueError("max_aa must be >= min_aa.")
        if self.max_ambiguous_codons < 0:
            raise ValueError("max_ambiguous_codons must be >= 0.")
        if not 0 <= self.min_kozak_score <= 1:
            raise ValueError("min_kozak_score must be in [0, 1].")
        if self.min_kozak_score > 0 and self.kozak_model is None:
            raise ValueError(
                "A Kozak model is required when min_kozak_score > 0."
            )
        if self.threads < 0:
            raise ValueError("threads must be >= 0.")
        conflict = self.keep_categories.intersection(self.remove_categories)
        if conflict:
            raise ValueError(
                "Category cannot be both retained and removed: "
                + ", ".join(sorted(conflict))
            )

    @staticmethod
    def _parse_set(value: str | None) -> set[str]:
        """Parse comma-separated values."""
        if value is None:
            return set()
        return {
            item.strip()
            for item in str(value).split(",")
            if item.strip()
        }

    @staticmethod
    def output_header(header: Sequence[str]) -> list[str]:
        """Append missing filter columns."""
        output = list(header)
        for field in OUTPUT_COLUMNS:
            if field not in output:
                output.append(field)
        return output

    @staticmethod
    def _parse_int(value: Any) -> int | None:
        """Parse an integer without silent fallback."""
        try:
            text = str(value).strip()
            return int(text) if text else None
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _parse_int_list(value: Any) -> list[int] | None:
        """Parse comma-separated integers."""
        text = str(value).strip().rstrip(",")
        if not text:
            return []
        try:
            return [int(item) for item in text.split(",") if item]
        except ValueError:
            return None

    def annotate_kozak(
        self,
        record: dict[str, str],
    ) -> KozakResult | None:
        """Add compact Kozak annotations to one dictionary record."""
        if self._fast_kozak is None:
            record["kozak_pwm_name"] = "none"
            record["kozak_pwm_score"] = "NA"
            record["kozak_pwm_level"] = "NA"
            record["kozak_valid_ratio"] = "NA"
            return None
        result = self._fast_kozak.score_result(
            sequence=record.get("kozak_seq", ""),
            start_index=record.get("kozak_start_index"),
            start_codon=record.get("start_codon", ""),
        )
        record["kozak_pwm_name"] = self._kozak_spec.name
        record["kozak_pwm_score"] = (
            "NA" if result.score is None else f"{result.score:.6f}"
        )
        record["kozak_pwm_level"] = result.level
        record["kozak_valid_ratio"] = f"{result.valid_ratio:.6f}"
        return result

    def _structure_reasons(
        self,
        record: Mapping[str, str],
        aa_length: int | None,
        completeness: str,
    ) -> tuple[list[str], bool]:
        """Check enabled structural metadata."""
        reasons: list[str] = []
        checked = False
        nt_length: int | None = None

        if "nt_length" in record:
            checked = True
            nt_length = self._parse_int(record.get("nt_length"))
            if nt_length is None or nt_length <= 0:
                reasons.append("invalid_nt_length")
            elif aa_length is not None:
                expected = (
                    aa_length * 3 + 3
                    if completeness == "complete"
                    else aa_length * 3
                )
                if nt_length != expected:
                    reasons.append("inconsistent_orf_length")

        if "stop_codon" in record:
            checked = True
            stop_codon = (
                str(record.get("stop_codon", ""))
                .strip()
                .upper()
                .replace("U", "T")
            )
            if completeness == "complete" and stop_codon not in STOP_CODONS:
                reasons.append("invalid_stop_codon")

        if "ambiguous_codon_count" in record:
            checked = True
            ambiguous = self._parse_int(record.get("ambiguous_codon_count"))
            if ambiguous is None or ambiguous < 0:
                reasons.append("invalid_ambiguous_codon_count")
            elif ambiguous > self.max_ambiguous_codons:
                reasons.append("too_many_ambiguous_codons")

        if self.validate_exon_blocks:
            has_start = "exon_starts" in record
            has_end = "exon_ends" in record
            if has_start or has_end:
                checked = True
                if not (has_start and has_end):
                    reasons.append("incomplete_exon_blocks")
                else:
                    starts = self._parse_int_list(record.get("exon_starts"))
                    ends = self._parse_int_list(record.get("exon_ends"))
                    if starts is None or ends is None:
                        reasons.append("invalid_exon_blocks")
                    elif len(starts) != len(ends) or not starts:
                        reasons.append("invalid_exon_blocks")
                    elif any(
                        end <= start for start, end in zip(starts, ends)
                    ):
                        reasons.append("invalid_exon_blocks")
                    else:
                        if "exon_count" in record:
                            exon_count = self._parse_int(
                                record.get("exon_count")
                            )
                            if exon_count != len(starts):
                                reasons.append("inconsistent_exon_count")
                        if nt_length is not None and nt_length > 0:
                            block_length = sum(
                                end - start
                                for start, end in zip(starts, ends)
                            )
                            if block_length != nt_length:
                                reasons.append("inconsistent_block_length")
        return reasons, checked

    def decide(
        self,
        record: Mapping[str, str],
        kozak_result: KozakResult | None = None,
    ) -> FilterDecision:
        """Evaluate all enabled filters for an in-memory record."""
        reasons: list[str] = []
        source_strand = str(record.get("source_strand", "")).strip().lower()
        priority = str(record.get("priority", "")).strip()
        completeness = str(record.get("completeness", "")).strip()
        category = str(record.get("category", "")).strip()
        start_codon = (
            str(record.get("start_codon", ""))
            .strip()
            .upper()
            .replace("U", "T")
        )
        aa_length = self._parse_int(record.get("aa_length"))

        if source_strand not in {"sense", "antisense"}:
            reasons.append("invalid_source_strand")
        elif self.require_sense and source_strand != "sense":
            reasons.append("non_sense_strand")

        if not priority:
            reasons.append("missing_priority")
        elif self.require_primary and priority != "primary":
            reasons.append(f"non_primary:{priority}")

        if not completeness:
            reasons.append("missing_completeness")
        elif self.require_complete and completeness != "complete":
            reasons.append(f"incomplete_orf:{completeness}")

        if not category:
            reasons.append("missing_category")
        else:
            if category in self.remove_categories:
                reasons.append(f"removed_category:{category}")
            if self.keep_categories and category not in self.keep_categories:
                reasons.append(f"category_not_allowed:{category}")

        if not start_codon:
            reasons.append("missing_start_codon")
        elif self.keep_start_codons and start_codon not in self.keep_start_codons:
            reasons.append(f"start_codon_not_allowed:{start_codon}")

        if aa_length is None:
            reasons.append("invalid_aa_length")
        else:
            if aa_length < self.min_aa:
                reasons.append("too_short")
            if aa_length > self.max_aa:
                reasons.append("too_long")

        structure_reasons, structure_checked = self._structure_reasons(
            record=record,
            aa_length=aa_length,
            completeness=completeness,
        )
        reasons.extend(structure_reasons)

        if self.kozak_model is not None and self.min_kozak_score > 0:
            if kozak_result is None or kozak_result.score is None:
                reasons.append("insufficient_kozak_context")
            elif kozak_result.score < self.min_kozak_score:
                reasons.append("weak_kozak_pwm")

        structure_status = (
            "FAIL"
            if structure_reasons
            else ("PASS" if structure_checked else "NA")
        )
        return FilterDecision(
            keep=not reasons,
            reasons=tuple(reasons),
            structure_status=structure_status,
        )

    def check_record(
        self,
        record: Mapping[str, str],
    ) -> tuple[bool, str]:
        """Return a backward-compatible flag and reason string."""
        decision = self.decide(record)
        return (
            decision.keep,
            "PASS" if decision.keep else ";".join(decision.reasons),
        )

    @staticmethod
    def _annotate_decision(
        record: dict[str, str],
        decision: FilterDecision,
    ) -> None:
        """Add filter decision fields."""
        record["structure_status"] = decision.structure_status
        record["filter_status"] = "PASS" if decision.keep else "FAIL"
        record["filter_reason"] = (
            "PASS" if decision.keep else ";".join(decision.reasons)
        )

    def filter_records(
        self,
        records: list[dict[str, str]],
        report_every: int = 0,
    ) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
        """Filter a mutable in-memory record list."""
        passed: list[dict[str, str]] = []
        removed: list[dict[str, str]] = []
        total = len(records)
        for index, record in enumerate(records, start=1):
            if report_every > 0 and index % report_every == 0:
                progress_print(f"records={index:,}/{total:,}.")
            kozak_result = self.annotate_kozak(record)
            decision = self.decide(record, kozak_result)
            self._annotate_decision(record, decision)
            (passed if decision.keep else removed).append(record)
        return passed, removed

    def add_kozak_pwm_fields(
        self,
        header: list[str],
        records: list[dict[str, str]],
    ) -> list[str]:
        """Annotate Kozak fields for compatibility callers."""
        for record in records:
            self.annotate_kozak(record)
        return self.output_header(header)

    def _config(self) -> _FilterConfig:
        """Build the immutable binary-engine configuration."""
        return _FilterConfig(
            keep_start_codons=frozenset(
                codon.encode("ascii") for codon in self.keep_start_codons
            ),
            min_aa=self.min_aa,
            max_aa=self.max_aa,
            keep_categories=frozenset(
                category.encode("utf-8")
                for category in self.keep_categories
            ),
            remove_categories=frozenset(
                category.encode("utf-8")
                for category in self.remove_categories
            ),
            require_sense=self.require_sense,
            require_primary=self.require_primary,
            require_complete=self.require_complete,
            max_ambiguous_codons=self.max_ambiguous_codons,
            min_kozak_score=self.min_kozak_score,
            validate_exon_blocks=self.validate_exon_blocks,
            write_all=self.write_all,
            kozak_spec=self._kozak_spec,
        )

    @staticmethod
    def _available_cpus() -> int:
        """Return CPUs available to the current process or SLURM allocation."""
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except (AttributeError, OSError):
            return max(1, os.cpu_count() or 1)

    def _requested_workers(self) -> int:
        """Resolve automatic worker selection."""
        if self.threads > 0:
            return self.threads
        if self.kozak_model is not None or self.validate_exon_blocks:
            return self._available_cpus()
        return 1

    @staticmethod
    def _split_ranges(
        path: Path,
        data_start: int,
        requested_workers: int,
    ) -> list[tuple[int, int]]:
        """Split a text file into newline-aligned byte ranges."""
        file_size = path.stat().st_size
        data_size = max(0, file_size - data_start)
        if (
            requested_workers <= 1
            or data_size < _MIN_PARALLEL_BYTES
        ):
            return [(data_start, file_size)]

        workers = min(
            requested_workers,
            max(1, data_size // _MIN_PARALLEL_BYTES),
        )
        if workers <= 1:
            return [(data_start, file_size)]

        boundaries = [data_start]
        with path.open("rb", buffering=_BUFFER_BYTES) as handle:
            for index in range(1, workers):
                target = data_start + data_size * index // workers
                handle.seek(target)
                handle.readline()
                boundary = handle.tell()
                if boundary > boundaries[-1] and boundary < file_size:
                    boundaries.append(boundary)
        boundaries.append(file_size)
        return [
            (left, right)
            for left, right in zip(boundaries, boundaries[1:])
            if right > left
        ]

    @staticmethod
    def _write_summary_file(
        path: Path,
        summary: FilterSummary,
    ) -> None:
        """Write the compact filtering summary."""
        with path.open(
            "w",
            encoding="utf-8",
            buffering=_BUFFER_BYTES,
        ) as handle:
            handle.write("section\titem\tcount\n")
            handle.write(f"overall\ttotal\t{summary.total}\n")
            handle.write(f"overall\tpassed\t{summary.passed}\n")
            handle.write(f"overall\tremoved\t{summary.removed}\n")
            handle.write(f"runtime\tworkers\t{summary.workers}\n")
            handle.write(
                "runtime\twrite_all\t"
                f"{int(summary.write_all)}\n"
            )
            handle.write(
                "runtime\tstrict_block_check\t"
                f"{int(summary.strict_block_check)}\n"
            )
            handle.write(
                "runtime\tkozak_scored\t"
                f"{int(summary.kozak_scored)}\n"
            )
            for category in sorted(summary.category_counts):
                handle.write(
                    f"category\t{category}\t"
                    f"{summary.category_counts[category]}\n"
                )
            for reason in sorted(summary.reason_counts):
                handle.write(
                    f"reason\t{reason}\t"
                    f"{summary.reason_counts[reason]}\n"
                )

    @staticmethod
    def _merge_shards(
        output_path: Path,
        header_line: bytes,
        shard_paths: Sequence[Path],
    ) -> None:
        """Merge ordered worker shards and prepend one header."""
        with output_path.open("wb", buffering=_BUFFER_BYTES) as output:
            output.write(header_line)
            for shard_path in shard_paths:
                with shard_path.open("rb", buffering=_BUFFER_BYTES) as source:
                    shutil.copyfileobj(source, output, length=_BUFFER_BYTES)

    def filter_file(
        self,
        input_path: str | Path,
        output_prefix: str | Path,
    ) -> FilterSummary:
        """Filter a scanner message table with the binary streaming engine."""
        source_path = Path(input_path)
        header = ORFTable.read_header(source_path)
        plan = _ColumnPlan.build(
            header,
            require_kozak=self.kozak_model is not None,
        )
        config = self._config()
        output_prefix_path = Path(output_prefix)
        output_prefix_path.parent.mkdir(parents=True, exist_ok=True)
        requested_workers = self._requested_workers()

        with source_path.open("rb", buffering=_BUFFER_BYTES) as source:
            raw_header = source.readline()
            data_start = source.tell()
        if not raw_header:
            raise ValueError(f"Empty ORF message table: {source_path}")

        ranges = self._split_ranges(
            path=source_path,
            data_start=data_start,
            requested_workers=requested_workers,
        )
        summary = FilterSummary(
            workers=len(ranges),
            write_all=self.write_all,
            strict_block_check=self.validate_exon_blocks,
            kozak_scored=self.kozak_model is not None,
        )
        output_header_line = (
            "\t".join(plan.output_header).encode("utf-8") + b"\n"
        )

        final_paths = {
            "passed": Path(f"{output_prefix_path}.passed.message.txt"),
            "removed": Path(f"{output_prefix_path}.removed.message.txt"),
            "summary": Path(f"{output_prefix_path}.filter_summary.txt"),
        }
        if self.write_all:
            final_paths["all"] = Path(
                f"{output_prefix_path}.all.message.txt"
            )

        with tempfile.TemporaryDirectory(
            prefix=f".{output_prefix_path.name}.filter.",
            dir=str(output_prefix_path.parent),
        ) as temporary_directory:
            temporary_root = Path(temporary_directory)
            if len(ranges) == 1:
                task = _RangeTask(
                    input_path=str(source_path),
                    start=ranges[0][0],
                    end=ranges[0][1],
                    header=tuple(header),
                    config=config,
                    passed_path=str(temporary_root / "passed.complete"),
                    removed_path=str(temporary_root / "removed.complete"),
                    all_path=(
                        str(temporary_root / "all.complete")
                        if self.write_all
                        else None
                    ),
                    write_header=True,
                    report_progress=True,
                )
                range_summary = _filter_range_worker(task)
                summary.merge(range_summary)
                temporary_outputs = {
                    "passed": Path(task.passed_path),
                    "removed": Path(task.removed_path),
                }
                if task.all_path is not None:
                    temporary_outputs["all"] = Path(task.all_path)
            else:
                tasks: list[_RangeTask] = []
                for index, (start, end) in enumerate(ranges):
                    tasks.append(
                        _RangeTask(
                            input_path=str(source_path),
                            start=start,
                            end=end,
                            header=tuple(header),
                            config=config,
                            passed_path=str(
                                temporary_root / f"passed.{index:04d}"
                            ),
                            removed_path=str(
                                temporary_root / f"removed.{index:04d}"
                            ),
                            all_path=(
                                str(temporary_root / f"all.{index:04d}")
                                if self.write_all
                                else None
                            ),
                            write_header=False,
                            report_progress=False,
                        )
                    )

                results: dict[int, _RangeSummary] = {}
                with ProcessPoolExecutor(
                    max_workers=len(tasks),
                ) as executor:
                    future_to_index = {
                        executor.submit(_filter_range_worker, task): index
                        for index, task in enumerate(tasks)
                    }
                    completed_total = 0
                    for future in as_completed(future_to_index):
                        index = future_to_index[future]
                        result = future.result()
                        results[index] = result
                        completed_total += result.total
                        progress_print(
                            "completed_shards={done:,}/{total:,}, "
                            "records={records:,}.".format(
                                done=len(results),
                                total=len(tasks),
                                records=completed_total,
                            )
                        )

                for index in range(len(tasks)):
                    summary.merge(results[index])

                temporary_outputs = {}
                passed_complete = temporary_root / "passed.complete"
                self._merge_shards(
                    output_path=passed_complete,
                    header_line=output_header_line,
                    shard_paths=[Path(task.passed_path) for task in tasks],
                )
                temporary_outputs["passed"] = passed_complete

                removed_complete = temporary_root / "removed.complete"
                self._merge_shards(
                    output_path=removed_complete,
                    header_line=output_header_line,
                    shard_paths=[Path(task.removed_path) for task in tasks],
                )
                temporary_outputs["removed"] = removed_complete

                if self.write_all:
                    all_complete = temporary_root / "all.complete"
                    self._merge_shards(
                        output_path=all_complete,
                        header_line=output_header_line,
                        shard_paths=[
                            Path(task.all_path)
                            for task in tasks
                            if task.all_path is not None
                        ],
                    )
                    temporary_outputs["all"] = all_complete

            summary_complete = temporary_root / "summary.complete"
            self._write_summary_file(summary_complete, summary)
            temporary_outputs["summary"] = summary_complete

            for name, final_path in final_paths.items():
                os.replace(temporary_outputs[name], final_path)

        return summary

    read_table = staticmethod(ORFTable.read_table)
    write_table = staticmethod(ORFTable.write_table)


def _parse_int_bytes(value: bytes) -> int | None:
    """Parse integer bytes."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _normalize_codon_bytes(value: bytes) -> bytes:
    """Return uppercase DNA codon bytes."""
    codon = value.strip().upper()
    return codon.replace(b"U", b"T") if b"U" in codon else codon


def _validate_blocks_bytes(
    fields: list[bytes],
    plan: _ColumnPlan,
    nt_length: int | None,
) -> list[str]:
    """Strictly validate exon blocks when explicitly requested."""
    reasons: list[str] = []
    has_start = plan.exon_starts >= 0
    has_end = plan.exon_ends >= 0
    if not (has_start or has_end):
        return reasons
    if not (has_start and has_end):
        return ["incomplete_exon_blocks"]

    starts_text = fields[plan.exon_starts].strip().rstrip(b",")
    ends_text = fields[plan.exon_ends].strip().rstrip(b",")
    if not starts_text or not ends_text:
        return ["invalid_exon_blocks"]
    start_fields = starts_text.split(b",")
    end_fields = ends_text.split(b",")
    if len(start_fields) != len(end_fields) or not start_fields:
        return ["invalid_exon_blocks"]

    block_length = 0
    try:
        for start_value, end_value in zip(start_fields, end_fields):
            start = int(start_value)
            end = int(end_value)
            if end <= start:
                return ["invalid_exon_blocks"]
            block_length += end - start
    except ValueError:
        return ["invalid_exon_blocks"]

    if plan.exon_count >= 0:
        exon_count = _parse_int_bytes(fields[plan.exon_count])
        if exon_count != len(start_fields):
            reasons.append("inconsistent_exon_count")
    if nt_length is not None and nt_length > 0 and block_length != nt_length:
        reasons.append("inconsistent_block_length")
    return reasons


def _format_output_row(
    original_line: bytes,
    fields: list[bytes],
    plan: _ColumnPlan,
    output_values: tuple[bytes, ...],
) -> bytes:
    """Append or replace filter columns and return one encoded row."""
    if plan.append_outputs:
        return original_line + b"\t" + b"\t".join(output_values) + b"\n"

    if len(fields) < plan.output_field_count:
        fields.extend([b""] * (plan.output_field_count - len(fields)))
    for position, value in zip(plan.output_positions, output_values):
        fields[position] = value
    return b"\t".join(fields) + b"\n"


def _evaluate_binary_record(
    fields: list[bytes],
    plan: _ColumnPlan,
    config: _FilterConfig,
    kozak_scorer: _FastKozakScorer | None,
) -> tuple[bool, list[str], bytes, tuple[bytes, ...], bytes]:
    """Evaluate one split binary row without constructing a dictionary."""
    reasons: list[str] = []
    structure_reasons: list[str] = []
    structure_checked = False

    source_strand = fields[plan.source_strand].strip().lower()
    priority = fields[plan.priority].strip()
    completeness = fields[plan.completeness].strip()
    category = fields[plan.category].strip()
    start_codon = _normalize_codon_bytes(fields[plan.start_codon])
    aa_length = _parse_int_bytes(fields[plan.aa_length])

    if source_strand not in {b"sense", b"antisense"}:
        reasons.append("invalid_source_strand")
    elif config.require_sense and source_strand != b"sense":
        reasons.append("non_sense_strand")

    if not priority:
        reasons.append("missing_priority")
    elif config.require_primary and priority != b"primary":
        reasons.append(
            "non_primary:" + priority.decode("utf-8", errors="replace")
        )

    if not completeness:
        reasons.append("missing_completeness")
    elif config.require_complete and completeness != b"complete":
        reasons.append(
            "incomplete_orf:"
            + completeness.decode("utf-8", errors="replace")
        )

    if not category:
        reasons.append("missing_category")
    else:
        category_text = category.decode("utf-8", errors="replace")
        if category in config.remove_categories:
            reasons.append("removed_category:" + category_text)
        if config.keep_categories and category not in config.keep_categories:
            reasons.append("category_not_allowed:" + category_text)

    if not start_codon:
        reasons.append("missing_start_codon")
    elif config.keep_start_codons and start_codon not in config.keep_start_codons:
        reasons.append(
            "start_codon_not_allowed:"
            + start_codon.decode("ascii", errors="replace")
        )

    if aa_length is None:
        reasons.append("invalid_aa_length")
    else:
        if aa_length < config.min_aa:
            reasons.append("too_short")
        if aa_length > config.max_aa:
            reasons.append("too_long")

    nt_length: int | None = None
    if plan.nt_length >= 0:
        structure_checked = True
        nt_length = _parse_int_bytes(fields[plan.nt_length])
        if nt_length is None or nt_length <= 0:
            structure_reasons.append("invalid_nt_length")
        elif aa_length is not None:
            expected = (
                aa_length * 3 + 3
                if completeness == b"complete"
                else aa_length * 3
            )
            if nt_length != expected:
                structure_reasons.append("inconsistent_orf_length")

    if plan.stop_codon >= 0:
        structure_checked = True
        stop_codon = _normalize_codon_bytes(fields[plan.stop_codon])
        if completeness == b"complete" and stop_codon not in STOP_CODON_BYTES:
            structure_reasons.append("invalid_stop_codon")

    if plan.ambiguous_codon_count >= 0:
        structure_checked = True
        ambiguous = _parse_int_bytes(fields[plan.ambiguous_codon_count])
        if ambiguous is None or ambiguous < 0:
            structure_reasons.append("invalid_ambiguous_codon_count")
        elif ambiguous > config.max_ambiguous_codons:
            structure_reasons.append("too_many_ambiguous_codons")

    if config.validate_exon_blocks:
        if plan.exon_starts >= 0 or plan.exon_ends >= 0:
            structure_checked = True
        structure_reasons.extend(
            _validate_blocks_bytes(
                fields=fields,
                plan=plan,
                nt_length=nt_length,
            )
        )
    reasons.extend(structure_reasons)

    if kozak_scorer is None:
        kozak_name = b"none"
        kozak_score_text = b"NA"
        kozak_level = b"NA"
        valid_ratio_text = b"NA"
        kozak_score = None
    else:
        sequence = fields[plan.kozak_seq]
        start_index = (
            fields[plan.kozak_start_index]
            if plan.kozak_start_index >= 0
            else b""
        )
        kozak_score, kozak_level, valid_ratio = kozak_scorer.score_raw(
            sequence,
            start_index,
            start_codon,
        )
        kozak_name = kozak_scorer.name_bytes
        kozak_score_text = (
            b"NA"
            if kozak_score is None
            else f"{kozak_score:.6f}".encode("ascii")
        )
        valid_ratio_text = f"{valid_ratio:.6f}".encode("ascii")
        if config.min_kozak_score > 0:
            if kozak_score is None:
                reasons.append("insufficient_kozak_context")
            elif kozak_score < config.min_kozak_score:
                reasons.append("weak_kozak_pwm")

    structure_status = (
        b"FAIL"
        if structure_reasons
        else (b"PASS" if structure_checked else b"NA")
    )
    keep = not reasons
    filter_status = b"PASS" if keep else b"FAIL"
    reason_text = (
        b"PASS"
        if keep
        else ";".join(reasons).encode("utf-8")
    )
    output_values = (
        kozak_name,
        kozak_score_text,
        kozak_level,
        valid_ratio_text,
        structure_status,
        filter_status,
        reason_text,
    )
    return keep, reasons, structure_status, output_values, category or b"NA"


def _flush_buffer(
    handle: BinaryIO,
    buffer: list[bytes],
) -> None:
    """Write and clear one line buffer."""
    if buffer:
        handle.write(b"".join(buffer))
        buffer.clear()


def _filter_range_worker(task: _RangeTask) -> _RangeSummary:
    """Filter one newline-aligned byte range and write independent shards."""
    plan = _ColumnPlan.build(
        task.header,
        require_kozak=task.config.kozak_spec is not None,
    )
    kozak_scorer = (
        _FastKozakScorer(task.config.kozak_spec)
        if task.config.kozak_spec is not None
        else None
    )
    header_line = "\t".join(plan.output_header).encode("utf-8") + b"\n"
    summary = _RangeSummary()

    handles: dict[str, BinaryIO] = {}
    try:
        handles["passed"] = open(
            task.passed_path,
            "wb",
            buffering=_BUFFER_BYTES,
        )
        handles["removed"] = open(
            task.removed_path,
            "wb",
            buffering=_BUFFER_BYTES,
        )
        if task.all_path is not None:
            handles["all"] = open(
                task.all_path,
                "wb",
                buffering=_BUFFER_BYTES,
            )
        if task.write_header:
            for handle in handles.values():
                handle.write(header_line)

        with open(task.input_path, "rb", buffering=_BUFFER_BYTES) as source:
            source.seek(task.start)
            position = task.start
            while position < task.end:
                raw_line = source.readline()
                if not raw_line:
                    break
                line_start = position
                position += len(raw_line)
                if line_start >= task.end:
                    break
                original_line = raw_line.rstrip(b"\r\n")
                if not original_line:
                    continue
                fields = original_line.split(b"\t")
                if len(fields) != plan.field_count:
                    raise ValueError(
                        "ORF field-count mismatch near byte offset "
                        f"{line_start}: expected {plan.field_count}, "
                        f"observed {len(fields)}."
                    )

                keep, reasons, _structure_status, output_values, category = (
                    _evaluate_binary_record(
                        fields=fields,
                        plan=plan,
                        config=task.config,
                        kozak_scorer=kozak_scorer,
                    )
                )
                output_row = _format_output_row(
                    original_line=original_line,
                    fields=fields,
                    plan=plan,
                    output_values=output_values,
                )
                target = "passed" if keep else "removed"
                handles[target].write(output_row)

                if "all" in handles:
                    handles["all"].write(output_row)

                summary.total += 1
                summary.category_counts[category] += 1
                if keep:
                    summary.passed += 1
                else:
                    summary.removed += 1
                    summary.reason_counts.update(reasons)

                if (
                    task.report_progress
                    and summary.total % 1_000_000 == 0
                ):
                    progress_print(
                        "records={total:,}, passed={passed:,}, "
                        "removed={removed:,}.".format(
                            total=summary.total,
                            passed=summary.passed,
                            removed=summary.removed,
                        )
                    )

    finally:
        for handle in handles.values():
            handle.close()
    return summary


# Backward-compatible aliases.
FastORFFilter = ORFFilter
read_table = ORFTable.read_table
write_table = ORFTable.write_table
