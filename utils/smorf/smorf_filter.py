#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Perform structural quality control and optional Kozak filtering of smORFs.
# Input: ORF message table and an optional Kozak model.
# Output: Passed, removed, all-record, and summary tables.

"""Streaming structural and rule-based filtering for candidate smORFs.

The filter has two responsibilities:

1. Basic ORF quality control using scanner metadata.
2. Optional Kozak flanking-context annotation and thresholding.

Ribo-seq evidence, proteomics, conservation, and cross-method consensus are
intentionally kept outside this module.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType
from typing import Any, TextIO

from utils.ribo.ArgsParser import progress_print
from .smorf_kozak import KozakModel, KozakResult


STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
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


@dataclass(frozen=True, slots=True)
class FilterDecision:
    """Store one filtering decision.

    Attributes:
        keep: Whether the ORF passes all enabled rules.
        reasons: Ordered failure reasons.
        structure_status: ``PASS``, ``FAIL``, or ``NA``.
    """

    keep: bool
    reasons: tuple[str, ...]
    structure_status: str


@dataclass(slots=True)
class FilterSummary:
    """Accumulate filtering statistics.

    Attributes:
        total: Number of processed ORFs.
        passed: Number of passed ORFs.
        removed: Number of removed ORFs.
        reason_counts: Counts by failure reason.
        category_counts: Counts by ORF category.
    """

    total: int = 0
    passed: int = 0
    removed: int = 0
    reason_counts: Counter[str] = field(default_factory=Counter)
    category_counts: Counter[str] = field(default_factory=Counter)

    def update(
        self,
        record: Mapping[str, str],
        decision: FilterDecision,
    ) -> None:
        """Update summary counts.

        Args:
            record: Annotated ORF record.
            decision: Filtering decision.
        """
        self.total += 1
        self.category_counts[
            str(record.get("category", "")).strip() or "NA"
        ] += 1

        if decision.keep:
            self.passed += 1
        else:
            self.removed += 1
            self.reason_counts.update(decision.reasons)


class ORFTable:
    """Read and write strict tab-delimited ORF tables."""

    @staticmethod
    def read_header(path: str | Path) -> list[str]:
        """Read and validate the input header.

        Args:
            path: ORF message table.

        Returns:
            Header fields.

        Raises:
            ValueError: If the header is empty or duplicated.
        """
        input_path = Path(path)
        with input_path.open("r", encoding="utf-8") as handle:
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
        """Validate required columns.

        Args:
            header: Input header.
            require_kozak: Require ``kozak_seq``.

        Raises:
            ValueError: If a required column is missing.
        """
        required = list(REQUIRED_COLUMNS)
        if require_kozak:
            required.append("kozak_seq")

        missing = [field for field in required if field not in header]
        if missing:
            raise ValueError(
                "ORF message table is missing column(s): "
                + ", ".join(missing)
            )

    @staticmethod
    def iter_table(
        path: str | Path,
        header: Sequence[str] | None = None,
    ) -> Iterator[dict[str, str]]:
        """Yield strictly parsed records.

        Args:
            path: Input table.
            header: Optional previously read header.

        Yields:
            ORF records.

        Raises:
            ValueError: If a row has the wrong field count.
        """
        input_path = Path(path)
        expected = list(header) if header is not None else None

        with input_path.open("r", encoding="utf-8") as handle:
            actual = handle.readline().rstrip("\n\r").split("\t")
            if expected is None:
                expected = actual
            elif actual != expected:
                raise ValueError("ORF header changed between read passes.")

            for line_number, raw_line in enumerate(handle, start=2):
                line = raw_line.rstrip("\n\r")
                if not line:
                    continue

                fields = line.split("\t")
                if len(fields) != len(expected):
                    raise ValueError(
                        f"ORF field-count mismatch at line {line_number}: "
                        f"expected {len(expected)}, observed {len(fields)}."
                    )
                yield dict(zip(expected, fields))

    @staticmethod
    def read_table(
        path: str | Path,
    ) -> tuple[list[str], list[dict[str, str]]]:
        """Read a complete table into memory.

        Args:
            path: Input table.

        Returns:
            Header and records.
        """
        header = ORFTable.read_header(path)
        ORFTable.validate_header(header)
        return header, list(ORFTable.iter_table(path, header))

    @staticmethod
    def write_table(
        path: str | Path,
        header: Sequence[str],
        records: Iterable[Mapping[str, Any]],
    ) -> None:
        """Write an ORF table.

        Args:
            path: Output table.
            header: Output fields.
            records: ORF records.
        """
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w", encoding="utf-8") as handle:
            handle.write("\t".join(header) + "\n")
            for record in records:
                handle.write(
                    "\t".join(
                        str(record.get(field, ""))
                        for field in header
                    )
                    + "\n"
                )


class _AtomicWriter:
    """Atomically stream filter outputs."""

    def __init__(
        self,
        prefix: str | Path,
        header: Sequence[str],
    ) -> None:
        output_prefix = Path(prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)

        self.header = list(header)
        self.final_paths = {
            "passed": Path(f"{output_prefix}.passed.message.txt"),
            "removed": Path(f"{output_prefix}.removed.message.txt"),
            "all": Path(f"{output_prefix}.all.message.txt"),
            "summary": Path(f"{output_prefix}.filter_summary.txt"),
        }
        self.temp_paths = {
            name: path.with_name(path.name + ".tmp")
            for name, path in self.final_paths.items()
        }
        self.handles: dict[str, TextIO] = {}

    def __enter__(self) -> "_AtomicWriter":
        """Open temporary output files."""
        self.handles = {
            name: path.open("w", encoding="utf-8")
            for name, path in self.temp_paths.items()
        }
        header_line = "\t".join(self.header) + "\n"
        for name in ("passed", "removed", "all"):
            self.handles[name].write(header_line)
        return self

    def write(
        self,
        target: str,
        record: Mapping[str, Any],
    ) -> None:
        """Write one record.

        Args:
            target: ``passed``, ``removed``, or ``all``.
            record: ORF record.
        """
        self.handles[target].write(
            "\t".join(
                str(record.get(field, ""))
                for field in self.header
            )
            + "\n"
        )

    def write_summary(self, summary: FilterSummary) -> None:
        """Write overall, category, and reason counts.

        Args:
            summary: Filter summary.
        """
        handle = self.handles["summary"]
        handle.write("section\titem\tcount\n")
        handle.write(f"overall\ttotal\t{summary.total}\n")
        handle.write(f"overall\tpassed\t{summary.passed}\n")
        handle.write(f"overall\tremoved\t{summary.removed}\n")

        for category in sorted(summary.category_counts):
            handle.write(
                f"category\t{category}\t"
                f"{summary.category_counts[category]}\n"
            )
        for reason in sorted(summary.reason_counts):
            handle.write(
                f"reason\t{reason}\t{summary.reason_counts[reason]}\n"
            )

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        """Commit successful outputs or discard temporary files."""
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()

        if exception_type is None:
            for name, temp_path in self.temp_paths.items():
                temp_path.replace(self.final_paths[name])
        else:
            for temp_path in self.temp_paths.values():
                temp_path.unlink(missing_ok=True)
        return False


class ORFFilter:
    """Filter ORFs using metadata, structural QC, and optional Kozak score.

    Args:
        keep_start_codons: Comma-separated retained start codons.
        min_aa: Minimum peptide length.
        max_aa: Maximum peptide length.
        keep_categories: Comma-separated category allowlist.
        remove_categories: Comma-separated category denylist.
        require_sense: Require sense-strand scanning orientation.
        require_primary: Require primary overlap priority.
        require_complete: Require a complete terminal stop.
        max_ambiguous_codons: Maximum accepted ambiguous codon count.
        kozak_model: Optional Kozak model.
        min_kozak_score: Minimum normalized Kozak score.

    Raises:
        ValueError: If filter settings are inconsistent.
    """

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

        conflict = self.keep_categories.intersection(
            self.remove_categories
        )
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
        """Append filter output columns.

        Args:
            header: Input fields.

        Returns:
            Extended output fields.
        """
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
        """Add compact Kozak annotations.

        Args:
            record: Mutable ORF record.

        Returns:
            Kozak result or ``None``.
        """
        if self.kozak_model is None:
            record["kozak_pwm_name"] = "none"
            record["kozak_pwm_score"] = "NA"
            record["kozak_pwm_level"] = "NA"
            record["kozak_valid_ratio"] = "NA"
            return None

        result = self.kozak_model.score(
            sequence=record.get("kozak_seq", ""),
            start_index=record.get("kozak_start_index"),
            start_codon=record.get("start_codon", ""),
        )
        record["kozak_pwm_name"] = self.kozak_model.name
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
        """Check optional ORF structural metadata.

        Args:
            record: ORF record.
            aa_length: Parsed peptide length.
            completeness: ORF completeness label.

        Returns:
            Structural failure reasons and whether any structural field was
            available.
        """
        reasons: list[str] = []
        checked = False

        nt_length = None
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
            ambiguous = self._parse_int(
                record.get("ambiguous_codon_count")
            )
            if ambiguous is None or ambiguous < 0:
                reasons.append("invalid_ambiguous_codon_count")
            elif ambiguous > self.max_ambiguous_codons:
                reasons.append("too_many_ambiguous_codons")

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
                elif any(end <= start for start, end in zip(starts, ends)):
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
        """Evaluate all enabled filters.

        Args:
            record: Annotated ORF record.
            kozak_result: Optional previously calculated Kozak result.

        Returns:
            Filtering decision.
        """
        reasons: list[str] = []

        source_strand = str(
            record.get("source_strand", "")
        ).strip().lower()
        priority = str(record.get("priority", "")).strip()
        completeness = str(
            record.get("completeness", "")
        ).strip()
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
            if (
                self.keep_categories
                and category not in self.keep_categories
            ):
                reasons.append(f"category_not_allowed:{category}")

        if not start_codon:
            reasons.append("missing_start_codon")
        elif (
            self.keep_start_codons
            and start_codon not in self.keep_start_codons
        ):
            reasons.append(
                f"start_codon_not_allowed:{start_codon}"
            )

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
        """Return a backward-compatible pass flag and reason string.

        Args:
            record: ORF record. Kozak fields should already be annotated when
                a Kozak threshold is enabled.

        Returns:
            Pass flag and ``PASS`` or semicolon-delimited failure reasons.
        """
        decision = self.decide(record)
        return (
            decision.keep,
            "PASS" if decision.keep else ";".join(decision.reasons),
        )

    def filter_records(
        self,
        records: list[dict[str, str]],
        report_every: int = 0,
    ) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
        """Filter an in-memory record list.

        Args:
            records: Mutable ORF records.
            report_every: Optional compatibility progress interval.

        Returns:
            Passed and removed records.
        """
        passed = []
        removed = []
        total = len(records)

        for index, record in enumerate(records, start=1):
            if report_every > 0 and index % report_every == 0:
                progress_print(
                    "records={done:,}/{total:,}.".format(
                        done=index,
                        total=total,
                    )
                )
            kozak_result = self.annotate_kozak(record)
            decision = self.decide(record, kozak_result)
            self._annotate_decision(record, decision)
            (passed if decision.keep else removed).append(record)

        return passed, removed

    @staticmethod
    def _annotate_decision(
        record: dict[str, str],
        decision: FilterDecision,
    ) -> None:
        """Add decision fields to one record."""
        record["structure_status"] = decision.structure_status
        record["filter_status"] = "PASS" if decision.keep else "FAIL"
        record["filter_reason"] = (
            "PASS"
            if decision.keep
            else ";".join(decision.reasons)
        )

    def filter_file(
        self,
        input_path: str | Path,
        output_prefix: str | Path,
    ) -> FilterSummary:
        """Stream-filter an ORF message table.

        Args:
            input_path: Input message table.
            output_prefix: Output prefix.

        Returns:
            Filter summary.
        """
        header = ORFTable.read_header(input_path)
        ORFTable.validate_header(
            header,
            require_kozak=self.kozak_model is not None,
        )
        output_header = self.output_header(header)
        summary = FilterSummary()

        with _AtomicWriter(output_prefix, output_header) as writer:
            for record in ORFTable.iter_table(input_path, header):
                kozak_result = self.annotate_kozak(record)
                decision = self.decide(record, kozak_result)
                self._annotate_decision(record, decision)

                writer.write("all", record)
                writer.write(
                    "passed" if decision.keep else "removed",
                    record,
                )
                summary.update(record, decision)

                if summary.total % 100_000 == 0:
                    progress_print(
                        "records={total:,}, passed={passed:,}, "
                        "removed={removed:,}.".format(
                            total=summary.total,
                            passed=summary.passed,
                            removed=summary.removed,
                        )
                    )

            writer.write_summary(summary)

        return summary

    def add_kozak_pwm_fields(
        self,
        header: list[str],
        records: list[dict[str, str]],
    ) -> list[str]:
        """Compatibility method for in-memory callers.

        Args:
            header: Input header.
            records: Mutable records.

        Returns:
            Extended header.
        """
        for record in records:
            self.annotate_kozak(record)
        return self.output_header(header)

    # Backward-compatible table helpers.
    read_table = staticmethod(ORFTable.read_table)
    write_table = staticmethod(ORFTable.write_table)
