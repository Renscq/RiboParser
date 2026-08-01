#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Load, build, score, and export lightweight Kozak PWM models.
# Input: Built-in models, custom PWM files, or annotated ORF contexts.
# Output: Anchor-aligned Kozak scores between 0 and 1.

"""Lightweight Kozak context scoring for smORF filtering.

The public API intentionally contains one model class and one result class.
Sequences are aligned by the first nucleotide of the start codon. The three
start-codon positions are excluded from scoring because start-codon identity is
already filtered independently by the cluster structural rules.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

BASES = ("A", "C", "G", "T")
BASE_SET = frozenset(BASES)
DEFAULT_START_INDEX = 6
DEFAULT_VALID_RATIO = 0.5

# Compact A/C/G/T probability rows.
_PLANT = (
    (0.45, 0.35, 0.10, 0.10),
    (0.10, 0.65, 0.15, 0.10),
    (0.70, 0.10, 0.15, 0.05),
    (0.55, 0.35, 0.05, 0.05),
    (0.55, 0.25, 0.10, 0.10),
    (0.65, 0.15, 0.10, 0.10),
    (0.97, 0.01, 0.01, 0.01),
    (0.01, 0.01, 0.01, 0.97),
    (0.01, 0.01, 0.97, 0.01),
    (0.05, 0.05, 0.85, 0.05),
    (0.10, 0.65, 0.15, 0.10),
)
_VERTEBRATE = (
    (0.05, 0.10, 0.80, 0.05),
    (0.05, 0.80, 0.10, 0.05),
    (0.05, 0.80, 0.10, 0.05),
    (0.45, 0.05, 0.45, 0.05),
    (0.05, 0.80, 0.10, 0.05),
    (0.05, 0.80, 0.10, 0.05),
    (0.97, 0.01, 0.01, 0.01),
    (0.01, 0.01, 0.01, 0.97),
    (0.01, 0.01, 0.97, 0.01),
    (0.05, 0.05, 0.85, 0.05),
)
_DROSOPHILA = (
    (0.75, 0.05, 0.10, 0.10),
    (0.10, 0.05, 0.10, 0.75),
    (0.45, 0.45, 0.05, 0.05),
    (0.75, 0.05, 0.10, 0.10),
    (0.75, 0.05, 0.10, 0.10),
    (0.45, 0.45, 0.05, 0.05),
    (0.97, 0.01, 0.01, 0.01),
    (0.01, 0.01, 0.01, 0.97),
    (0.01, 0.01, 0.97, 0.01),
    (0.75, 0.05, 0.10, 0.10),
    (0.45, 0.45, 0.05, 0.05),
    (0.05, 0.75, 0.10, 0.10),
)
_YEAST = (
    (0.80, 0.05, 0.05, 0.10),
    (0.80, 0.05, 0.05, 0.10),
    (0.80, 0.05, 0.05, 0.10),
    (0.80, 0.05, 0.05, 0.10),
    (0.80, 0.05, 0.05, 0.10),
    (0.80, 0.05, 0.05, 0.10),
    (0.97, 0.01, 0.01, 0.01),
    (0.01, 0.01, 0.01, 0.97),
    (0.01, 0.01, 0.97, 0.01),
    (0.05, 0.05, 0.05, 0.85),
    (0.05, 0.75, 0.10, 0.10),
    (0.05, 0.10, 0.10, 0.75),
)

BUILTIN_MODELS = {
    "plant": _PLANT,
    "vertebrate": _VERTEBRATE,
    "drosophila": _DROSOPHILA,
    "yeast": _YEAST,
}

# Backward-compatible exported names.
BUILTIN_KOZAK_PWMS = BUILTIN_MODELS
BUILTIN_KOZAK_CONSENSUS = {
    "plant": "ACAACAATGGC",
    "vertebrate": "GCCRCCATGG",
    "drosophila": "ATMAAMATGAMC",
    "yeast": "AAAAAAATGTCT",
}


@dataclass(frozen=True, slots=True)
class KozakResult:
    """Store one Kozak scoring result.

    Attributes:
        score: Normalized score, or ``None`` for insufficient context.
        level: ``strong``, ``moderate``, ``weak``, or ``NA``.
        valid_ratio: Fraction of informative positions containing A/C/G/T.
    """

    score: float | None
    level: str
    valid_ratio: float


class KozakModel:
    """Represent a start-anchor-aligned Kozak PWM.

    Args:
        rows: PWM rows in A/C/G/T order.
        name: Human-readable model name.
        source: ``builtin``, ``custom``, or ``annotated``.
        start_index: First start-codon nucleotide in the model.
        training_count: Number of contexts used for an annotated model.

    Raises:
        ValueError: If the PWM or anchor is invalid.
    """

    def __init__(
        self,
        rows: Sequence[Sequence[float]],
        name: str,
        source: str,
        start_index: int = DEFAULT_START_INDEX,
        training_count: int = 0,
    ) -> None:
        self.rows = self._normalize_rows(rows)
        self.name = str(name)
        self.source = str(source)
        self.start_index = int(start_index)
        self.training_count = int(training_count)

        if not 0 <= self.start_index <= len(self.rows) - 3:
            raise ValueError("Kozak start index must identify a complete start codon.")

    @staticmethod
    def _normalize_rows(
        rows: Sequence[Sequence[float]],
    ) -> tuple[tuple[float, float, float, float], ...]:
        """Validate and normalize PWM rows.

        Args:
            rows: A/C/G/T weights.

        Returns:
            Normalized immutable rows.

        Raises:
            ValueError: If rows are empty or contain invalid values.
        """
        if not rows:
            raise ValueError("Kozak PWM is empty.")

        normalized = []
        for row_number, row in enumerate(rows, start=1):
            if len(row) != 4:
                raise ValueError(f"Kozak PWM row {row_number} must contain four values.")
            values = tuple(float(value) for value in row)
            if any(not math.isfinite(value) or value < 0 for value in values):
                raise ValueError(f"Invalid Kozak PWM value in row {row_number}.")

            total = sum(values)
            if total <= 0:
                raise ValueError(f"Kozak PWM row {row_number} has zero total weight.")
            normalized.append(tuple(value / total for value in values))
        return tuple(normalized)

    @classmethod
    def from_builtin(cls, name: str) -> KozakModel:
        """Load one built-in Kozak model.

        Args:
            name: Built-in model name.

        Returns:
            Kozak model.

        Raises:
            ValueError: If the model name is unknown.
        """
        key = str(name).lower()
        if key not in BUILTIN_MODELS:
            raise ValueError(
                f"Unknown Kozak model: {name}. Available: {', '.join(sorted(BUILTIN_MODELS))}"
            )
        return cls(
            rows=BUILTIN_MODELS[key],
            name=key,
            source="builtin",
        )

    @classmethod
    def from_pwm_file(
        cls,
        path: str | Path,
    ) -> KozakModel:
        """Load a tab-delimited custom PWM.

        The file must contain A, C, G, and T columns. A ``position`` column is
        optional; when present, position 0 identifies the start anchor.

        Args:
            path: PWM file.

        Returns:
            Custom Kozak model.

        Raises:
            ValueError: If the PWM table is malformed.
        """
        pwm_path = Path(path)
        with pwm_path.open("r", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n\r").split("\t")
            if not header or header == [""]:
                raise ValueError(f"Empty Kozak PWM file: {pwm_path}")

            names = [field.strip().upper() for field in header]
            indices = {name: index for index, name in enumerate(names)}
            missing = [base for base in BASES if base not in indices]
            if missing:
                raise ValueError("Kozak PWM is missing column(s): " + ", ".join(missing))

            position_index = indices.get("POSITION")
            positions = []
            rows = []

            for line_number, raw_line in enumerate(handle, start=2):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) != len(header):
                    raise ValueError(f"PWM field-count mismatch at line {line_number}.")
                try:
                    rows.append(tuple(float(fields[indices[base]]) for base in BASES))
                    if position_index is not None:
                        positions.append(int(fields[position_index]))
                except ValueError as error:
                    raise ValueError(f"Invalid PWM value at line {line_number}.") from error

        if not rows:
            raise ValueError(f"Empty Kozak PWM file: {pwm_path}")

        start_index = DEFAULT_START_INDEX
        if positions:
            if 0 not in positions:
                raise ValueError("PWM position column must contain position 0.")
            start_index = positions.index(0)

        return cls(
            rows=rows,
            name=pwm_path.stem,
            source="custom",
            start_index=start_index,
        )

    @classmethod
    def from_annotated_records(
        cls,
        records: Iterable[Mapping[str, Any]],
        minimum_records: int = 100,
    ) -> KozakModel:
        """Build a species-specific model from annotated ORFs.

        Only complete, primary, sense-strand ``annotated_ORF`` or
        ``annotated_mORF`` records are used. Duplicate genes and duplicate
        contexts are removed.

        Args:
            records: ORF message records.
            minimum_records: Minimum unique training contexts.

        Returns:
            Annotated Kozak model.

        Raises:
            ValueError: If too few usable records remain.
        """
        selected: list[tuple[str, int]] = []
        used_genes: set[str] = set()
        used_contexts: set[tuple[str, int]] = set()

        for record in records:
            if str(record.get("category", "")) not in {
                "annotated_ORF",
                "annotated_mORF",
            }:
                continue
            if str(record.get("source_strand", "sense")) != "sense":
                continue
            if str(record.get("priority", "primary")) != "primary":
                continue
            if str(record.get("completeness", "complete")) != "complete":
                continue

            sequence = cls._clean_sequence(record.get("kozak_seq", ""))
            if not sequence:
                continue

            try:
                start_index = cls.resolve_start_index(
                    sequence=sequence,
                    provided=record.get("kozak_start_index"),
                    start_codon=str(record.get("start_codon", "")),
                )
            except ValueError:
                continue

            gene_id = str(record.get("gene_id", "")).strip()
            context_key = (sequence, start_index)
            if context_key in used_contexts:
                continue
            if gene_id and gene_id != "." and gene_id in used_genes:
                continue

            used_contexts.add(context_key)
            if gene_id and gene_id != ".":
                used_genes.add(gene_id)
            selected.append(context_key)

        if len(selected) < int(minimum_records):
            raise ValueError(
                f"Too few annotated Kozak contexts: {len(selected)} < {int(minimum_records)}"
            )

        upstream = Counter(index for _, index in selected).most_common(1)[0][0]
        downstream = Counter(len(sequence) - index - 3 for sequence, index in selected).most_common(
            1
        )[0][0]
        length = upstream + 3 + downstream

        counts = [[1.0, 1.0, 1.0, 1.0] for _ in range(length)]
        for sequence, start_index in selected:
            aligned = cls._align(
                sequence=sequence,
                sequence_start=start_index,
                target_length=length,
                target_start=upstream,
            )
            for position, base in enumerate(aligned):
                if base in BASE_SET:
                    counts[position][BASES.index(base)] += 1.0

        return cls(
            rows=counts,
            name="annotated",
            source="annotated",
            start_index=upstream,
            training_count=len(selected),
        )

    @staticmethod
    def _clean_sequence(value: Any) -> str:
        """Return uppercase DNA context."""
        return str(value).strip().upper().replace("U", "T")

    @staticmethod
    def resolve_start_index(
        sequence: str,
        provided: Any = None,
        start_codon: str = "",
    ) -> int:
        """Resolve the start-codon index.

        Args:
            sequence: Kozak context.
            provided: Optional explicit scanner index.
            start_codon: Expected start codon.

        Returns:
            Resolved zero-based index.

        Raises:
            ValueError: If no valid anchor can be found.
        """
        sequence = KozakModel._clean_sequence(sequence)
        expected = KozakModel._clean_sequence(start_codon)

        try:
            explicit = int(provided)
        except (TypeError, ValueError):
            explicit = None

        if explicit is not None and 0 <= explicit <= len(sequence) - 3:
            if not expected or sequence[explicit : explicit + 3] == expected:
                return explicit

        if len(expected) == 3:
            candidates = [
                index
                for index in range(len(sequence) - 2)
                if sequence[index : index + 3] == expected
            ]
            if candidates:
                return min(
                    candidates,
                    key=lambda index: (
                        abs(index - DEFAULT_START_INDEX),
                        index,
                    ),
                )

        if len(sequence) - 3 >= DEFAULT_START_INDEX:
            return DEFAULT_START_INDEX
        raise ValueError("Cannot resolve Kozak start-codon index.")

    @staticmethod
    def _align(
        sequence: str,
        sequence_start: int,
        target_length: int,
        target_start: int,
    ) -> str:
        """Align a context by its start-codon anchor."""
        output = ["N"] * target_length
        offset = target_start - sequence_start
        for source_index, base in enumerate(sequence):
            target_index = source_index + offset
            if 0 <= target_index < target_length:
                output[target_index] = base
        return "".join(output)

    def score(
        self,
        sequence: str,
        start_index: Any = None,
        start_codon: str = "",
    ) -> KozakResult:
        """Score Kozak flanking context.

        Args:
            sequence: Kozak context.
            start_index: Optional scanner-provided anchor.
            start_codon: Expected start codon.

        Returns:
            Kozak score and context-validity metadata.
        """
        sequence = self._clean_sequence(sequence)
        try:
            resolved = self.resolve_start_index(
                sequence=sequence,
                provided=start_index,
                start_codon=start_codon,
            )
        except ValueError:
            return KozakResult(None, "NA", 0.0)

        aligned = self._align(
            sequence=sequence,
            sequence_start=resolved,
            target_length=len(self.rows),
            target_start=self.start_index,
        )
        excluded = {
            self.start_index,
            self.start_index + 1,
            self.start_index + 2,
        }

        values = []
        informative = 0
        for position, (base, row) in enumerate(zip(aligned, self.rows)):
            if position in excluded:
                continue

            row_min = min(row)
            row_max = max(row)
            if math.isclose(row_min, row_max, abs_tol=1e-15):
                continue

            informative += 1
            if base not in BASE_SET:
                continue

            probability = row[BASES.index(base)]
            value = (math.log(probability) - math.log(row_min)) / (
                math.log(row_max) - math.log(row_min)
            )
            values.append(max(0.0, min(1.0, value)))

        valid_ratio = len(values) / informative if informative else 0.0
        if not values or valid_ratio < DEFAULT_VALID_RATIO:
            return KozakResult(None, "NA", valid_ratio)

        score = sum(values) / len(values)
        if score >= 0.75:
            level = "strong"
        elif score >= 0.50:
            level = "moderate"
        else:
            level = "weak"
        return KozakResult(score, level, valid_ratio)

    def export(self, path: str | Path) -> None:
        """Export the PWM with start-relative coordinates.

        Args:
            path: Output file.
        """
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w", encoding="utf-8") as handle:
            handle.write("position\tA\tC\tG\tT\n")
            for index, row in enumerate(self.rows):
                relative = index - self.start_index
                handle.write(f"{relative}\t" + "\t".join(f"{value:.6f}" for value in row) + "\n")


# Backward-compatible class alias.
KozakPWM = KozakModel
