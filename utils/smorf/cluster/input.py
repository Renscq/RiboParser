#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Read scanner tables and compile column positions.
# Input: smorf_scanner message table.
# Output: Validated rows and binary column plan.

"""Read scanner tables and compile column positions."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from .config import (
    _BUFFER_BYTES,
    ANNOTATED_OVERLAP_COLUMNS,
    CLUSTER_COLUMNS,
    CLUSTER_REQUIRED_COLUMNS,
    FILTER_COLUMNS,
)

_REQUIRED_COLUMNS = frozenset(name.decode("ascii") for name in CLUSTER_REQUIRED_COLUMNS)
_KOZAK_COLUMNS = frozenset({"kozak_seq", "kozak_start_index"})


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
    def from_header(cls, header: Sequence[bytes]) -> _ColumnPlan:
        """Build and validate one column plan."""
        input_header = tuple(header)
        lookup = {name: index for index, name in enumerate(input_header)}
        missing = [name for name in CLUSTER_REQUIRED_COLUMNS if name not in lookup]
        if missing:
            raise ValueError(
                "ORF message table is missing clustering column(s): "
                + ", ".join(item.decode("utf-8", errors="replace") for item in missing)
            )

        filter_header = list(input_header)
        filter_lookup = dict(lookup)
        for name in (*FILTER_COLUMNS, *ANNOTATED_OVERLAP_COLUMNS):
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


class ORFTable:
    """Read and validate strict tab-delimited scanner message tables."""

    @staticmethod
    def read_header(path: str | Path) -> list[str]:
        """Read a non-empty header and reject duplicate column names."""
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
        duplicates = sorted(name for name, count in Counter(header).items() if count > 1)
        if duplicates:
            raise ValueError("Duplicate ORF column(s): " + ", ".join(duplicates))
        return header

    @staticmethod
    def validate_header(
        header: Sequence[str],
        require_kozak: bool = False,
    ) -> None:
        """Validate columns required by clustering and optional Kozak training."""
        required = set(_REQUIRED_COLUMNS)
        if require_kozak:
            required.update(_KOZAK_COLUMNS)
        missing = sorted(required.difference(header))
        if missing:
            raise ValueError("ORF message table is missing column(s): " + ", ".join(missing))

    @staticmethod
    def iter_table(
        path: str | Path,
        header: Sequence[str] | None = None,
    ) -> Iterator[dict[str, str]]:
        """Stream validated rows as mappings without loading the full table."""
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
                        f"{line_number}: expected {expected_count}, observed {len(fields)}."
                    )
                yield dict(zip(expected, fields))
