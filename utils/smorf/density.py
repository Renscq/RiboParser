#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Parse, normalize, sort, and query sparse WIG or bedGraph P-site density.
# Input: Plain or gzip-compressed WIG/bedGraph density files.
# Output: Chromosome-level sparse absolute-density objects.

"""Sparse P-site density readers with automatic bedGraph order recovery.

All density values are converted to absolute values. This supports negative
minus-strand bedGraph conventions without changing the biological strand
assigned in the density list.

bedGraph files are checked for chromosome grouping and coordinate order before
analysis. Files with separated chromosome blocks or decreasing coordinates are
automatically sorted to a temporary file with GNU ``sort``. WIG files retain
their native step-block order requirement.
"""

from __future__ import annotations

import gzip
import math
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import numpy as np

SUPPORTED_DENSITY_FORMATS = frozenset({"auto", "wig", "bedgraph"})


@dataclass(frozen=True, slots=True)
class ChromDensity:
    """Store one chromosome as sorted non-overlapping sparse intervals.

    Attributes:
        chrom: Chromosome name.
        starts: Zero-based interval starts.
        ends: Zero-based interval ends.
        values: Non-negative density values.
        max_end: Maximum observed interval end.
    """

    chrom: str
    starts: np.ndarray
    ends: np.ndarray
    values: np.ndarray
    max_end: int

    def query(self, start: int, end: int) -> np.ndarray:
        """Return density for a genomic half-open interval.

        Args:
            start: Zero-based inclusive query start.
            end: Zero-based exclusive query end.

        Returns:
            Float32 density vector of length ``end - start``.

        Raises:
            ValueError: If the query interval is invalid.
        """
        start = int(start)
        end = int(end)
        if start < 0 or end < start:
            raise ValueError(f"Invalid density query interval: [{start}, {end}).")
        if end == start:
            return np.zeros(0, dtype=np.float32)

        output = np.zeros(end - start, dtype=np.float32)
        if self.starts.size == 0:
            return output

        left = int(np.searchsorted(self.ends, start, side="right"))
        right = int(np.searchsorted(self.starts, end, side="left"))

        for index in range(left, right):
            overlap_start = max(start, int(self.starts[index]))
            overlap_end = min(end, int(self.ends[index]))
            if overlap_end <= overlap_start:
                continue
            output[overlap_start - start : overlap_end - start] = self.values[index]

        return output


@dataclass(slots=True)
class PreparedDensity:
    """Store a density path prepared for one analysis pass.

    Attributes:
        path: Original or temporary sorted path.
        temporary: Whether ``path`` should be deleted after use.
        was_sorted: Whether automatic sorting was performed.
    """

    path: str
    temporary: bool
    was_sorted: bool

    def cleanup(self) -> None:
        """Delete an automatically generated temporary file."""
        if self.temporary:
            Path(self.path).unlink(missing_ok=True)


def smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open a plain or gzip-compressed file.

    Args:
        path: Input path.
        mode: File-open mode.

    Returns:
        File handle.
    """
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, mode)
    return file_path.open(mode)


def infer_density_format(
    path: str | Path,
    user_format: str | None = None,
) -> str:
    """Infer WIG or bedGraph format.

    Args:
        path: Density file path.
        user_format: Optional explicit format.

    Returns:
        ``wig`` or ``bedgraph``.

    Raises:
        ValueError: If the explicit or inferred format is unsupported.
    """
    if user_format and user_format != "auto":
        format_name = str(user_format).lower()
        if format_name not in SUPPORTED_DENSITY_FORMATS:
            raise ValueError(f"Unsupported density format: {format_name}")
        return format_name

    lower = str(path).lower()
    if lower.endswith((".bedgraph", ".bdg", ".bedgraph.gz", ".bdg.gz")):
        return "bedgraph"
    if lower.endswith((".wig", ".wiggle", ".wig.gz", ".wiggle.gz")):
        return "wig"
    raise ValueError(
        "Cannot infer density format from filename. Use --density-format wig or bedgraph."
    )


def _iter_bedgraph_data_lines(
    path: str | Path,
) -> Iterator[tuple[int, str]]:
    """Yield bedGraph data lines with line numbers."""
    with smart_open(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            text = raw_line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            yield line_number, text


def bedgraph_needs_sort(path: str | Path) -> bool:
    """Check chromosome grouping and coordinate order.

    Args:
        path: bedGraph path.

    Returns:
        ``True`` when automatic sorting is required.

    Raises:
        ValueError: If a bedGraph line is malformed.
    """
    current_chrom: str | None = None
    closed_chromosomes: set[str] = set()
    last_start = -1
    last_end = -1

    for line_number, text in _iter_bedgraph_data_lines(path):
        fields = text.split()
        if len(fields) < 4:
            raise ValueError(f"bedGraph line {line_number} has fewer than four fields.")
        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as error:
            raise ValueError(f"Invalid bedGraph coordinates at line {line_number}.") from error

        if current_chrom is None:
            current_chrom = chrom
            last_start = start
            last_end = end
            continue

        if chrom != current_chrom:
            closed_chromosomes.add(current_chrom)
            if chrom in closed_chromosomes:
                return True
            current_chrom = chrom
            last_start = start
            last_end = end
            continue

        if start < last_start or (start == last_start and end < last_end):
            return True
        last_start = start
        last_end = end

    return False


def _stream_bedgraph_to_sort(
    path: str | Path,
    process: subprocess.Popen[bytes],
) -> None:
    """Stream data-only bedGraph lines to GNU sort."""
    if process.stdin is None:
        raise RuntimeError("GNU sort stdin was not created.")

    try:
        for _, text in _iter_bedgraph_data_lines(path):
            process.stdin.write(text.encode("utf-8") + b"\n")
        process.stdin.close()
    except BaseException:
        process.kill()
        raise


def sort_bedgraph(
    path: str | Path,
    work_directory: str | Path,
) -> PreparedDensity:
    """Sort a bedGraph by chromosome, start, and end.

    Args:
        path: Original bedGraph path.
        work_directory: Directory for the temporary sorted file.

    Returns:
        Prepared temporary density path.

    Raises:
        RuntimeError: If GNU sort is unavailable or fails.
    """
    sort_program = shutil.which("sort")
    if sort_program is None:
        raise RuntimeError(
            "The bedGraph is not grouped by chromosome and GNU sort "
            "is unavailable. Sort it with: "
            "LC_ALL=C sort -k1,1 -k2,2n -k3,3n input.bedgraph"
        )

    work_path = Path(work_directory)
    work_path.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_path = tempfile.mkstemp(
        prefix="smorf_evidence.",
        suffix=".sorted.bedGraph",
        dir=work_path,
        text=False,
    )
    os.close(descriptor)

    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    command = [
        sort_program,
        "-T",
        str(work_path),
        "-k1,1",
        "-k2,2n",
        "-k3,3n",
    ]

    try:
        with Path(temporary_path).open("wb") as output_handle:
            process = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=output_handle,
                stderr=subprocess.PIPE,
                env=environment,
            )
            _stream_bedgraph_to_sort(path, process)
            if process.stderr is not None:
                stderr = process.stderr.read().decode(
                    "utf-8",
                    errors="replace",
                )
                process.stderr.close()
            else:
                stderr = ""
            return_code = process.wait()

        if return_code != 0:
            raise RuntimeError("Automatic bedGraph sorting failed: " + stderr.strip())

        return PreparedDensity(
            path=temporary_path,
            temporary=True,
            was_sorted=True,
        )
    except BaseException:
        Path(temporary_path).unlink(missing_ok=True)
        raise


def prepare_density(
    path: str | Path,
    file_format: str,
    work_directory: str | Path,
) -> PreparedDensity:
    """Prepare a density file for one ordered analysis pass.

    Args:
        path: Density path.
        file_format: Explicit or automatic format.
        work_directory: Temporary-sort directory.

    Returns:
        Original or temporary prepared path.
    """
    format_name = infer_density_format(path, file_format)
    if format_name != "bedgraph":
        return PreparedDensity(
            path=str(path),
            temporary=False,
            was_sorted=False,
        )

    if not bedgraph_needs_sort(path):
        return PreparedDensity(
            path=str(path),
            temporary=False,
            was_sorted=False,
        )
    return sort_bedgraph(path, work_directory)


def _validate_interval(
    chrom: str,
    start: int,
    end: int,
    value: float,
    line_number: int,
) -> None:
    """Validate one density interval."""
    if not chrom:
        raise ValueError(f"Empty chromosome at density line {line_number}.")
    if start < 0 or end <= start:
        raise ValueError(f"Invalid density interval at line {line_number}: {chrom}:{start}-{end}")
    if not math.isfinite(value):
        raise ValueError(f"Invalid density value at line {line_number}: {value}")


def _build_chrom_density(
    chrom: str,
    records: list[tuple[int, int, float]],
) -> ChromDensity:
    """Build a sorted and merged chromosome density object."""
    if not records:
        return ChromDensity(
            chrom=chrom,
            starts=np.zeros(0, dtype=np.int64),
            ends=np.zeros(0, dtype=np.int64),
            values=np.zeros(0, dtype=np.float32),
            max_end=0,
        )

    records.sort(key=lambda item: (item[0], item[1]))
    merged: list[list[float | int]] = []

    for start, end, value in records:
        value = abs(float(value))
        if merged and start < int(merged[-1][1]):
            raise ValueError(
                f"Overlapping density intervals on {chrom}: "
                f"[{int(merged[-1][0])}, {int(merged[-1][1])}) and "
                f"[{start}, {end})."
            )
        if (
            merged
            and start == int(merged[-1][1])
            and math.isclose(
                float(merged[-1][2]),
                value,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ):
            merged[-1][1] = end
        else:
            merged.append([start, end, value])

    starts = np.asarray(
        [int(record[0]) for record in merged],
        dtype=np.int64,
    )
    ends = np.asarray(
        [int(record[1]) for record in merged],
        dtype=np.int64,
    )
    values = np.asarray(
        [float(record[2]) for record in merged],
        dtype=np.float32,
    )

    return ChromDensity(
        chrom=chrom,
        starts=starts,
        ends=ends,
        values=values,
        max_end=int(ends[-1]),
    )


def _parse_wig_attrs(line: str) -> dict[str, str]:
    """Parse fixedStep or variableStep WIG attributes."""
    attributes = {}
    for item in line.split()[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def _yield_chromosome(
    chrom: str | None,
    records: list[tuple[int, int, float]],
    seen: set[str],
) -> ChromDensity | None:
    """Finalize one chromosome and validate WIG/block ordering."""
    if chrom is None:
        return None
    if chrom in seen:
        raise ValueError(f"Density chromosome appears in multiple separated blocks: {chrom}.")
    seen.add(chrom)
    return _build_chrom_density(chrom, records)


def _iter_bedgraph(path: str | Path) -> Iterator[ChromDensity]:
    """Yield chromosome density from an ordered bedGraph."""
    current_chrom: str | None = None
    records: list[tuple[int, int, float]] = []
    seen: set[str] = set()

    for line_number, text in _iter_bedgraph_data_lines(path):
        fields = text.split()
        if len(fields) < 4:
            raise ValueError(f"bedGraph line {line_number} has fewer than four fields.")

        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
            value = abs(float(fields[3]))
        except ValueError as error:
            raise ValueError(f"Invalid bedGraph numeric value at line {line_number}.") from error

        _validate_interval(
            chrom=chrom,
            start=start,
            end=end,
            value=value,
            line_number=line_number,
        )

        if current_chrom is None:
            current_chrom = chrom
        elif chrom != current_chrom:
            result = _yield_chromosome(
                current_chrom,
                records,
                seen,
            )
            if result is not None:
                yield result
            current_chrom = chrom
            records = []

        records.append((start, end, value))

    result = _yield_chromosome(current_chrom, records, seen)
    if result is not None:
        yield result


def _iter_wig(path: str | Path) -> Iterator[ChromDensity]:
    """Yield chromosome density from a WIG file."""
    current_chrom: str | None = None
    records: list[tuple[int, int, float]] = []
    seen: set[str] = set()

    mode: str | None = None
    position: int | None = None
    step = 1
    span = 1

    with smart_open(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            text = raw_line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue

            if text.startswith("fixedStep"):
                attributes = _parse_wig_attrs(text)
                chrom = attributes.get("chrom", "")
                if not chrom:
                    raise ValueError(f"fixedStep header lacks chrom at line {line_number}.")
                if current_chrom is None:
                    current_chrom = chrom
                elif chrom != current_chrom:
                    result = _yield_chromosome(
                        current_chrom,
                        records,
                        seen,
                    )
                    if result is not None:
                        yield result
                    current_chrom = chrom
                    records = []

                try:
                    position = int(attributes.get("start", "1")) - 1
                    step = int(attributes.get("step", "1"))
                    span = int(attributes.get("span", "1"))
                except ValueError as error:
                    raise ValueError(f"Invalid fixedStep header at line {line_number}.") from error
                if position < 0 or step < 1 or span < 1:
                    raise ValueError(f"Invalid fixedStep coordinates at line {line_number}.")
                mode = "fixed"
                continue

            if text.startswith("variableStep"):
                attributes = _parse_wig_attrs(text)
                chrom = attributes.get("chrom", "")
                if not chrom:
                    raise ValueError(f"variableStep header lacks chrom at line {line_number}.")
                if current_chrom is None:
                    current_chrom = chrom
                elif chrom != current_chrom:
                    result = _yield_chromosome(
                        current_chrom,
                        records,
                        seen,
                    )
                    if result is not None:
                        yield result
                    current_chrom = chrom
                    records = []

                try:
                    span = int(attributes.get("span", "1"))
                except ValueError as error:
                    raise ValueError(
                        f"Invalid variableStep header at line {line_number}."
                    ) from error
                if span < 1:
                    raise ValueError(f"Invalid variableStep span at line {line_number}.")
                mode = "variable"
                position = None
                continue

            if current_chrom is None or mode is None:
                raise ValueError(f"WIG value occurs before a step header at line {line_number}.")

            if mode == "fixed":
                if position is None:
                    raise ValueError(f"Missing fixedStep position at line {line_number}.")
                try:
                    value = abs(float(text.split()[0]))
                except ValueError as error:
                    raise ValueError(f"Invalid WIG value at line {line_number}.") from error
                start = position
                end = position + span
                position += step
            else:
                fields = text.split()
                if len(fields) < 2:
                    raise ValueError(f"Invalid variableStep line {line_number}.")
                try:
                    start = int(fields[0]) - 1
                    value = abs(float(fields[1]))
                except ValueError as error:
                    raise ValueError(
                        f"Invalid variableStep value at line {line_number}."
                    ) from error
                end = start + span

            _validate_interval(
                chrom=current_chrom,
                start=start,
                end=end,
                value=value,
                line_number=line_number,
            )
            records.append((start, end, value))

    result = _yield_chromosome(current_chrom, records, seen)
    if result is not None:
        yield result


def iter_density_chromosomes(
    path: str | Path,
    file_format: str = "auto",
) -> Iterator[ChromDensity]:
    """Yield sparse chromosome densities in one ordered file pass.

    Args:
        path: Prepared WIG or bedGraph path.
        file_format: ``auto``, ``wig``, or ``bedgraph``.

    Yields:
        Chromosome density objects.
    """
    format_name = infer_density_format(path, file_format)
    if format_name == "bedgraph":
        yield from _iter_bedgraph(path)
    else:
        yield from _iter_wig(path)
