#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Read and validate smORF tables, density lists, and genePred annotation.
# Input: ORF message tables, density-list tables, and optional genePred files.
# Output: Validated pandas tables and annotation mappings.

"""Input and output helpers for smORF Ribo-seq evidence analysis."""

from __future__ import annotations

import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import pandas as pd

from .smorf_riboseq_constants import (
    DensityTrack,
    SUPPORTED_DENSITY_FORMATS,
    VALID_STRANDS,
)


REQUIRED_ORF_COLUMNS = (
    "orf_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "genomic_start",
    "genomic_end",
    "nt_length",
)


@dataclass(frozen=True, slots=True)
class GenePredRecord:
    """Store one genePred transcript model.

    Attributes:
        name: genePred record name.
        chrom: Chromosome name.
        strand: ``+`` or ``-``.
        tx_start: Genomic transcript start.
        tx_end: Genomic transcript end.
        exon_starts: Exon starts in ascending genomic order.
        exon_ends: Exon ends in ascending genomic order.
    """

    name: str
    chrom: str
    strand: str
    tx_start: int
    tx_end: int
    exon_starts: tuple[int, ...]
    exon_ends: tuple[int, ...]


def eprint(message: str) -> None:
    """Print one progress message to stderr.

    Args:
        message: Message text.
    """
    print(message, file=sys.stderr, flush=True)


def smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open a plain or gzip-compressed text file.

    Args:
        path: Input path.
        mode: Open mode.

    Returns:
        File handle.
    """
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, mode)
    return file_path.open(mode)


def parse_comma_ints(value: object) -> list[int]:
    """Parse a comma-separated integer field.

    Args:
        value: Table value.

    Returns:
        Parsed integers.

    Raises:
        ValueError: If any token is invalid.
    """
    if pd.isna(value):
        return []

    text = str(value).strip().rstrip(",")
    if not text or text == ".":
        return []

    try:
        return [int(token) for token in text.split(",") if token]
    except ValueError as error:
        raise ValueError(
            f"Invalid comma-separated integer field: {value}"
        ) from error


def normalize_density_strand(value: object) -> str:
    """Normalize a density strand label.

    Args:
        value: Raw strand value.

    Returns:
        ``+``, ``-``, or ``.``.

    Raises:
        ValueError: If the value is unknown.
    """
    text = "" if pd.isna(value) else str(value).strip().lower()
    aliases = {
        "+": "+",
        "plus": "+",
        "forward": "+",
        "fwd": "+",
        "sense": "+",
        "-": "-",
        "minus": "-",
        "reverse": "-",
        "rev": "-",
        "antisense": "-",
        ".": ".",
        "both": ".",
        "all": ".",
        "unstranded": ".",
        "none": ".",
        "": ".",
    }
    if text not in aliases:
        raise ValueError(
            f"Invalid density strand: {value}. "
            "Use +, -, or unstranded."
        )
    return aliases[text]


def infer_density_strand(path: str | Path) -> str:
    """Infer density strand from a file name.

    Recognized plus labels include ``plus``, ``forward``, and ``fwd``.
    Recognized minus labels include ``minus``, ``reverse``, and ``rev``.
    Files without an unambiguous label are treated as unstranded.

    Args:
        path: Density file path.

    Returns:
        ``+``, ``-``, or ``.``.
    """
    name = Path(path).name.lower()
    tokens = [
        token
        for token in re.split(r"[^a-z0-9]+", name)
        if token
    ]
    plus = any(
        token in {"plus", "forward", "fwd", "pos", "positive"}
        for token in tokens
    )
    minus = any(
        token in {"minus", "reverse", "rev", "neg", "negative"}
        for token in tokens
    )

    if plus and not minus:
        return "+"
    if minus and not plus:
        return "-"
    return "."


def _normalize_format(value: object) -> str:
    """Normalize a density format value."""
    text = "auto" if pd.isna(value) else str(value).strip().lower()
    if not text:
        text = "auto"
    if text not in SUPPORTED_DENSITY_FORMATS:
        raise ValueError(
            f"Unsupported density format: {value}"
        )
    return text


def read_density_list(args: object) -> list[DensityTrack]:
    """Read and validate density-track definitions.

    ``--density`` may be repeated. Strand is inferred from each file name when
    direct input is used. The density-list ``strand`` column remains supported
    and takes precedence over file-name inference.

    A sample may use either one unstranded track or a plus/minus pair.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Validated density tracks in input order.

    Raises:
        ValueError: If definitions are missing, duplicated, or conflicting.
        FileNotFoundError: If a density file does not exist.
    """
    direct_density = list(getattr(args, "density", None) or [])
    legacy_plus = getattr(args, "density_plus", None)
    legacy_minus = getattr(args, "density_minus", None)
    if legacy_plus:
        direct_density.append(legacy_plus)
    if legacy_minus:
        direct_density.append(legacy_minus)

    has_direct = bool(direct_density)
    density_list = getattr(args, "density_list", None)
    if density_list and has_direct:
        raise ValueError(
            "--density-list cannot be combined with direct --density files."
        )

    tracks: list[DensityTrack] = []
    if density_list:
        table = pd.read_csv(
            density_list,
            sep="\t",
            comment="#",
            dtype=str,
            keep_default_na=False,
        )
        table.columns = [
            str(column).strip().lower()
            for column in table.columns
        ]
        missing = {"sample", "path"} - set(table.columns)
        if missing:
            raise ValueError(
                "Density list is missing column(s): "
                + ", ".join(sorted(missing))
            )
        if "format" not in table.columns:
            table["format"] = "auto"

        has_strand_column = "strand" in table.columns
        for row_number, row in enumerate(
            table.itertuples(index=False),
            start=2,
        ):
            row_data = row._asdict()
            sample = str(row_data["sample"]).strip()
            path = str(row_data["path"]).strip()
            if not sample or not path:
                raise ValueError(
                    f"Empty sample or path in density-list line "
                    f"{row_number}."
                )

            raw_strand = (
                row_data.get("strand", "")
                if has_strand_column
                else ""
            )
            strand = (
                normalize_density_strand(raw_strand)
                if str(raw_strand).strip()
                else infer_density_strand(path)
            )
            tracks.append(
                DensityTrack(
                    sample=sample,
                    strand=strand,
                    path=path,
                    file_format=_normalize_format(
                        row_data.get("format", "auto")
                    ),
                )
            )
    else:
        sample = str(getattr(args, "sample", "sample1")).strip()
        file_format = _normalize_format(
            getattr(args, "density_format", "auto")
        )
        for path in direct_density:
            if path == legacy_plus:
                strand = "+"
            elif path == legacy_minus:
                strand = "-"
            else:
                strand = infer_density_strand(path)
            tracks.append(
                DensityTrack(
                    sample=sample,
                    strand=strand,
                    path=str(path),
                    file_format=file_format,
                )
            )

    if not tracks:
        raise ValueError(
            "Provide --density-list or at least one --density file."
        )

    seen_pairs: set[tuple[str, str]] = set()
    sample_strands: dict[str, set[str]] = {}
    seen_paths: set[str] = set()

    for track in tracks:
        if track.strand not in VALID_STRANDS:
            raise ValueError(
                f"Invalid density-track strand: {track.strand}"
            )

        path = Path(track.path)
        if not path.is_file():
            raise FileNotFoundError(path)

        canonical_path = str(path.resolve())
        if canonical_path in seen_paths:
            raise ValueError(
                f"Density file is listed more than once: {track.path}"
            )
        seen_paths.add(canonical_path)

        pair = (track.sample, track.strand)
        if pair in seen_pairs:
            raise ValueError(
                "Duplicate density sample/strand definition: "
                f"{track.sample}/{track.strand}"
            )
        seen_pairs.add(pair)
        sample_strands.setdefault(track.sample, set()).add(track.strand)

    for sample, strands in sample_strands.items():
        if "." in strands and len(strands) > 1:
            raise ValueError(
                f"Sample {sample} mixes unstranded and strand-specific "
                "density tracks."
            )
        if strands not in ({'.'}, {'+', '-'}):
            raise ValueError(
                f"Sample {sample} requires either one unstranded track "
                "or a complete plus/minus pair. Add a strand column to "
                "the density list when file names do not contain plus/minus."
            )

    return tracks


def _coerce_integer_column(
    table: pd.DataFrame,
    column: str,
) -> None:
    """Strictly convert one table column to integer."""
    numeric = pd.to_numeric(table[column], errors="coerce")
    invalid = numeric.isna()
    if invalid.any():
        examples = table.loc[invalid, column].head(5).tolist()
        raise ValueError(
            f"Invalid integer values in {column}: {examples}"
        )
    table[column] = numeric.astype("int64")


def read_orf_table(
    path: str | Path,
    coord_mode: str = "0based-half-open",
) -> pd.DataFrame:
    """Read and strictly validate a filtered smORF table.

    Rows with ``filter_status`` are restricted to ``PASS``. Invalid numeric
    fields, duplicated ORF IDs, unsupported strands, and inconsistent genomic
    bounds raise errors instead of being silently dropped.

    Args:
        path: Filtered ORF message table.
        coord_mode: Input coordinate mode.

    Returns:
        Validated ORF table.

    Raises:
        ValueError: If the table is malformed.
    """
    with smart_open(path, "rt") as handle:
        header_line = handle.readline().rstrip("\n\r")
    raw_header = header_line.split("\t") if header_line else []
    duplicated_columns = [
        column
        for column, count in pd.Series(raw_header).value_counts().items()
        if count > 1
    ]
    if duplicated_columns:
        raise ValueError(
            "Duplicate ORF table column(s): "
            + ", ".join(sorted(duplicated_columns))
        )

    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    missing = [
        column
        for column in REQUIRED_ORF_COLUMNS
        if column not in table.columns
    ]
    if missing:
        raise ValueError(
            "ORF table is missing required column(s): "
            + ", ".join(missing)
        )

    if "filter_status" in table.columns:
        table = table.loc[
            table["filter_status"].astype(str).str.upper().eq("PASS")
        ].copy()

    if table.empty:
        raise ValueError("No PASS ORF remains for evidence analysis.")

    for column in ("genomic_start", "genomic_end", "nt_length"):
        _coerce_integer_column(table, column)

    if "aa_length" in table.columns:
        _coerce_integer_column(table, "aa_length")

    if coord_mode == "1based-closed":
        table["genomic_start"] -= 1
    elif coord_mode != "0based-half-open":
        raise ValueError(f"Unsupported coordinate mode: {coord_mode}")

    invalid_interval = (
        (table["genomic_start"] < 0)
        | (table["genomic_end"] <= table["genomic_start"])
        | (table["nt_length"] <= 0)
    )
    if invalid_interval.any():
        examples = table.loc[
            invalid_interval,
            ["orf_id", "genomic_start", "genomic_end", "nt_length"],
        ].head(5)
        raise ValueError(
            "Invalid ORF interval or length:\n"
            + examples.to_string(index=False)
        )

    invalid_strand = ~table["strand"].isin(["+", "-"])
    if invalid_strand.any():
        examples = table.loc[
            invalid_strand,
            ["orf_id", "strand"],
        ].head(5)
        raise ValueError(
            "Invalid ORF genomic strand:\n"
            + examples.to_string(index=False)
        )

    duplicate_ids = table["orf_id"].duplicated(keep=False)
    if duplicate_ids.any():
        examples = table.loc[
            duplicate_ids,
            "orf_id",
        ].head(10).tolist()
        raise ValueError(
            f"Duplicate ORF IDs are not supported: {examples}"
        )

    table["orf_id"] = table["orf_id"].astype(str)
    table["gene_id"] = table["gene_id"].astype(str)
    table["transcript_id"] = table["transcript_id"].astype(str)
    table["chrom"] = table["chrom"].astype(str)
    table["strand"] = table["strand"].astype(str)

    return table.reset_index(drop=True)


def peek_genepred_names(
    path: str | Path,
    maximum_records: int = 200,
) -> list[str]:
    """Read a small genePred name sample without parsing records.

    Args:
        path: genePred path.
        maximum_records: Maximum non-comment names returned.

    Returns:
        Sampled record names in file order.
    """
    names: list[str] = []
    with smart_open(path, "rt") as handle:
        for raw_line in handle:
            text = raw_line.strip()
            if not text or text.startswith("#"):
                continue
            names.append(text.split("\t", 1)[0].strip())
            if len(names) >= int(maximum_records):
                break
    return names


def read_genepred(
    path: str | Path | None,
    coord_mode: str = "0based-half-open",
    required_names: set[str] | None = None,
) -> dict[str, GenePredRecord]:
    """Read selected genePred records keyed by name.

    When ``required_names`` is provided, non-matching rows are skipped after
    reading only the first field. This prevents scanner-generated genePred files
    containing millions of ORFs from being retained entirely in memory.

    For transcript-aware release analysis, the genePred names must match the
    ORF table ``transcript_id`` values. The scanner ORF genePred output normally
    does not satisfy this requirement.

    Args:
        path: Optional genePred path.
        coord_mode: Input coordinate mode.
        required_names: Optional record-name allowlist.

    Returns:
        Selected genePred records keyed by name.

    Raises:
        ValueError: If a selected genePred record is malformed or duplicated.
    """
    if path is None:
        return {}

    wanted = required_names
    records: dict[str, GenePredRecord] = {}

    with smart_open(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            text = raw_line.strip()
            if not text or text.startswith("#"):
                continue

            name, separator, remainder = text.partition("\t")
            name = name.strip()
            if not separator:
                raise ValueError(
                    f"genePred line {line_number} has fewer than 10 columns."
                )
            if wanted is not None and name not in wanted:
                continue

            fields = [name, *remainder.split("\t")]
            if len(fields) < 10:
                raise ValueError(
                    f"genePred line {line_number} has fewer than "
                    "10 columns."
                )

            chrom = fields[1].strip()
            strand = fields[2].strip()
            if name in records:
                raise ValueError(
                    f"Duplicate selected genePred record name: {name}"
                )
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid genePred strand at line {line_number}."
                )

            try:
                tx_start = int(fields[3])
                tx_end = int(fields[4])
                exon_count = int(fields[7])
                starts = parse_comma_ints(fields[8])
                ends = parse_comma_ints(fields[9])
            except ValueError as error:
                raise ValueError(
                    f"Invalid genePred numeric field at line "
                    f"{line_number}."
                ) from error

            if coord_mode == "1based-closed":
                tx_start -= 1
                starts = [value - 1 for value in starts]

            if (
                tx_start < 0
                or tx_end <= tx_start
                or len(starts) != exon_count
                or len(ends) != exon_count
                or not starts
            ):
                raise ValueError(
                    f"Invalid genePred structure at line {line_number}."
                )
            if any(
                start < tx_start
                or end > tx_end
                or end <= start
                for start, end in zip(starts, ends)
            ):
                raise ValueError(
                    f"Invalid exon block at genePred line "
                    f"{line_number}."
                )
            if any(
                next_start < current_end
                for current_end, next_start in zip(
                    ends[:-1],
                    starts[1:],
                )
            ):
                raise ValueError(
                    f"Overlapping or unsorted exons at genePred "
                    f"line {line_number}."
                )

            records[name] = GenePredRecord(
                name=name,
                chrom=chrom,
                strand=strand,
                tx_start=tx_start,
                tx_end=tx_end,
                exon_starts=tuple(starts),
                exon_ends=tuple(ends),
            )

    return records


def get_orf_blocks(
    row: pd.Series,
    genepred_records: dict[str, GenePredRecord],
) -> tuple[list[int], list[int], str]:
    """Resolve ORF exon blocks.

    Priority:
        1. ``exon_starts`` and ``exon_ends`` from the ORF table.
        2. genePred record keyed by ``orf_id``.
        3. single genomic interval fallback.

    Args:
        row: ORF table row.
        genepred_records: Optional genePred mapping.

    Returns:
        Starts, ends, and block source.

    Raises:
        ValueError: If provided block fields are malformed.
    """
    row_fields = (
        set(row.index)
        if hasattr(row, "index")
        else set(row.keys())
    )
    has_starts = "exon_starts" in row_fields
    has_ends = "exon_ends" in row_fields

    if has_starts or has_ends:
        if not (has_starts and has_ends):
            raise ValueError(
                f"ORF {row['orf_id']} has incomplete exon block fields."
            )
        starts = parse_comma_ints(row["exon_starts"])
        ends = parse_comma_ints(row["exon_ends"])
        if len(starts) != len(ends) or not starts:
            raise ValueError(
                f"ORF {row['orf_id']} has invalid exon blocks."
            )
        return starts, ends, "orf_table"

    orf_id = str(row["orf_id"])
    if orf_id in genepred_records:
        record = genepred_records[orf_id]
        return (
            list(record.exon_starts),
            list(record.exon_ends),
            "orf_genepred",
        )

    return (
        [int(row["genomic_start"])],
        [int(row["genomic_end"])],
        "interval",
    )


def write_output(
    table: pd.DataFrame,
    output: str | Path,
) -> None:
    """Write a tab-delimited table.

    Args:
        table: Output table.
        output: Plain or gzip-compressed output path.
    """
    output_path = str(output)
    if output_path.lower().endswith(".gz"):
        table.to_csv(
            output_path,
            sep="\t",
            index=False,
            compression="gzip",
        )
    else:
        table.to_csv(
            output_path,
            sep="\t",
            index=False,
        )


def _split_output_suffix(path: str | Path) -> tuple[str, str]:
    """Split a table output into stem and compound suffix."""
    text = str(path)
    suffixes = (
        ".tsv.gz",
        ".txt.gz",
        ".tab.gz",
        ".tsv",
        ".txt",
        ".tab",
    )
    lower = text.lower()
    for suffix in suffixes:
        if lower.endswith(suffix):
            return text[: -len(suffix)], text[-len(suffix):]
    return text, ".tsv"


def sanitize_sample_name(sample: str) -> str:
    """Convert a sample identifier to a safe file-name component.

    Args:
        sample: Density-list sample identifier.

    Returns:
        Safe sample component.
    """
    cleaned = re.sub(r"[\\/:\s]+", "_", str(sample).strip())
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", cleaned)
    cleaned = cleaned.strip("._")
    return cleaned or "sample"


def sample_output_path(
    output: str | Path,
    sample: str,
) -> str:
    """Derive one sample-level output path.

    Args:
        output: User output template.
        sample: Density-list sample identifier.

    Returns:
        Output path with ``.<sample>`` inserted before the suffix.
    """
    stem, suffix = _split_output_suffix(output)
    return f"{stem}.{sanitize_sample_name(sample)}{suffix}"


def summary_output_path(
    output: str | Path,
    sample: str | None = None,
) -> str:
    """Derive a summary output path.

    Args:
        output: User output template or sample output path.
        sample: Optional sample identifier.

    Returns:
        Sample-specific summary path.
    """
    sample_path = (
        sample_output_path(output, sample)
        if sample is not None
        else str(output)
    )
    stem, suffix = _split_output_suffix(sample_path)
    return f"{stem}.summary{suffix}"
