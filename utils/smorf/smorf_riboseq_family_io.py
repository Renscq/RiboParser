#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.007
# Function: Read smORF family inputs with high-performance coordinate recovery.
# Input: smorf_cluster outputs, source ORF table, density list, and genePred.
# Output: Validated family/member tables and sample/group density definitions.

"""Input helpers for family-aware smORF Ribo-seq evidence analysis."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TextIO

import pandas as pd

from utils.ribo.ArgsParser import message_print

from .smorf_riboseq_constants import DensityTrack

FAMILY_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "family_id",
    "orf_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "category",
    "genomic_start",
    "genomic_end",
    "nt_length",
)
MEMBER_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "family_id",
    "primary_orf_id",
    "representative_orf_id",
    "member_orf_id",
    "family_role",
    "collapse_reason",
)
ORF_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "orf_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "category",
    "genomic_start",
    "genomic_end",
    "nt_length",
    "exon_starts",
    "exon_ends",
)


@dataclass(frozen=True, slots=True)
class DensitySample:
    """Store all density tracks for one biological sample."""

    sample: str
    group: str
    tracks: tuple[DensityTrack, ...]


@dataclass(frozen=True, slots=True)
class FamilyInput:
    """Store full representative ORF records and indexed member mappings."""

    representatives: pd.DataFrame
    members: pd.DataFrame
    members_indexed: pd.DataFrame
    family_primary: pd.DataFrame


def smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain or gzip-compressed text."""
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8")
    return file_path.open(mode, encoding="utf-8")


def _normalize_columns(table: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with stripped lower-case column names."""
    table = table.copy()
    table.columns = [str(column).strip().lower() for column in table.columns]
    return table


def _require_columns(
    table: pd.DataFrame,
    required: tuple[str, ...],
    label: str,
) -> None:
    """Validate required table columns."""
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(
            f"{label} is missing column(s): " + ", ".join(missing)
        )


def _read_table(path: str | Path) -> pd.DataFrame:
    """Read one tab-delimited plain or gzip-compressed table."""
    return _normalize_columns(
        pd.read_csv(
            path,
            sep="\t",
            comment="#",
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
    )


def _parse_int_list(value: object) -> tuple[int, ...]:
    """Parse a comma-separated integer list."""
    text = str(value).strip().rstrip(",")
    if not text or text == ".":
        return ()
    return tuple(int(token) for token in text.split(",") if token)


def _validate_orf_rows(table: pd.DataFrame, label: str) -> pd.DataFrame:
    """Validate coordinates and normalize key ORF fields."""
    _require_columns(table, ORF_REQUIRED_COLUMNS, label)
    output = table.copy()
    for column in ("genomic_start", "genomic_end", "nt_length"):
        output[column] = pd.to_numeric(output[column], errors="raise").astype(
            "int64"
        )
    if output["orf_id"].duplicated().any():
        duplicated = output.loc[
            output["orf_id"].duplicated(), "orf_id"
        ].iloc[0]
        raise ValueError(f"Duplicate ORF identifier in {label}: {duplicated}")
    invalid = output["genomic_end"].le(output["genomic_start"])
    if invalid.any():
        raise ValueError(
            f"Invalid genomic interval in {label}: "
            + str(output.loc[invalid, "orf_id"].iloc[0])
        )
    return output


def _source_representatives(
    source_path: str | Path,
    representative_ids: set[str],
) -> pd.DataFrame:
    """Stream a large source ORF table with direct indexed field extraction.

    This avoids pandas chunk construction for every row of a potentially very
    large scanner message table. Only fields needed by the evidence engine are
    retained.
    """
    required = set(ORF_REQUIRED_COLUMNS)
    optional = {
        "start_codon",
        "stop_codon",
        "aa_length",
        "completeness",
        "source_strand",
        "priority",
        "frame",
    }
    selected_rows: list[dict[str, str]] = []
    found: set[str] = set()
    scanned = 0

    with smart_open(source_path, "rt") as handle:
        header: list[str] | None = None
        index: dict[str, int] = {}
        keep_columns: tuple[str, ...] = ()
        orf_index = -1
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            if header is None:
                header = [token.strip().lower() for token in raw_line.rstrip("\n").split("\t")]
                index = {column: position for position, column in enumerate(header)}
                missing = required - set(index)
                if missing:
                    raise ValueError(
                        "ORF source table is missing column(s): "
                        + ", ".join(sorted(missing))
                    )
                orf_index = index["orf_id"]
                keep_columns = tuple(
                    column
                    for column in header
                    if column in required or column in optional
                )
                continue

            fields = raw_line.rstrip("\n").split("\t")
            scanned += 1
            if orf_index >= len(fields):
                raise ValueError(
                    f"Malformed ORF source line {scanned + 1:,}: missing orf_id."
                )
            orf_id = fields[orf_index]
            if orf_id not in representative_ids:
                if scanned % 2_000_000 == 0:
                    message_print(
                        f"Scanned ORF source: rows={scanned:,}, "
                        f"representatives={len(found):,}/{len(representative_ids):,}."
                    )
                continue

            selected_rows.append(
                {
                    column: (fields[index[column]] if index[column] < len(fields) else "")
                    for column in keep_columns
                }
            )
            found.add(orf_id)
            if len(found) == len(representative_ids):
                break

    missing = representative_ids - found
    if missing:
        preview = ", ".join(sorted(missing)[:10])
        raise ValueError(
            f"ORF source is missing {len(missing):,} representative ORF(s): "
            f"{preview}"
        )
    message_print(
        f"Recovered representative coordinates: {len(found):,} ORFs "
        f"from {scanned:,} source rows."
    )
    return pd.DataFrame.from_records(selected_rows)


def read_family_input(
    family_table: str | Path,
    family_members: str | Path,
    orf_source: str | Path | None,
) -> FamilyInput:
    """Read smorf_cluster outputs and recover full representative coordinates.

    ``family.members.txt`` is intentionally compact. When it does not contain
    full exon-block coordinates, ``orf_source`` must point to the ORF message
    table used as the input of ``smorf_cluster``.
    """
    primary = _read_table(family_table)
    members = _read_table(family_members)
    _require_columns(primary, FAMILY_REQUIRED_COLUMNS, "family table")
    _require_columns(members, MEMBER_REQUIRED_COLUMNS, "family member table")

    if primary["family_id"].duplicated().any():
        family_id = primary.loc[
            primary["family_id"].duplicated(), "family_id"
        ].iloc[0]
        raise ValueError(f"Duplicate family identifier: {family_id}")
    member_families = set(members["family_id"])
    missing_family = member_families - set(primary["family_id"])
    if missing_family:
        raise ValueError(
            "Member table contains unknown family: "
            + sorted(missing_family)[0]
        )

    representative_ids = set(members["representative_orf_id"])
    representative_ids.update(primary["orf_id"])
    coordinate_columns = set(ORF_REQUIRED_COLUMNS)
    if coordinate_columns.issubset(members.columns):
        representatives = members.loc[
            members["member_orf_id"].eq(members["representative_orf_id"])
        ].copy()
        representatives = representatives.rename(
            columns={"member_orf_id": "orf_id"}
        )
        representatives = representatives.drop_duplicates("orf_id")
    else:
        multi_representative = len(representative_ids - set(primary["orf_id"]))
        if multi_representative and orf_source is None:
            raise ValueError(
                "The compact family member table lacks ORF coordinates. "
                "Provide --orf-source with the ORF message table used by "
                "smorf_cluster."
            )
        if orf_source is None:
            representatives = primary.copy()
        else:
            representatives = _source_representatives(
                source_path=orf_source,
                representative_ids=representative_ids,
            )

    representatives = _validate_orf_rows(
        representatives,
        "representative ORF table",
    )
    representatives = representatives.loc[
        representatives["orf_id"].isin(representative_ids)
    ].copy()

    rep_to_family = (
        members[["family_id", "representative_orf_id"]]
        .drop_duplicates()
        .rename(columns={"representative_orf_id": "orf_id"})
    )
    representatives = representatives.merge(
        rep_to_family,
        on="orf_id",
        how="left",
        validate="one_to_one",
    )
    primary_lookup = primary.set_index("family_id")["orf_id"]
    representatives["structural_primary"] = representatives["family_id"].map(
        primary_lookup
    ).eq(representatives["orf_id"])
    representatives["exon_starts_parsed"] = representatives[
        "exon_starts"
    ].map(_parse_int_list)
    representatives["exon_ends_parsed"] = representatives[
        "exon_ends"
    ].map(_parse_int_list)
    invalid_blocks = pd.Series(
        [
            (
                not starts
                or len(starts) != len(ends)
                or sum(end - start for start, end in zip(starts, ends))
                != int(nt_length)
            )
            for starts, ends, nt_length in zip(
                representatives["exon_starts_parsed"],
                representatives["exon_ends_parsed"],
                representatives["nt_length"],
            )
        ],
        index=representatives.index,
        dtype=bool,
    )
    if invalid_blocks.any():
        raise ValueError(
            "Invalid exon blocks for representative ORF: "
            + str(representatives.loc[invalid_blocks, "orf_id"].iloc[0])
        )

    members = members.copy()
    members["is_exact_duplicate"] = members["member_orf_id"].ne(
        members["representative_orf_id"]
    )
    members_indexed = members.set_index("family_id", drop=False)
    return FamilyInput(
        representatives=representatives,
        members=members,
        members_indexed=members_indexed,
        family_primary=primary,
    )


def _normalize_strand(value: object) -> str:
    """Normalize a density strand label."""
    text = str(value).strip().lower()
    aliases = {
        "+": "+",
        "plus": "+",
        "forward": "+",
        "fwd": "+",
        "-": "-",
        "minus": "-",
        "reverse": "-",
        "rev": "-",
        ".": ".",
        "": ".",
        "both": ".",
        "all": ".",
        "unstranded": ".",
    }
    if text not in aliases:
        raise ValueError(f"Invalid density strand: {value}")
    return aliases[text]


def _infer_strand(path: str | Path) -> str:
    """Infer plus/minus strand from density file name."""
    tokens = [
        token
        for token in re.split(r"[^a-z0-9]+", Path(path).name.lower())
        if token
    ]
    plus = any(token in {"plus", "forward", "fwd", "pos"} for token in tokens)
    minus = any(token in {"minus", "reverse", "rev", "neg"} for token in tokens)
    if plus and not minus:
        return "+"
    if minus and not plus:
        return "-"
    return "."


def read_density_design(args: object) -> tuple[DensitySample, ...]:
    """Read sample/group density definitions.

    The density list accepts ``sample``, ``path``, optional ``strand``,
    optional ``format``, and optional ``group`` columns.
    """
    density_list = getattr(args, "density_list", None)
    direct_density = list(getattr(args, "density", None) or [])
    samples: dict[str, dict[str, object]] = {}

    if density_list:
        table = _read_table(density_list)
        _require_columns(table, ("sample", "path"), "density list")
        for row_number, row in enumerate(table.to_dict("records"), start=2):
            sample = str(row["sample"]).strip()
            path = str(row["path"]).strip()
            if not sample or not path:
                raise ValueError(
                    f"Empty sample or path in density-list line {row_number}."
                )
            group = str(row.get("group", sample)).strip() or sample
            strand_raw = str(row.get("strand", "")).strip()
            strand = _normalize_strand(strand_raw) if strand_raw else _infer_strand(path)
            file_format = str(row.get("format", "auto")).strip().lower() or "auto"
            if file_format not in {"auto", "wig", "bedgraph"}:
                raise ValueError(f"Invalid density format: {file_format}")
            if not Path(path).is_file():
                raise FileNotFoundError(path)
            entry = samples.setdefault(
                sample,
                {"group": group, "tracks": []},
            )
            if entry["group"] != group:
                raise ValueError(
                    f"Sample {sample} is assigned to multiple groups."
                )
            entry["tracks"].append(
                DensityTrack(
                    sample=sample,
                    strand=strand,
                    path=path,
                    file_format=file_format,
                )
            )
    else:
        sample = str(getattr(args, "sample", "sample1")).strip()
        group = str(getattr(args, "group", sample)).strip() or sample
        file_format = str(
            getattr(args, "density_format", "auto")
        ).strip().lower()
        for path in direct_density:
            if not Path(path).is_file():
                raise FileNotFoundError(path)
            samples.setdefault(sample, {"group": group, "tracks": []})[
                "tracks"
            ].append(
                DensityTrack(
                    sample=sample,
                    strand=_infer_strand(path),
                    path=str(path),
                    file_format=file_format,
                )
            )

    if not samples:
        raise ValueError("Provide --density-list or at least one --density file.")

    output: list[DensitySample] = []
    seen_paths: set[str] = set()
    for sample, metadata in samples.items():
        tracks = tuple(metadata["tracks"])
        strands = [track.strand for track in tracks]
        if strands == ["."]:
            pass
        elif len(strands) == 2 and set(strands) == {"+", "-"}:
            pass
        else:
            raise ValueError(
                f"Sample {sample} requires exactly one unstranded track or "
                "exactly one plus/minus pair."
            )
        for track in tracks:
            canonical = str(Path(track.path).resolve())
            if canonical in seen_paths:
                raise ValueError(f"Density file is listed more than once: {track.path}")
            seen_paths.add(canonical)
        output.append(
            DensitySample(
                sample=sample,
                group=str(metadata["group"]),
                tracks=tracks,
            )
        )
    return tuple(output)


def resolve_analysis_mode(requested: str, sample_count: int) -> str:
    """Resolve automatic analysis mode."""
    mode = str(requested).strip().lower()
    aliases = {
        "individual": "separate",
        "sample": "separate",
        "merge": "pooled",
        "merged": "pooled",
        "single": "single",
    }
    mode = aliases.get(mode, mode)
    if mode == "auto":
        return "single" if sample_count == 1 else "hybrid"
    if mode not in {"single", "separate", "pooled", "hybrid"}:
        raise ValueError(f"Unknown analysis mode: {requested}")
    if sample_count == 1:
        return "single"
    if mode == "single":
        raise ValueError("single mode requires exactly one sample.")
    return mode
