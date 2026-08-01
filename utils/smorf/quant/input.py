#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Read reliable smORFs and density designs.
# Input: smORF genePred and density-list files.
# Output: Validated ORF records and sample tracks.

"""Read reliable smORFs and density designs."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from utils.smorf.density import (
    infer_density_format,
)

from .models import (
    _BUFFER_BYTES,
    METADATA_COLUMNS,
    DensityTrack,
    QuantBlock,
    QuantORF,
    SampleTracks,
)


def _smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain or gzip-compressed text."""
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(
            file_path,
            mode,
            encoding=None if "b" in mode else "utf-8",
            newline=None if "b" in mode else "",
        )
    return file_path.open(
        mode,
        encoding=None if "b" in mode else "utf-8",
        newline=None if "b" in mode else "",
        buffering=_BUFFER_BYTES,
    )


def _parse_int_list(value: str, line_number: int, column: str) -> tuple[int, ...]:
    """Parse a comma-separated genePred integer list."""
    text = str(value).strip().rstrip(",")
    if not text:
        return ()
    try:
        return tuple(int(item) for item in text.split(",") if item)
    except ValueError as error:
        raise ValueError(f"Invalid {column} at genePred line {line_number}.") from error


def _intersect_cds_blocks(
    exon_starts: Sequence[int],
    exon_ends: Sequence[int],
    cds_start: int,
    cds_end: int,
) -> list[tuple[int, int]]:
    """Return genomic exon intersections with the CDS interval."""
    return [
        (max(start, cds_start), min(end, cds_end))
        for start, end in zip(exon_starts, exon_ends)
        if min(end, cds_end) > max(start, cds_start)
    ]


def _trim_oriented_blocks(
    blocks: Sequence[tuple[int, int]],
    strand: str,
    retained_length: int,
) -> tuple[QuantBlock, ...]:
    """Trim CDS blocks at the transcript 3-prime end."""
    oriented = list(blocks if strand == "+" else reversed(blocks))
    remaining = int(retained_length)
    transcript_offset = 0
    output: list[QuantBlock] = []

    for start, end in oriented:
        if remaining <= 0:
            break
        block_length = end - start
        take = min(block_length, remaining)
        if strand == "+":
            kept_start = start
            kept_end = start + take
        else:
            kept_start = end - take
            kept_end = end
        output.append(
            QuantBlock(
                start=kept_start,
                end=kept_end,
                transcript_offset=transcript_offset,
            )
        )
        transcript_offset += take
        remaining -= take

    if remaining != 0 or transcript_offset != retained_length:
        raise ValueError("Failed to trim CDS blocks to the requested length.")
    return tuple(output)


def _terminal_stop_is_complete(fields: Sequence[str]) -> bool:
    """Return whether a genePred record has a complete terminal stop."""
    if len(fields) >= 14:
        status = fields[13].strip().lower()
        if status:
            return status in {"cmpl", "complete", "full"}
    return True


def read_smorf_genepred(
    path: str | Path,
    *,
    include_stop: bool = False,
) -> tuple[QuantORF, ...]:
    """Read and validate reliable smORF genePred records.

    The quantified region is the spliced CDS. A complete terminal stop codon
    is removed unless ``include_stop`` is enabled.
    """
    records: list[QuantORF] = []
    seen_ids: set[str] = set()

    with _smart_open(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            text = raw_line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if len(fields) < 10:
                raise ValueError(f"genePred line {line_number} has fewer than 10 columns.")

            orf_id = fields[0].strip()
            chrom = fields[1].strip()
            strand = fields[2].strip()
            if not orf_id or not chrom:
                raise ValueError(f"Empty ORF ID or chromosome at genePred line {line_number}.")
            if orf_id in seen_ids:
                raise ValueError(f"Duplicate ORF ID in genePred: {orf_id}")
            if strand not in {"+", "-"}:
                raise ValueError(f"Invalid strand at genePred line {line_number}: {strand}")

            try:
                tx_start = int(fields[3])
                tx_end = int(fields[4])
                cds_start = int(fields[5])
                cds_end = int(fields[6])
                exon_count = int(fields[7])
            except ValueError as error:
                raise ValueError(
                    f"Invalid numeric field at genePred line {line_number}."
                ) from error

            exon_starts = _parse_int_list(fields[8], line_number, "exonStarts")
            exon_ends = _parse_int_list(fields[9], line_number, "exonEnds")
            if (
                tx_start < 0
                or tx_end <= tx_start
                or cds_start < tx_start
                or cds_end > tx_end
                or cds_end <= cds_start
                or exon_count < 1
                or len(exon_starts) != exon_count
                or len(exon_ends) != exon_count
            ):
                raise ValueError(f"Invalid genePred structure at line {line_number}: {orf_id}")
            previous_end = -1
            for start, end in zip(exon_starts, exon_ends):
                if start < tx_start or end > tx_end or end <= start or start < previous_end:
                    raise ValueError(
                        f"Invalid exon structure at genePred line {line_number}: {orf_id}"
                    )
                previous_end = end

            cds_blocks = _intersect_cds_blocks(
                exon_starts,
                exon_ends,
                cds_start,
                cds_end,
            )
            raw_coding_length = sum(end - start for start, end in cds_blocks)
            if raw_coding_length <= 0 or raw_coding_length % 3 != 0:
                raise ValueError(
                    f"Spliced CDS length is not a positive multiple of three "
                    f"for {orf_id}: {raw_coding_length}"
                )

            remove_stop = not include_stop and _terminal_stop_is_complete(fields)
            coding_length = raw_coding_length - (3 if remove_stop else 0)
            if coding_length <= 0 or coding_length % 3 != 0:
                raise ValueError(f"Quantified CDS length is invalid for {orf_id}: {coding_length}")

            gene_id = fields[11].strip() if len(fields) >= 12 and fields[11].strip() else orf_id
            records.append(
                QuantORF(
                    order=len(records),
                    orf_id=orf_id,
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    tx_start=tx_start,
                    tx_end=tx_end,
                    cds_start=cds_start,
                    cds_end=cds_end,
                    exon_count=exon_count,
                    coding_nt_length=coding_length,
                    blocks=_trim_oriented_blocks(
                        cds_blocks,
                        strand,
                        coding_length,
                    ),
                )
            )
            seen_ids.add(orf_id)

    if not records:
        raise ValueError("The smORF genePred contains no quantifiable records.")
    return tuple(records)


def _normalize_density_strand(value: str, line_number: int) -> str:
    """Normalize one density-list strand label."""
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
        "both": ".",
        "all": ".",
        "unstranded": ".",
    }
    if text not in aliases:
        raise ValueError(f"Invalid density strand at line {line_number}: {value}")
    return aliases[text]


def read_density_design(path: str | Path) -> tuple[SampleTracks, ...]:
    """Read a header-based density design in stable input order."""
    design_path = Path(path).expanduser().resolve()
    with _smart_open(design_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original_header = tuple(reader.fieldnames or ())
        if not original_header:
            raise ValueError("Density list requires a tab-delimited header.")
        header = {str(name).strip().lower(): str(name) for name in original_header}
        sample_column = next(
            (header[name] for name in ("sample", "name") if name in header),
            None,
        )
        path_column = next(
            (header[name] for name in ("path", "file", "density") if name in header),
            None,
        )
        strand_column = header.get("strand")
        format_column = header.get("format")
        if sample_column is None or path_column is None or strand_column is None:
            raise ValueError(
                "Density list requires sample/name, strand, and path/file/density columns."
            )

        sample_order: list[str] = []
        tracks_by_sample: dict[str, list[DensityTrack]] = defaultdict(list)
        seen_paths: set[str] = set()
        seen_pairs: set[tuple[str, str]] = set()

        for line_number, row in enumerate(reader, start=2):
            sample = str(row.get(sample_column, "")).strip()
            density_text = str(row.get(path_column, "")).strip()
            strand = _normalize_density_strand(
                str(row.get(strand_column, "")),
                line_number,
            )
            if not sample or not density_text:
                raise ValueError(
                    f"Density-list line {line_number} contains an empty sample or path."
                )
            if sample in METADATA_COLUMNS:
                raise ValueError(f"Density sample name conflicts with an output column: {sample}")

            density_path = Path(density_text).expanduser()
            if not density_path.is_absolute():
                density_path = design_path.parent / density_path
            density_path = density_path.resolve()
            if not density_path.is_file():
                raise FileNotFoundError(density_path)
            canonical = str(density_path)
            if canonical in seen_paths:
                raise ValueError(f"Density file is listed more than once: {density_path}")
            pair = (sample, strand)
            if pair in seen_pairs:
                raise ValueError(f"Duplicate density sample/strand: {sample}/{strand}")

            user_format = (
                str(row.get(format_column, "auto")).strip().lower()
                if format_column is not None
                else "auto"
            )
            file_format = infer_density_format(
                density_path,
                user_format=user_format or "auto",
            )
            if sample not in tracks_by_sample:
                sample_order.append(sample)
            tracks_by_sample[sample].append(
                DensityTrack(
                    sample=sample,
                    strand=strand,
                    path=canonical,
                    file_format=file_format,
                )
            )
            seen_paths.add(canonical)
            seen_pairs.add(pair)

    if not sample_order:
        raise ValueError("Density list contains no data rows.")

    output: list[SampleTracks] = []
    for sample in sample_order:
        tracks = tuple(tracks_by_sample[sample])
        strands = {track.strand for track in tracks}
        if strands not in ({"."}, {"+", "-"}):
            raise ValueError(
                f"Sample {sample} requires one unstranded track or a complete plus/minus pair."
            )
        output.append(SampleTracks(sample=sample, tracks=tracks))
    return tuple(output)
