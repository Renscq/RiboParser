#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-31
# Version: 0.2.8.24-dev.001
# Function: Quantify reliable smORF P-site density from genomic density tracks.
# Input: Reliable smORF genePred and a sample/strand density design table.
# Output: Wide smORF-by-sample raw P-site density count matrix.

"""High-performance smORF P-site density quantification.

Workflow
--------
1. Parse reliable smORF genePred or genePredExt records.
2. Build transcript-oriented CDS blocks and optionally remove terminal stops.
3. Validate the multi-sample density design and strand assignments.
4. Parse each WIG/bedGraph track once and quantify genomic overlaps directly.
5. Write one atomic, input-order-preserving smORF count matrix.

The implementation quantifies sparse density intervals without constructing a
dense nucleotide profile for every ORF. Frame assignment is calculated from
the spliced transcript offset, so it remains correct across exon junctions and
on the minus strand.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TextIO

import numpy as np

from utils.ribo.ArgsParser import progress_print

from .smorf_riboseq_density import (
    ChromDensity,
    infer_density_format,
    iter_density_chromosomes,
    prepare_density,
)


OUTPUT_SUFFIX: Final[str] = ".density_quant.txt"
METADATA_COLUMNS: Final[tuple[str, ...]] = (
    "orf_id",
    "gene_id",
    "chrom",
    "strand",
    "tx_start",
    "tx_end",
    "cds_start",
    "cds_end",
    "exon_count",
    "coding_nt_length",
    "coding_codon_count",
)
VALID_FRAMES: Final[frozenset[str]] = frozenset({"all", "0", "1", "2"})
_BUFFER_BYTES: Final[int] = 8 * 1024 * 1024
_WORKER_RECORDS_BY_CHROM: Mapping[
    str,
    Mapping[str, tuple["QuantORF", ...]],
] | None = None
_WORKER_RECORD_COUNT: int = 0


@dataclass(frozen=True, slots=True)
class QuantBlock:
    """Store one quantified genomic CDS block and its transcript offset."""

    start: int
    end: int
    transcript_offset: int


@dataclass(frozen=True, slots=True)
class QuantORF:
    """Store compact ORF geometry required for density quantification."""

    order: int
    orf_id: str
    gene_id: str
    chrom: str
    strand: str
    tx_start: int
    tx_end: int
    cds_start: int
    cds_end: int
    exon_count: int
    coding_nt_length: int
    blocks: tuple[QuantBlock, ...]

    @property
    def coding_codon_count(self) -> int:
        """Return the quantified CDS length in codons."""
        return self.coding_nt_length // 3


@dataclass(frozen=True, slots=True)
class DensityTrack:
    """Describe one strand-specific or unstranded density track."""

    sample: str
    strand: str
    path: str
    file_format: str


@dataclass(frozen=True, slots=True)
class SampleTracks:
    """Store all density tracks assigned to one biological sample."""

    sample: str
    tracks: tuple[DensityTrack, ...]


@dataclass(frozen=True, slots=True)
class QuantConfig:
    """Store public smORF quantification settings."""

    genepred: str
    density_list: str
    output_prefix: str
    frame: str = "all"
    include_stop: bool = False
    threads: int = 1


@dataclass(frozen=True, slots=True)
class SampleQuantResult:
    """Store one sample's count vector and chromosome validation state."""

    sample: str
    counts: np.ndarray
    observed_target_chromosomes: int


@dataclass(frozen=True, slots=True)
class QuantResult:
    """Store final smORF quantification outputs and summary values."""

    output: str
    orf_count: int
    sample_count: int
    effective_workers: int
    frame: str
    include_stop: bool


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
        raise ValueError(
            f"Invalid {column} at genePred line {line_number}."
        ) from error


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
                raise ValueError(
                    f"genePred line {line_number} has fewer than 10 columns."
                )

            orf_id = fields[0].strip()
            chrom = fields[1].strip()
            strand = fields[2].strip()
            if not orf_id or not chrom:
                raise ValueError(
                    f"Empty ORF ID or chromosome at genePred line {line_number}."
                )
            if orf_id in seen_ids:
                raise ValueError(f"Duplicate ORF ID in genePred: {orf_id}")
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid strand at genePred line {line_number}: {strand}"
                )

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

            exon_starts = _parse_int_list(
                fields[8], line_number, "exonStarts"
            )
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
                raise ValueError(
                    f"Invalid genePred structure at line {line_number}: {orf_id}"
                )
            previous_end = -1
            for start, end in zip(exon_starts, exon_ends):
                if (
                    start < tx_start
                    or end > tx_end
                    or end <= start
                    or start < previous_end
                ):
                    raise ValueError(
                        f"Invalid exon structure at genePred line {line_number}: "
                        f"{orf_id}"
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

            remove_stop = (
                not include_stop and _terminal_stop_is_complete(fields)
            )
            coding_length = raw_coding_length - (3 if remove_stop else 0)
            if coding_length <= 0 or coding_length % 3 != 0:
                raise ValueError(
                    f"Quantified CDS length is invalid for {orf_id}: "
                    f"{coding_length}"
                )

            gene_id = (
                fields[11].strip()
                if len(fields) >= 12 and fields[11].strip()
                else orf_id
            )
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
        raise ValueError(
            f"Invalid density strand at line {line_number}: {value}"
        )
    return aliases[text]


def read_density_design(path: str | Path) -> tuple[SampleTracks, ...]:
    """Read a header-based density design in stable input order."""
    design_path = Path(path).expanduser().resolve()
    with _smart_open(design_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original_header = tuple(reader.fieldnames or ())
        if not original_header:
            raise ValueError(
                "Density list requires a tab-delimited header."
            )
        header = {
            str(name).strip().lower(): str(name)
            for name in original_header
        }
        sample_column = next(
            (header[name] for name in ("sample", "name") if name in header),
            None,
        )
        path_column = next(
            (
                header[name]
                for name in ("path", "file", "density")
                if name in header
            ),
            None,
        )
        strand_column = header.get("strand")
        format_column = header.get("format")
        if sample_column is None or path_column is None or strand_column is None:
            raise ValueError(
                "Density list requires sample/name, strand, and "
                "path/file/density columns."
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
                    f"Density-list line {line_number} contains an empty "
                    "sample or path."
                )
            if sample in METADATA_COLUMNS:
                raise ValueError(
                    f"Density sample name conflicts with an output column: {sample}"
                )

            density_path = Path(density_text).expanduser()
            if not density_path.is_absolute():
                density_path = design_path.parent / density_path
            density_path = density_path.resolve()
            if not density_path.is_file():
                raise FileNotFoundError(density_path)
            canonical = str(density_path)
            if canonical in seen_paths:
                raise ValueError(
                    f"Density file is listed more than once: {density_path}"
                )
            pair = (sample, strand)
            if pair in seen_pairs:
                raise ValueError(
                    f"Duplicate density sample/strand: {sample}/{strand}"
                )

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
                f"Sample {sample} requires one unstranded track or a complete "
                "plus/minus pair."
            )
        output.append(SampleTracks(sample=sample, tracks=tracks))
    return tuple(output)


def _count_positions_in_frame(start: int, end: int, frame: int) -> int:
    """Count transcript positions in ``[start, end)`` assigned to one frame."""
    return (end + 2 - frame) // 3 - (start + 2 - frame) // 3


def _quantify_orf(
    density: ChromDensity,
    record: QuantORF,
    frame: str,
) -> float:
    """Quantify one ORF directly from sparse density intervals."""
    total = 0.0
    selected_frame = None if frame == "all" else int(frame)

    for block in record.blocks:
        left = int(np.searchsorted(density.ends, block.start, side="right"))
        right = int(np.searchsorted(density.starts, block.end, side="left"))
        for index in range(left, right):
            value = abs(float(density.values[index]))
            if value <= 0.0 or not math.isfinite(value):
                continue
            overlap_start = max(block.start, int(density.starts[index]))
            overlap_end = min(block.end, int(density.ends[index]))
            if overlap_end <= overlap_start:
                continue

            if selected_frame is None:
                position_count = overlap_end - overlap_start
            else:
                if record.strand == "+":
                    tx_start = (
                        block.transcript_offset
                        + overlap_start
                        - block.start
                    )
                    tx_end = (
                        block.transcript_offset
                        + overlap_end
                        - block.start
                    )
                else:
                    tx_start = (
                        block.transcript_offset
                        + block.end
                        - overlap_end
                    )
                    tx_end = (
                        block.transcript_offset
                        + block.end
                        - overlap_start
                    )
                position_count = _count_positions_in_frame(
                    tx_start,
                    tx_end,
                    selected_frame,
                )
            total += value * position_count
    return total


def _records_by_chromosome(
    records: Sequence[QuantORF],
) -> dict[str, dict[str, tuple[QuantORF, ...]]]:
    """Group ORFs by chromosome and strand."""
    temporary: dict[str, dict[str, list[QuantORF]]] = defaultdict(
        lambda: {"+": [], "-": []}
    )
    for record in records:
        temporary[record.chrom][record.strand].append(record)
    return {
        chrom: {
            strand: tuple(strand_records)
            for strand, strand_records in mapping.items()
        }
        for chrom, mapping in temporary.items()
    }


def _quantify_sample(
    sample_tracks: SampleTracks,
    records_by_chrom: Mapping[str, Mapping[str, tuple[QuantORF, ...]]],
    record_count: int,
    frame: str,
    work_directory: str,
) -> SampleQuantResult:
    """Quantify one sample while parsing each assigned track once."""
    counts = np.zeros(record_count, dtype=np.float64)
    observed_target_chromosomes: set[str] = set()

    for track_number, track in enumerate(sample_tracks.tracks, start=1):
        track_work = (
            Path(work_directory)
            / f"track_{track_number}_{track.strand.replace('.', 'all')}"
        )
        prepared = prepare_density(
            track.path,
            track.file_format,
            track_work,
        )
        try:
            iterator_format = (
                "bedgraph" if prepared.was_sorted else track.file_format
            )
            for density in iter_density_chromosomes(
                prepared.path,
                iterator_format,
            ):
                chrom_records = records_by_chrom.get(density.chrom)
                if chrom_records is None:
                    continue
                observed_target_chromosomes.add(density.chrom)
                strands = ("+", "-") if track.strand == "." else (track.strand,)
                for strand in strands:
                    for record in chrom_records[strand]:
                        counts[record.order] = _quantify_orf(
                            density,
                            record,
                            frame,
                        )
        finally:
            prepared.cleanup()

    return SampleQuantResult(
        sample=sample_tracks.sample,
        counts=counts,
        observed_target_chromosomes=len(observed_target_chromosomes),
    )


def _initialize_sample_worker(
    records_by_chrom: Mapping[str, Mapping[str, tuple[QuantORF, ...]]],
    record_count: int,
) -> None:
    """Initialize one worker with the read-only ORF geometry index."""
    global _WORKER_RECORDS_BY_CHROM, _WORKER_RECORD_COUNT
    _WORKER_RECORDS_BY_CHROM = records_by_chrom
    _WORKER_RECORD_COUNT = int(record_count)


def _sample_worker(
    sample_tracks: SampleTracks,
    frame: str,
    work_directory: str,
) -> SampleQuantResult:
    """Process-pool entry for one sample."""
    if _WORKER_RECORDS_BY_CHROM is None or _WORKER_RECORD_COUNT < 1:
        raise RuntimeError("smorf_quant worker was not initialized.")
    return _quantify_sample(
        sample_tracks,
        _WORKER_RECORDS_BY_CHROM,
        _WORKER_RECORD_COUNT,
        frame,
        work_directory,
    )


def _format_density(value: float) -> str:
    """Format one non-negative finite density value compactly."""
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Invalid quantified density value: {value}")
    rounded = round(value)
    if math.isclose(value, rounded, rel_tol=0.0, abs_tol=1e-10):
        return str(int(rounded))
    return format(value, ".10g")


def output_path_from_prefix(prefix: str | Path) -> Path:
    """Derive the standard count-matrix path from an output prefix."""
    return Path(f"{prefix}{OUTPUT_SUFFIX}")


def _write_matrix(
    output_path: Path,
    records: Sequence[QuantORF],
    sample_order: Sequence[str],
    counts_by_sample: Mapping[str, np.ndarray],
) -> None:
    """Write a complete wide count matrix atomically."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + ".tmp")
    try:
        with temporary.open(
            "w",
            encoding="utf-8",
            newline="",
            buffering=_BUFFER_BYTES,
        ) as handle:
            handle.write("\t".join((*METADATA_COLUMNS, *sample_order)) + "\n")
            for record in records:
                row = [
                    record.orf_id,
                    record.gene_id,
                    record.chrom,
                    record.strand,
                    str(record.tx_start),
                    str(record.tx_end),
                    str(record.cds_start),
                    str(record.cds_end),
                    str(record.exon_count),
                    str(record.coding_nt_length),
                    str(record.coding_codon_count),
                ]
                row.extend(
                    _format_density(counts_by_sample[sample][record.order])
                    for sample in sample_order
                )
                handle.write("\t".join(row) + "\n")
        os.replace(temporary, output_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


class SmORFQuantifier:
    """Quantify reliable smORFs across genomic P-site density samples."""

    def __init__(self, config: QuantConfig) -> None:
        """Initialize the quantifier and validate scalar configuration."""
        frame = str(config.frame).lower()
        if frame not in VALID_FRAMES:
            raise ValueError(
                f"Unsupported frame: {config.frame}. Use all, 0, 1, or 2."
            )
        if int(config.threads) < 1:
            raise ValueError("threads must be >= 1.")
        self.config = QuantConfig(
            genepred=str(config.genepred),
            density_list=str(config.density_list),
            output_prefix=str(config.output_prefix),
            frame=frame,
            include_stop=bool(config.include_stop),
            threads=int(config.threads),
        )

    def run(self) -> QuantResult:
        """Run annotation parsing, sample quantification, and atomic output."""
        records = read_smorf_genepred(
            self.config.genepred,
            include_stop=self.config.include_stop,
        )
        samples = read_density_design(self.config.density_list)
        records_by_chrom = _records_by_chromosome(records)
        sample_order = tuple(sample.sample for sample in samples)
        workers = max(1, min(self.config.threads, len(samples)))
        counts_by_sample: dict[str, np.ndarray] = {}

        with tempfile.TemporaryDirectory(prefix="smorf_quant.") as work_root:
            if workers == 1:
                for number, sample in enumerate(samples, start=1):
                    result = _quantify_sample(
                        sample,
                        records_by_chrom,
                        len(records),
                        self.config.frame,
                        str(Path(work_root) / f"sample_{number}"),
                    )
                    if result.observed_target_chromosomes == 0:
                        raise ValueError(
                            f"No density chromosome matched the smORF genePred "
                            f"for sample {result.sample}."
                        )
                    counts_by_sample[result.sample] = result.counts
                    progress_print(
                        f"quantified sample {number:,}/{len(samples):,}: "
                        f"{result.sample}"
                    )
            else:
                with ProcessPoolExecutor(
                    max_workers=workers,
                    initializer=_initialize_sample_worker,
                    initargs=(records_by_chrom, len(records)),
                ) as executor:
                    futures = {
                        executor.submit(
                            _sample_worker,
                            sample,
                            self.config.frame,
                            str(Path(work_root) / f"sample_{number}"),
                        ): sample.sample
                        for number, sample in enumerate(samples, start=1)
                    }
                    completed = 0
                    for future in as_completed(futures):
                        result = future.result()
                        if result.observed_target_chromosomes == 0:
                            raise ValueError(
                                f"No density chromosome matched the smORF "
                                f"genePred for sample {result.sample}."
                            )
                        counts_by_sample[result.sample] = result.counts
                        completed += 1
                        progress_print(
                            f"quantified sample {completed:,}/{len(samples):,}: "
                            f"{result.sample}"
                        )

        output_path = output_path_from_prefix(self.config.output_prefix)
        _write_matrix(
            output_path,
            records,
            sample_order,
            counts_by_sample,
        )
        return QuantResult(
            output=str(output_path),
            orf_count=len(records),
            sample_count=len(samples),
            effective_workers=workers,
            frame=self.config.frame,
            include_stop=self.config.include_stop,
        )


def run_smorf_quant(args: object) -> QuantResult:
    """Build a quantification config from argparse arguments and run it."""
    quantifier = SmORFQuantifier(
        QuantConfig(
            genepred=str(getattr(args, "genepred")),
            density_list=str(getattr(args, "density_list")),
            output_prefix=str(getattr(args, "output")),
            frame=str(getattr(args, "frame", "all")),
            include_stop=bool(getattr(args, "include_stop", False)),
            threads=int(getattr(args, "threads", 1)),
        )
    )
    return quantifier.run()
