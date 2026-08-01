#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Count frame-aware P-sites for reliable smORFs.
# Input: ORF records and sparse density tracks.
# Output: Per-sample raw count vectors.

"""Count frame-aware P-sites for reliable smORFs."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from utils.smorf.density import (
    ChromDensity,
    iter_density_chromosomes,
    prepare_density,
)

from .models import (
    QuantORF,
    SampleQuantResult,
    SampleTracks,
)

_WORKER_RECORDS_BY_CHROM: (
    Mapping[
        str,
        Mapping[str, tuple[QuantORF, ...]],
    ]
    | None
) = None

_WORKER_RECORD_COUNT: int = 0


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
                    tx_start = block.transcript_offset + overlap_start - block.start
                    tx_end = block.transcript_offset + overlap_end - block.start
                else:
                    tx_start = block.transcript_offset + block.end - overlap_end
                    tx_end = block.transcript_offset + block.end - overlap_start
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
    temporary: dict[str, dict[str, list[QuantORF]]] = defaultdict(lambda: {"+": [], "-": []})
    for record in records:
        temporary[record.chrom][record.strand].append(record)
    return {
        chrom: {strand: tuple(strand_records) for strand, strand_records in mapping.items()}
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
            Path(work_directory) / f"track_{track_number}_{track.strand.replace('.', 'all')}"
        )
        prepared = prepare_density(
            track.path,
            track.file_format,
            track_work,
        )
        try:
            iterator_format = "bedgraph" if prepared.was_sorted else track.file_format
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
