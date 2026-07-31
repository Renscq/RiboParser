#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.001
# Function: Evaluate smORF Ribo-seq evidence with sparse profiles and resumable output.
# Input: smorf_cluster family table, P-site density tracks, and optional genePred.
# Output: Compact sample evidence tables, summaries, and resumable checkpoints.

"""Fast sample-level smORF Ribo-seq evidence analysis.

The implementation uses five performance rules:

1. Compile ORF genomic geometry once before processing samples.
2. Map sparse P-site intervals directly to transcript codons and frames.
3. Apply the abundance/coverage core gate before expensive shape and release
   calculations.
4. Write compact evidence rows in buffered batches.
5. Checkpoint after every completed chromosome so interrupted runs can resume.
"""

from __future__ import annotations

import gzip
import json
import math
import multiprocessing
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np

from utils.ribo.ArgsParser import message_print, warning_print

from .smorf_riboseq_constants import (
    DensityTrack,
    EVIDENCE_LEVELS,
    EvidenceThresholds,
    STOP_CODONS,
)
from .smorf_riboseq_density import (
    ChromDensity,
    iter_density_chromosomes,
    prepare_density,
)
from .smorf_riboseq_io import (
    GenePredRecord,
    get_orf_blocks,
    read_density_list,
    read_genepred,
    read_orf_table,
    sample_output_path,
)
from .smorf_riboseq_metrics import (
    calculate_pausing,
    classify_translation_evidence,
    safe_divide,
)


COMPACT_OUTPUT_COLUMNS = (
    "sample",
    "density_strand",
    "orf_id",
    "family_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "category",
    "genomic_start",
    "genomic_end",
    "nt_length",
    "coding_nt_length",
    "coding_codon_count",
    "profile_status",
    "rpf_sum",
    "rpf_per_codon",
    "covered_codon",
    "covered_codon_ratio",
    "max_density",
    "frame0_density",
    "frame1_density",
    "frame2_density",
    "frame0_ratio",
    "frame1_ratio",
    "frame2_ratio",
    "frame0_vs_alt_ratio",
    "periodicity_score",
    "periodicity_label",
    "release_ratio",
    "release_label",
    "coverage_shape",
    "evidence_components",
    "evidence_reason",
    "translation_evidence",
)

FULL_OUTPUT_COLUMNS = (
    "sample",
    "density_strand",
    "orf_id",
    "family_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "category",
    "genomic_start",
    "genomic_end",
    "nt_length",
    "coding_nt_length",
    "coding_codon_count",
    "block_source",
    "profile_status",
    "profile_reason",
    "rpf_sum",
    "rpf_mean",
    "rpf_per_codon",
    "covered_nt",
    "coverage_ratio",
    "covered_codon",
    "covered_codon_ratio",
    "max_density",
    "frame0_density",
    "frame1_density",
    "frame2_density",
    "frame0_ratio",
    "frame1_ratio",
    "frame2_ratio",
    "frame0_vs_alt_ratio",
    "periodicity_score",
    "periodicity_label",
    "periodicity_evaluable",
    "start_codon_mean",
    "body_codon_mean",
    "pre_stop_codon_mean",
    "start_pause_ratio",
    "stop_pause_ratio",
    "start_pause_label",
    "stop_pause_label",
    "pausing_label",
    "pausing_evaluable",
    "post_stop_mean",
    "release_ratio",
    "release_drop_score",
    "release_label",
    "release_evaluable",
    "release_context",
    "post_stop_nt",
    "codon_gini",
    "max_to_mean_ratio",
    "top10_fraction",
    "coverage_shape",
    "shape_evaluable",
    "evidence_components",
    "evidence_reason",
    "translation_evidence",
)

CHECKPOINT_VERSION = 1
WRITE_BATCH_SIZE = 2048
CHECKPOINT_ORF_INTERVAL = 50_000
COPY_BUFFER_SIZE = 16 * 1024 * 1024


@dataclass(frozen=True, slots=True)
class GeometryBlock:
    """Store one genomic block and its transcript-oriented offset."""

    start: int
    end: int
    offset: int


@dataclass(frozen=True, slots=True)
class ORFGeometry:
    """Store sample-independent ORF metadata and compiled genomic geometry."""

    input_index: int
    orf_id: str
    family_id: str
    gene_id: str
    transcript_id: str
    chrom: str
    strand: str
    category: str
    genomic_start: int
    genomic_end: int
    nt_length: int
    coding_nt_length: int
    coding_codon_count: int
    block_source: str
    blocks: tuple[GeometryBlock, ...]
    downstream_blocks: tuple[GeometryBlock, ...]
    release_context: str


@dataclass(frozen=True, slots=True)
class GeometryIndex:
    """Store compiled ORFs and chromosome/strand lookup indices."""

    geometries: tuple[ORFGeometry, ...]
    all_by_chrom: Mapping[str, tuple[int, ...]]
    plus_by_chrom: Mapping[str, tuple[int, ...]]
    minus_by_chrom: Mapping[str, tuple[int, ...]]


@dataclass(frozen=True, slots=True)
class SparseCoreProfile:
    """Store core abundance, coverage, frame, and codon-profile metrics."""

    codon_profile: np.ndarray
    rpf_sum: float
    covered_nt: int
    covered_codon: int
    max_density: float
    frame0_density: float
    frame1_density: float
    frame2_density: float


@dataclass(frozen=True, slots=True)
class SampleEvidenceResult:
    """Store one sample's output and evidence counts."""

    sample: str
    output: str
    written_rows: int
    total_orfs: int
    no_evidence: int
    low_confidence: int
    medium_confidence: int
    high_confidence: int
    sorted_track_count: int
    resumed: bool
    partial_output: str


@dataclass(frozen=True, slots=True)
class EvidenceRunResult:
    """Store all sample-specific output information."""

    samples: tuple[SampleEvidenceResult, ...]
    threads: int


@dataclass(slots=True)
class EvidenceCounters:
    """Store cumulative evidence counts for checkpointing."""

    processed: int = 0
    written: int = 0
    no_evidence: int = 0
    low_confidence: int = 0
    medium_confidence: int = 0
    high_confidence: int = 0
    invalid_profile: int = 0

    def add_label(self, label: str) -> None:
        """Increment one evidence label."""
        self.processed += 1
        if label == "NoEvidence":
            self.no_evidence += 1
        elif label == "LowConfidence":
            self.low_confidence += 1
        elif label == "MediumConfidence":
            self.medium_confidence += 1
        elif label == "HighConfidence":
            self.high_confidence += 1
        else:
            self.invalid_profile += 1

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> "EvidenceCounters":
        """Build counters from checkpoint JSON."""
        return cls(
            processed=int(values.get("processed", 0)),
            written=int(values.get("written", 0)),
            no_evidence=int(values.get("no_evidence", 0)),
            low_confidence=int(values.get("low_confidence", 0)),
            medium_confidence=int(values.get("medium_confidence", 0)),
            high_confidence=int(values.get("high_confidence", 0)),
            invalid_profile=int(values.get("invalid_profile", 0)),
        )


def build_thresholds(args: object) -> EvidenceThresholds:
    """Build thresholds from a preset and optional overrides."""
    return EvidenceThresholds.from_mode(
        getattr(args, "evidence_mode", "balanced")
    ).with_overrides(
        min_rpf_sum=getattr(args, "min_rpf_sum", None),
        min_covered_codon=getattr(args, "min_covered_codon", None),
        min_codon_coverage=getattr(args, "min_codon_coverage", None),
        moderate_periodicity=getattr(args, "moderate_periodicity", None),
        strong_periodicity=getattr(args, "strong_periodicity", None),
    )


def _parse_blocks(
    row: Mapping[str, object],
    genepred_records: dict[str, GenePredRecord],
) -> tuple[list[int], list[int], str]:
    """Resolve and validate ORF blocks once."""
    starts, ends, source = get_orf_blocks(row, genepred_records)
    if len(starts) != len(ends) or not starts:
        raise ValueError(f"ORF {row['orf_id']} has invalid exon blocks.")

    blocks = sorted(
        ((int(start), int(end)) for start, end in zip(starts, ends)),
        key=lambda item: item[0],
    )
    previous_end: int | None = None
    for start, end in blocks:
        if start < 0 or end <= start:
            raise ValueError(
                f"ORF {row['orf_id']} has invalid block [{start}, {end})."
            )
        if previous_end is not None and start < previous_end:
            raise ValueError(
                f"ORF {row['orf_id']} has overlapping exon blocks."
            )
        previous_end = end
    return (
        [start for start, _ in blocks],
        [end for _, end in blocks],
        source,
    )


def _coding_nt_length(row: Mapping[str, object], profile_length: int) -> int:
    """Calculate coding length excluding a complete terminal stop codon."""
    declared_length = min(int(row["nt_length"]), int(profile_length))
    stop_codon = (
        str(row.get("stop_codon", ""))
        .strip()
        .upper()
        .replace("U", "T")
    )
    completeness = str(row.get("completeness", "")).strip().lower()
    complete_labels = {"", "complete", "cmpl", "full", "true", "yes"}
    if (
        stop_codon in STOP_CODONS
        and completeness in complete_labels
        and declared_length >= 6
    ):
        declared_length -= 3
    return declared_length - declared_length % 3


def _compile_blocks(
    starts: Sequence[int],
    ends: Sequence[int],
    strand: str,
) -> tuple[GeometryBlock, ...]:
    """Compile blocks in transcript 5-prime-to-3-prime order."""
    genomic_blocks = list(zip(starts, ends))
    ordered = genomic_blocks if strand == "+" else list(reversed(genomic_blocks))
    compiled: list[GeometryBlock] = []
    offset = 0
    for start, end in ordered:
        compiled.append(
            GeometryBlock(start=int(start), end=int(end), offset=offset)
        )
        offset += int(end) - int(start)
    return tuple(compiled)


def _plus_downstream_intervals(
    transcript: GenePredRecord,
    boundary: int,
    nt_window: int,
) -> list[tuple[int, int]]:
    """Return plus-strand transcript intervals after an ORF boundary."""
    intervals: list[tuple[int, int]] = []
    remaining = int(nt_window)
    started = False
    for exon_start, exon_end in zip(
        transcript.exon_starts,
        transcript.exon_ends,
    ):
        if not started:
            if exon_start <= boundary <= exon_end:
                started = True
                start = max(boundary, exon_start)
            elif boundary < exon_start:
                started = True
                start = exon_start
            else:
                continue
        else:
            start = exon_start
        if start >= exon_end:
            continue
        end = min(exon_end, start + remaining)
        intervals.append((start, end))
        remaining -= end - start
        if remaining <= 0:
            break
    return intervals


def _minus_downstream_intervals(
    transcript: GenePredRecord,
    boundary: int,
    nt_window: int,
) -> list[tuple[int, int]]:
    """Return minus-strand transcript intervals after an ORF boundary."""
    intervals: list[tuple[int, int]] = []
    remaining = int(nt_window)
    started = False
    blocks = list(zip(transcript.exon_starts, transcript.exon_ends))
    for exon_start, exon_end in reversed(blocks):
        if not started:
            if exon_start <= boundary <= exon_end:
                started = True
                end = min(boundary, exon_end)
            elif boundary > exon_end:
                started = True
                end = exon_end
            else:
                continue
        else:
            end = exon_end
        if end <= exon_start:
            continue
        start = max(exon_start, end - remaining)
        intervals.append((start, end))
        remaining -= end - start
        if remaining <= 0:
            break
    return intervals


def _compile_downstream_blocks(
    row: Mapping[str, object],
    starts: Sequence[int],
    ends: Sequence[int],
    transcript: GenePredRecord | None,
    strand: str,
    nt_window: int,
) -> tuple[tuple[GeometryBlock, ...], str]:
    """Compile transcript-aware downstream geometry once."""
    if nt_window <= 0:
        return (), "disabled"
    if transcript is None:
        return (), "unavailable"
    if transcript.strand != strand:
        raise ValueError(
            f"Transcript/ORF strand mismatch for {row['orf_id']}."
        )

    boundary = max(ends) if strand == "+" else min(starts)
    intervals = (
        _plus_downstream_intervals(transcript, boundary, nt_window)
        if strand == "+"
        else _minus_downstream_intervals(transcript, boundary, nt_window)
    )
    if not intervals:
        return (), "transcript"

    genomic_intervals = sorted(intervals, key=lambda item: item[0])
    block_starts = [start for start, _ in genomic_intervals]
    block_ends = [end for _, end in genomic_intervals]
    return _compile_blocks(block_starts, block_ends, strand), "transcript"


def compile_orf_geometry(
    orf_table: object,
    genepred_records: dict[str, GenePredRecord],
    post_stop_codons: int,
) -> GeometryIndex:
    """Compile all ORF geometry once before sample processing."""
    columns = tuple(str(column) for column in orf_table.columns)
    row_iterator = orf_table.itertuples(index=False, name=None)
    geometries: list[ORFGeometry] = []
    all_groups: dict[str, list[int]] = {}
    plus_groups: dict[str, list[int]] = {}
    minus_groups: dict[str, list[int]] = {}

    for input_index, values in enumerate(row_iterator):
        row = dict(zip(columns, values))
        strand = str(row["strand"])
        if strand not in {"+", "-"}:
            raise ValueError(
                f"ORF {row['orf_id']} has invalid strand: {strand}"
            )

        starts, ends, block_source = _parse_blocks(
            row,
            genepred_records,
        )
        actual_length = sum(
            int(end) - int(start)
            for start, end in zip(starts, ends)
        )
        declared_length = int(row["nt_length"])
        if actual_length != declared_length:
            raise ValueError(
                f"ORF {row['orf_id']} exon length {actual_length} does not "
                f"match nt_length {declared_length}."
            )

        coding_length = _coding_nt_length(row, actual_length)
        blocks = _compile_blocks(starts, ends, strand)
        transcript = genepred_records.get(str(row["transcript_id"]))
        downstream_blocks, release_context = _compile_downstream_blocks(
            row=row,
            starts=starts,
            ends=ends,
            transcript=transcript,
            strand=strand,
            nt_window=int(post_stop_codons) * 3,
        )
        geometry = ORFGeometry(
            input_index=input_index,
            orf_id=str(row["orf_id"]),
            family_id=str(row.get("family_id", row["orf_id"])),
            gene_id=str(row["gene_id"]),
            transcript_id=str(row["transcript_id"]),
            chrom=str(row["chrom"]),
            strand=strand,
            category=str(row.get("category", ".")),
            genomic_start=int(row["genomic_start"]),
            genomic_end=int(row["genomic_end"]),
            nt_length=declared_length,
            coding_nt_length=coding_length,
            coding_codon_count=coding_length // 3,
            block_source=block_source,
            blocks=blocks,
            downstream_blocks=downstream_blocks,
            release_context=release_context,
        )
        geometry_index = len(geometries)
        geometries.append(geometry)
        all_groups.setdefault(geometry.chrom, []).append(geometry_index)
        if strand == "+":
            plus_groups.setdefault(geometry.chrom, []).append(geometry_index)
        else:
            minus_groups.setdefault(geometry.chrom, []).append(geometry_index)

    return GeometryIndex(
        geometries=tuple(geometries),
        all_by_chrom={
            chrom: tuple(indices) for chrom, indices in all_groups.items()
        },
        plus_by_chrom={
            chrom: tuple(indices) for chrom, indices in plus_groups.items()
        },
        minus_by_chrom={
            chrom: tuple(indices) for chrom, indices in minus_groups.items()
        },
    )


def _prefix_frame_count(length: int, frame: int) -> int:
    """Count positions with one modulo-three frame in [0, length)."""
    quotient, remainder = divmod(max(0, int(length)), 3)
    return quotient + int(frame < remainder)


def _range_frame_count(start: int, end: int, frame: int) -> int:
    """Count positions with one modulo-three frame in [start, end)."""
    return (
        _prefix_frame_count(end, frame)
        - _prefix_frame_count(start, frame)
    )


def _collect_sparse_segments(
    density: ChromDensity,
    geometry: ORFGeometry,
    *,
    downstream: bool = False,
) -> list[tuple[int, int, float]]:
    """Collect transcript-offset sparse segments for one ORF."""
    blocks = geometry.downstream_blocks if downstream else geometry.blocks
    usable_length = (
        sum(block.end - block.start for block in blocks)
        if downstream
        else geometry.coding_nt_length
    )
    segments: list[tuple[int, int, float]] = []
    if usable_length <= 0 or density.starts.size == 0:
        return segments

    for block in blocks:
        block_length = block.end - block.start
        usable_in_block = min(
            block_length,
            usable_length - block.offset,
        )
        if usable_in_block <= 0:
            break

        left = int(
            np.searchsorted(density.ends, block.start, side="right")
        )
        right = int(
            np.searchsorted(density.starts, block.end, side="left")
        )
        for interval_index in range(left, right):
            value = abs(float(density.values[interval_index]))
            if value <= 0 or not math.isfinite(value):
                continue
            overlap_start = max(
                block.start,
                int(density.starts[interval_index]),
            )
            overlap_end = min(
                block.end,
                int(density.ends[interval_index]),
            )
            if overlap_end <= overlap_start:
                continue

            if geometry.strand == "+":
                tx_start = block.offset + overlap_start - block.start
                tx_end = block.offset + overlap_end - block.start
            else:
                tx_start = block.offset + block.end - overlap_end
                tx_end = block.offset + block.end - overlap_start

            tx_start = max(0, tx_start)
            tx_end = min(usable_length, tx_end)
            if tx_end > tx_start:
                segments.append((tx_start, tx_end, value))

    return segments


def _add_segment_to_codons(
    codon_profile: np.ndarray,
    start: int,
    end: int,
    value: float,
) -> None:
    """Add one constant nucleotide segment directly to codon sums."""
    if end <= start:
        return
    first_codon = start // 3
    last_codon = (end - 1) // 3

    if first_codon == last_codon:
        codon_profile[first_codon] += (end - start) * value
        return

    first_end = (first_codon + 1) * 3
    codon_profile[first_codon] += (first_end - start) * value

    last_start = last_codon * 3
    codon_profile[last_codon] += (end - last_start) * value

    if last_codon > first_codon + 1:
        codon_profile[first_codon + 1:last_codon] += 3.0 * value


def _accumulate_sparse_core(
    geometry: ORFGeometry,
    density: ChromDensity,
) -> SparseCoreProfile:
    """Map sparse P-sites directly to codons and reading frames."""
    codon_count = geometry.coding_codon_count
    zero_profile = np.zeros(0, dtype=np.float64)
    if codon_count == 0:
        return SparseCoreProfile(
            codon_profile=zero_profile,
            rpf_sum=0.0,
            covered_nt=0,
            covered_codon=0,
            max_density=0.0,
            frame0_density=0.0,
            frame1_density=0.0,
            frame2_density=0.0,
        )

    segments = _collect_sparse_segments(density, geometry)
    if not segments:
        # Avoid allocating a codon-length zero vector for NoEvidence ORFs.
        return SparseCoreProfile(
            codon_profile=zero_profile,
            rpf_sum=0.0,
            covered_nt=0,
            covered_codon=0,
            max_density=0.0,
            frame0_density=0.0,
            frame1_density=0.0,
            frame2_density=0.0,
        )

    codon_profile = np.zeros(codon_count, dtype=np.float64)
    rpf_sum = 0.0
    covered_nt = 0
    max_density = 0.0
    frame_density = [0.0, 0.0, 0.0]

    for start, end, value in segments:
        segment_length = end - start
        rpf_sum += segment_length * value
        covered_nt += segment_length
        max_density = max(max_density, value)
        for frame in range(3):
            frame_density[frame] += (
                _range_frame_count(start, end, frame) * value
            )
        _add_segment_to_codons(
            codon_profile,
            start,
            end,
            value,
        )

    return SparseCoreProfile(
        codon_profile=codon_profile,
        rpf_sum=float(rpf_sum),
        covered_nt=int(covered_nt),
        covered_codon=int(np.count_nonzero(codon_profile > 0)),
        max_density=float(max_density),
        frame0_density=float(frame_density[0]),
        frame1_density=float(frame_density[1]),
        frame2_density=float(frame_density[2]),
    )


def _core_metrics(
    core: SparseCoreProfile,
    geometry: ORFGeometry,
    thresholds: EvidenceThresholds,
) -> tuple[dict[str, object], dict[str, object], bool]:
    """Build quantification and periodicity metrics and evaluate core gate."""
    nt_count = geometry.coding_nt_length
    codon_count = geometry.coding_codon_count
    rpf_sum = core.rpf_sum
    covered_codon_ratio = safe_divide(core.covered_codon, codon_count)
    total = (
        core.frame0_density
        + core.frame1_density
        + core.frame2_density
    )
    frame0_ratio = safe_divide(core.frame0_density, total)
    frame1_ratio = safe_divide(core.frame1_density, total)
    frame2_ratio = safe_divide(core.frame2_density, total)
    alternate_mean = (
        core.frame1_density + core.frame2_density
    ) / 2.0
    frame0_vs_alt = safe_divide(
        core.frame0_density + thresholds.pseudocount,
        alternate_mean + thresholds.pseudocount,
    )
    periodicity_score = max(
        0.0,
        min(
            1.0,
            safe_divide(frame0_ratio - 1.0 / 3.0, 2.0 / 3.0),
        ),
    )
    periodicity_evaluable = (
        total >= thresholds.min_rpf_sum
        and core.covered_codon >= thresholds.min_covered_codon
    )
    if not periodicity_evaluable:
        periodicity_label = "NA"
    elif frame0_ratio >= thresholds.strong_periodicity:
        periodicity_label = "Strong"
    elif frame0_ratio >= thresholds.moderate_periodicity:
        periodicity_label = "Moderate"
    else:
        periodicity_label = "Weak"

    quant = {
        "rpf_sum": rpf_sum,
        "rpf_mean": safe_divide(rpf_sum, nt_count),
        "rpf_per_codon": safe_divide(rpf_sum, codon_count),
        "covered_nt": core.covered_nt,
        "coverage_ratio": safe_divide(core.covered_nt, nt_count),
        "covered_codon": core.covered_codon,
        "covered_codon_ratio": covered_codon_ratio,
        "max_density": core.max_density,
    }
    periodicity = {
        "frame0_density": core.frame0_density,
        "frame1_density": core.frame1_density,
        "frame2_density": core.frame2_density,
        "frame0_ratio": frame0_ratio,
        "frame1_ratio": frame1_ratio,
        "frame2_ratio": frame2_ratio,
        "frame0_vs_alt_ratio": frame0_vs_alt,
        "periodicity_score": periodicity_score,
        "periodicity_label": periodicity_label,
        "periodicity_evaluable": periodicity_evaluable,
    }
    core_pass = (
        rpf_sum >= thresholds.min_rpf_sum
        and core.covered_codon >= thresholds.min_covered_codon
        and covered_codon_ratio >= thresholds.min_codon_coverage
    )
    return quant, periodicity, core_pass


def _empty_pausing() -> dict[str, object]:
    """Return unevaluated pausing metrics."""
    return {
        "start_codon_mean": 0.0,
        "body_codon_mean": 0.0,
        "pre_stop_codon_mean": 0.0,
        "start_pause_ratio": 0.0,
        "stop_pause_ratio": 0.0,
        "start_pause_label": "NA",
        "stop_pause_label": "NA",
        "pausing_label": "NA",
        "pausing_evaluable": False,
    }


def _empty_release(context: str) -> dict[str, object]:
    """Return unevaluated release metrics."""
    return {
        "post_stop_mean": 0.0,
        "release_ratio": 0.0,
        "release_drop_score": 0.0,
        "release_label": "NA",
        "release_evaluable": False,
        "release_context": context,
        "post_stop_nt": 0,
    }


def _empty_shape() -> dict[str, object]:
    """Return unevaluated coverage-shape metrics."""
    return {
        "codon_gini": 0.0,
        "max_to_mean_ratio": 0.0,
        "top10_fraction": 0.0,
        "coverage_shape": "NA",
        "shape_evaluable": False,
    }


def _coverage_shape(
    codon_profile: np.ndarray,
    thresholds: EvidenceThresholds,
) -> dict[str, object]:
    """Calculate shape metrics with one sort operation."""
    values = np.asarray(codon_profile, dtype=np.float64)
    count = int(values.size)
    total = float(values.sum())
    if count < thresholds.min_shape_codons or total <= 0:
        return _empty_shape()

    covered_ratio = safe_divide(int(np.count_nonzero(values > 0)), count)
    mean_value = total / count
    max_value = float(values.max())
    max_to_mean = safe_divide(max_value, mean_value)
    sorted_values = np.sort(values)
    top_count = max(1, int(math.ceil(count * 0.10)))
    top_fraction = safe_divide(
        float(sorted_values[-top_count:].sum()),
        total,
    )

    indices = np.arange(1, count + 1, dtype=np.float64)
    gini = (
        2.0 * float(np.sum(indices * sorted_values)) / (count * total)
        - (count + 1.0) / count
    )
    gini = float(max(0.0, min(1.0, gini)))

    if (
        covered_ratio >= thresholds.uniform_coverage_ratio
        and gini <= thresholds.uniform_gini
        and max_to_mean <= thresholds.uniform_max_to_mean
    ):
        label = "Uniform"
    elif (
        max_to_mean >= thresholds.skewed_max_to_mean
        or top_fraction >= thresholds.skewed_top_fraction
    ):
        label = "Skewed"
    else:
        label = "Intermediate"

    return {
        "codon_gini": gini,
        "max_to_mean_ratio": max_to_mean,
        "top10_fraction": top_fraction,
        "coverage_shape": label,
        "shape_evaluable": True,
    }


def _release_metrics(
    geometry: ORFGeometry,
    density: ChromDensity,
    coding_codon_profile: np.ndarray,
    post_stop_codons: int,
    thresholds: EvidenceThresholds,
) -> dict[str, object]:
    """Calculate release only after the ORF passes the core gate."""
    if (
        geometry.release_context != "transcript"
        or not geometry.downstream_blocks
        or len(coding_codon_profile) < 3
        or post_stop_codons <= 0
    ):
        return _empty_release(geometry.release_context)

    downstream_nt = sum(
        block.end - block.start
        for block in geometry.downstream_blocks
    )
    usable_nt = min(downstream_nt, int(post_stop_codons) * 3)
    usable_nt -= usable_nt % 3
    if usable_nt < 3:
        return _empty_release(geometry.release_context)

    segments = _collect_sparse_segments(
        density,
        geometry,
        downstream=True,
    )
    downstream_codon = np.zeros(
        usable_nt // 3,
        dtype=np.float64,
    )
    for start, end, value in segments:
        clipped_end = min(end, usable_nt)
        if clipped_end > start:
            _add_segment_to_codons(
                downstream_codon,
                start,
                clipped_end,
                value,
            )

    pre_stop_mean = float(np.mean(coding_codon_profile[-3:]))
    post_stop_mean = float(np.mean(downstream_codon))
    pseudo = thresholds.pseudocount
    release_ratio = safe_divide(
        pre_stop_mean + pseudo,
        post_stop_mean + pseudo,
    )
    drop_score = max(
        0.0,
        min(
            1.0,
            1.0 - safe_divide(post_stop_mean, pre_stop_mean),
        ),
    )
    if pre_stop_mean <= 0:
        label = "Weak"
    elif release_ratio >= thresholds.strong_release:
        label = "Strong"
    elif release_ratio >= thresholds.moderate_release:
        label = "Moderate"
    else:
        label = "Weak"

    return {
        "post_stop_mean": post_stop_mean,
        "release_ratio": release_ratio,
        "release_drop_score": drop_score,
        "release_label": label,
        "release_evaluable": True,
        "release_context": geometry.release_context,
        "post_stop_nt": usable_nt,
    }


def _base_result(
    geometry: ORFGeometry,
    track: DensityTrack,
) -> dict[str, object]:
    """Build static output fields from compiled geometry."""
    return {
        "sample": track.sample,
        "density_strand": track.strand,
        "orf_id": geometry.orf_id,
        "family_id": geometry.family_id,
        "gene_id": geometry.gene_id,
        "transcript_id": geometry.transcript_id,
        "chrom": geometry.chrom,
        "strand": geometry.strand,
        "category": geometry.category,
        "genomic_start": geometry.genomic_start,
        "genomic_end": geometry.genomic_end,
        "nt_length": geometry.nt_length,
        "coding_nt_length": geometry.coding_nt_length,
        "coding_codon_count": geometry.coding_codon_count,
        "block_source": geometry.block_source,
        "profile_status": "OK",
        "profile_reason": "PASS",
    }


def score_one_geometry(
    geometry: ORFGeometry,
    density: ChromDensity,
    track: DensityTrack,
    thresholds: EvidenceThresholds,
    post_stop_codons: int,
) -> dict[str, object]:
    """Score one compiled ORF using sparse P-site accumulation."""
    base = _base_result(geometry, track)
    try:
        core = _accumulate_sparse_core(geometry, density)
        quant, periodicity, core_pass = _core_metrics(
            core,
            geometry,
            thresholds,
        )
        base.update(quant)
        base.update(periodicity)

        if core.rpf_sum <= 0:
            base.update(_empty_pausing())
            base.update(_empty_release(geometry.release_context))
            base.update(_empty_shape())
            base["evidence_components"] = "none"
            base["evidence_reason"] = "zero_rpf"
            base["translation_evidence"] = "NoEvidence"
            return base

        if not core_pass:
            failures: list[str] = []
            if core.rpf_sum < thresholds.min_rpf_sum:
                failures.append("low_rpf_sum")
            if core.covered_codon < thresholds.min_covered_codon:
                failures.append("few_covered_codons")
            if (
                safe_divide(
                    core.covered_codon,
                    geometry.coding_codon_count,
                )
                < thresholds.min_codon_coverage
            ):
                failures.append("low_codon_coverage")
            base.update(_empty_pausing())
            base.update(_empty_release(geometry.release_context))
            base.update(_empty_shape())
            base["evidence_components"] = "abundance_only"
            base["evidence_reason"] = ";".join(failures)
            base["translation_evidence"] = "LowConfidence"
            return base

        # Expensive supporting metrics are calculated only after core gate.
        pausing = calculate_pausing(core.codon_profile, thresholds)
        release = _release_metrics(
            geometry=geometry,
            density=density,
            coding_codon_profile=core.codon_profile,
            post_stop_codons=post_stop_codons,
            thresholds=thresholds,
        )
        shape = _coverage_shape(core.codon_profile, thresholds)
        evidence, components, reason = classify_translation_evidence(
            quant=quant,
            periodicity=periodicity,
            release=release,
            coverage_shape=shape,
            thresholds=thresholds,
        )
        base.update(pausing)
        base.update(release)
        base.update(shape)
        base["evidence_components"] = components
        base["evidence_reason"] = reason
        base["translation_evidence"] = evidence
        return base
    except Exception as error:
        base.update(
            {
                "profile_status": "InvalidProfile",
                "profile_reason": str(error),
                "rpf_sum": 0.0,
                "rpf_mean": 0.0,
                "rpf_per_codon": 0.0,
                "covered_nt": 0,
                "coverage_ratio": 0.0,
                "covered_codon": 0,
                "covered_codon_ratio": 0.0,
                "max_density": 0.0,
                "frame0_density": 0.0,
                "frame1_density": 0.0,
                "frame2_density": 0.0,
                "frame0_ratio": 0.0,
                "frame1_ratio": 0.0,
                "frame2_ratio": 0.0,
                "frame0_vs_alt_ratio": 0.0,
                "periodicity_score": 0.0,
                "periodicity_label": "NA",
                "periodicity_evaluable": False,
            }
        )
        base.update(_empty_pausing())
        base.update(_empty_release("invalid"))
        base.update(_empty_shape())
        base["evidence_components"] = "none"
        base["evidence_reason"] = "invalid_profile"
        base["translation_evidence"] = "NoEvidence"
        return base


def _format_value(value: object) -> str:
    """Format one TSV value without csv.DictWriter overhead."""
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return format(value, ".10g")
    if isinstance(value, (np.floating,)):
        numeric = float(value)
        return "NA" if not math.isfinite(numeric) else format(numeric, ".10g")
    if isinstance(value, (bool, np.bool_)):
        return "True" if bool(value) else "False"
    text = str(value)
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _partial_output_path(final_path: Path) -> Path:
    """Return an uncompressed, inspectable partial-output path."""
    name = final_path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
    if name.lower().endswith(".tsv"):
        name = name[:-4]
    elif name.lower().endswith(".txt"):
        name = name[:-4]
    return final_path.with_name(name + ".partial.tsv")


def _summary_output_path(final_path: Path) -> Path:
    """Return one sample summary path."""
    name = final_path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
    for suffix in (".tsv", ".txt"):
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    return final_path.with_name(name + ".summary.txt")


def _input_signature(path: str | Path) -> dict[str, object]:
    """Build a stable path-size-mtime signature."""
    resolved = Path(path).expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


class _CheckpointWriter:
    """Write batches and persist chromosome-level resume checkpoints."""

    def __init__(
        self,
        final_path: str | Path,
        columns: Sequence[str],
        signature: Mapping[str, object],
        restart: bool,
    ) -> None:
        self.final_path = Path(final_path)
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        self.partial_path = _partial_output_path(self.final_path)
        self.state_path = Path(str(self.final_path) + ".checkpoint.json")
        self.columns = tuple(columns)
        self.signature = dict(signature)
        self.restart = bool(restart)
        self.handle = None
        self.completed_units: set[str] = set()
        self.progress_unit: str | None = None
        self.progress_index = 0
        self.counters = EvidenceCounters()
        self.resumed = False

    def __enter__(self) -> "_CheckpointWriter":
        """Open a new or resumable partial output."""
        if self.restart:
            self.partial_path.unlink(missing_ok=True)
            self.state_path.unlink(missing_ok=True)

        if self.partial_path.exists() != self.state_path.exists():
            warning_print(
                "Incomplete evidence checkpoint pair detected and restarted."
            )
            self.partial_path.unlink(missing_ok=True)
            self.state_path.unlink(missing_ok=True)

        if self.partial_path.exists() and self.state_path.exists():
            state = json.loads(self.state_path.read_text(encoding="utf-8"))
            if (
                int(state.get("version", -1)) == CHECKPOINT_VERSION
                and state.get("signature") == self.signature
                and tuple(state.get("columns", [])) == self.columns
            ):
                offset = int(state["offset"])
                with self.partial_path.open("r+b") as handle:
                    handle.truncate(offset)
                self.completed_units = set(
                    str(value)
                    for value in state.get("completed_units", [])
                )
                progress_unit = state.get("progress_unit")
                self.progress_unit = (
                    None
                    if progress_unit in {None, ""}
                    else str(progress_unit)
                )
                self.progress_index = int(
                    state.get("progress_index", 0)
                )
                self.counters = EvidenceCounters.from_mapping(
                    state.get("counters", {})
                )
                self.resumed = True
            else:
                warning_print(
                    "Existing evidence checkpoint does not match the current "
                    "inputs or parameters and will be restarted."
                )
                self.partial_path.unlink(missing_ok=True)
                self.state_path.unlink(missing_ok=True)

        if not self.partial_path.exists():
            with self.partial_path.open("wb") as handle:
                header = "\t".join(self.columns) + "\n"
                handle.write(header.encode("utf-8"))

        self.handle = self.partial_path.open("ab", buffering=1024 * 1024)
        return self

    def write_rows(self, rows: Sequence[Mapping[str, object]]) -> None:
        """Write one batch of compact or full TSV rows."""
        if not rows:
            return
        if self.handle is None:
            raise RuntimeError("Checkpoint writer is not active.")
        text = "\n".join(
            "\t".join(
                _format_value(row.get(column, ""))
                for column in self.columns
            )
            for row in rows
        ) + "\n"
        self.handle.write(text.encode("utf-8"))
        self.counters.written += len(rows)

    def _save_checkpoint(self) -> None:
        """Flush output and atomically persist the current resume state."""
        if self.handle is None:
            raise RuntimeError("Checkpoint writer is not active.")
        self.handle.flush()
        os.fsync(self.handle.fileno())
        offset = int(self.handle.tell())
        state = {
            "version": CHECKPOINT_VERSION,
            "signature": self.signature,
            "columns": list(self.columns),
            "offset": offset,
            "completed_units": sorted(self.completed_units),
            "progress_unit": self.progress_unit,
            "progress_index": self.progress_index,
            "counters": asdict(self.counters),
        }
        temporary = self.state_path.with_name(
            self.state_path.name + ".tmp"
        )
        temporary.write_text(
            json.dumps(state, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        os.replace(temporary, self.state_path)

    def checkpoint_progress(
        self,
        unit: str,
        next_index: int,
    ) -> None:
        """Persist an in-progress chromosome cursor."""
        self.progress_unit = str(unit)
        self.progress_index = int(next_index)
        self._save_checkpoint()

    def checkpoint(self, unit: str) -> None:
        """Mark one chromosome complete and persist the state."""
        self.completed_units.add(str(unit))
        self.progress_unit = None
        self.progress_index = 0
        self._save_checkpoint()

    def finalize(self) -> None:
        """Commit the partial file to the requested final format."""
        if self.handle is not None:
            self.handle.flush()
            os.fsync(self.handle.fileno())
            self.handle.close()
            self.handle = None

        if str(self.final_path).lower().endswith(".gz"):
            temporary = self.final_path.with_name(
                self.final_path.name + ".tmp"
            )
            with self.partial_path.open("rb") as source:
                with gzip.open(temporary, "wb", compresslevel=6) as target:
                    shutil.copyfileobj(
                        source,
                        target,
                        length=COPY_BUFFER_SIZE,
                    )
            os.replace(temporary, self.final_path)
            self.partial_path.unlink(missing_ok=True)
        else:
            os.replace(self.partial_path, self.final_path)

        self.state_path.unlink(missing_ok=True)

    def close_preserving_partial(self) -> None:
        """Close without deleting the resumable partial output."""
        if self.handle is not None:
            self.handle.flush()
            os.fsync(self.handle.fileno())
            self.handle.close()
            self.handle = None

    def __exit__(self, exception_type, exception, traceback) -> bool:
        """Preserve partial output when interrupted or failed."""
        if exception_type is not None:
            self.close_preserving_partial()
        return False


def _minimum_evidence_rank(value: str) -> int:
    """Convert a minimum evidence label to a numeric rank."""
    return {
        "all": 0,
        "low": 1,
        "medium": 2,
        "high": 3,
    }[str(value).lower()]


def _geometry_groups_for_track(
    geometry_index: GeometryIndex,
    track: DensityTrack,
) -> Mapping[str, tuple[int, ...]]:
    """Return chromosome lookup appropriate for a density strand."""
    if track.strand == "+":
        return geometry_index.plus_by_chrom
    if track.strand == "-":
        return geometry_index.minus_by_chrom
    return geometry_index.all_by_chrom


def _sample_groups(
    tracks: Sequence[DensityTrack],
) -> list[tuple[str, tuple[DensityTrack, ...]]]:
    """Group density tracks by sample while preserving input order."""
    grouped: dict[str, list[DensityTrack]] = {}
    for track in tracks:
        grouped.setdefault(track.sample, []).append(track)
    return [
        (sample, tuple(sample_tracks))
        for sample, sample_tracks in grouped.items()
    ]


def _validate_sample_output_names(
    output_template: str | Path,
    sample_groups: Sequence[tuple[str, tuple[DensityTrack, ...]]],
) -> None:
    """Ensure sample names do not produce output collisions."""
    observed: dict[str, str] = {}
    for sample, _ in sample_groups:
        output_path = sample_output_path(output_template, sample)
        if output_path in observed:
            raise ValueError(
                "Sample names produce the same output path: "
                f"{observed[output_path]} and {sample}"
            )
        observed[output_path] = sample


def _write_summary(
    path: Path,
    sample: str,
    counters: EvidenceCounters,
    sorted_track_count: int,
    resumed: bool,
    partial_path: Path,
) -> None:
    """Write one concise sample summary."""
    lines = [
        ("sample", sample),
        ("processed_orfs", counters.processed),
        ("written_rows", counters.written),
        ("NoEvidence", counters.no_evidence),
        ("LowConfidence", counters.low_confidence),
        ("MediumConfidence", counters.medium_confidence),
        ("HighConfidence", counters.high_confidence),
        ("invalid_profile", counters.invalid_profile),
        ("auto_sorted_tracks", sorted_track_count),
        ("resumed", resumed),
        ("partial_output_during_run", str(partial_path)),
    ]
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in lines:
            handle.write(f"{key}\t{value}\n")
    os.replace(temporary, path)


def _build_run_signature(
    sample: str,
    tracks: Sequence[DensityTrack],
    args_signature: Mapping[str, object],
    columns: Sequence[str],
) -> dict[str, object]:
    """Build a checkpoint signature for one sample."""
    return {
        "sample": sample,
        "tracks": [
            {
                "sample": track.sample,
                "strand": track.strand,
                "format": track.file_format,
                **_input_signature(track.path),
            }
            for track in tracks
        ],
        "args": dict(args_signature),
        "columns": list(columns),
    }


_WORKER_GEOMETRY_INDEX: GeometryIndex | None = None
_WORKER_THRESHOLDS: EvidenceThresholds | None = None
_WORKER_OUTPUT_TEMPLATE: str | None = None
_WORKER_POST_STOP_CODONS = 10
_WORKER_MINIMUM_RANK = 1
_WORKER_KEEP_NO_EVIDENCE = False
_WORKER_OUTPUT_DETAIL = "compact"
_WORKER_RESTART = False
_WORKER_ORF_SIGNATURE: dict[str, object] | None = None


def _configure_worker(
    geometry_index: GeometryIndex,
    thresholds: EvidenceThresholds,
    output_template: str,
    post_stop_codons: int,
    minimum_rank: int,
    keep_no_evidence: bool,
    output_detail: str,
    restart: bool,
    orf_signature: Mapping[str, object],
) -> None:
    """Configure process-global read-only objects shared through fork."""
    global _WORKER_GEOMETRY_INDEX
    global _WORKER_THRESHOLDS
    global _WORKER_OUTPUT_TEMPLATE
    global _WORKER_POST_STOP_CODONS
    global _WORKER_MINIMUM_RANK
    global _WORKER_KEEP_NO_EVIDENCE
    global _WORKER_OUTPUT_DETAIL
    global _WORKER_RESTART
    global _WORKER_ORF_SIGNATURE

    _WORKER_GEOMETRY_INDEX = geometry_index
    _WORKER_THRESHOLDS = thresholds
    _WORKER_OUTPUT_TEMPLATE = str(output_template)
    _WORKER_POST_STOP_CODONS = int(post_stop_codons)
    _WORKER_MINIMUM_RANK = int(minimum_rank)
    _WORKER_KEEP_NO_EVIDENCE = bool(keep_no_evidence)
    _WORKER_OUTPUT_DETAIL = str(output_detail)
    _WORKER_RESTART = bool(restart)
    _WORKER_ORF_SIGNATURE = dict(orf_signature)


def _process_sample(
    sample_task: tuple[str, tuple[DensityTrack, ...]],
) -> SampleEvidenceResult:
    """Process one sample with chromosome-level checkpointing."""
    if (
        _WORKER_GEOMETRY_INDEX is None
        or _WORKER_THRESHOLDS is None
        or _WORKER_OUTPUT_TEMPLATE is None
        or _WORKER_ORF_SIGNATURE is None
    ):
        raise RuntimeError("smORF evidence worker is not configured.")

    sample, tracks = sample_task
    output = Path(
        sample_output_path(_WORKER_OUTPUT_TEMPLATE, sample)
    )
    output_directory = output.parent
    columns = (
        FULL_OUTPUT_COLUMNS
        if _WORKER_OUTPUT_DETAIL == "full"
        else COMPACT_OUTPUT_COLUMNS
    )
    signature_args = {
        "orf_table": _WORKER_ORF_SIGNATURE,
        "thresholds": asdict(_WORKER_THRESHOLDS),
        "post_stop_codons": _WORKER_POST_STOP_CODONS,
        "minimum_rank": _WORKER_MINIMUM_RANK,
        "keep_no_evidence": _WORKER_KEEP_NO_EVIDENCE,
        "output_detail": _WORKER_OUTPUT_DETAIL,
    }
    signature = _build_run_signature(
        sample=sample,
        tracks=tracks,
        args_signature=signature_args,
        columns=columns,
    )

    sorted_track_count = 0
    message_print(
        f"Start sample={sample}, tracks={len(tracks)}."
    )

    with _CheckpointWriter(
        final_path=output,
        columns=columns,
        signature=signature,
        restart=_WORKER_RESTART,
    ) as writer:
        if writer.resumed:
            message_print(
                "Resume sample={sample}: completed_units={units:,}, "
                "processed={processed:,}, written={written:,}.".format(
                    sample=sample,
                    units=len(writer.completed_units),
                    processed=writer.counters.processed,
                    written=writer.counters.written,
                )
            )

        for track_index, track in enumerate(tracks):
            prepared = prepare_density(
                path=track.path,
                file_format=track.file_format,
                work_directory=output_directory,
            )
            if prepared.was_sorted:
                sorted_track_count += 1
                warning_print(
                    "Density file required coordinate sorting: "
                    f"sample={sample}, file={track.path}"
                )
            prepared_track = replace(
                track,
                path=prepared.path,
                file_format=(
                    "bedgraph"
                    if prepared.was_sorted
                    else track.file_format
                ),
            )
            groups = _geometry_groups_for_track(
                _WORKER_GEOMETRY_INDEX,
                prepared_track,
            )
            seen_chromosomes: set[str] = set()
            try:
                for chrom_density in iter_density_chromosomes(
                    prepared_track.path,
                    prepared_track.file_format,
                ):
                    chrom = chrom_density.chrom
                    seen_chromosomes.add(chrom)
                    geometry_indices = groups.get(chrom)
                    if geometry_indices is None:
                        continue
                    unit = f"{track_index}\t{chrom}"
                    if unit in writer.completed_units:
                        continue

                    start_position = (
                        writer.progress_index
                        if writer.progress_unit == unit
                        else 0
                    )
                    if start_position > len(geometry_indices):
                        raise RuntimeError(
                            "Checkpoint ORF cursor exceeds chromosome size."
                        )

                    rows: list[dict[str, object]] = []
                    for position in range(
                        start_position,
                        len(geometry_indices),
                    ):
                        geometry_index = geometry_indices[position]
                        geometry = (
                            _WORKER_GEOMETRY_INDEX
                            .geometries[geometry_index]
                        )
                        result = score_one_geometry(
                            geometry=geometry,
                            density=chrom_density,
                            track=prepared_track,
                            thresholds=_WORKER_THRESHOLDS,
                            post_stop_codons=_WORKER_POST_STOP_CODONS,
                        )
                        label = str(result["translation_evidence"])
                        writer.counters.add_label(label)

                        include = (
                            (
                                label == "NoEvidence"
                                and _WORKER_KEEP_NO_EVIDENCE
                            )
                            or (
                                label != "NoEvidence"
                                and EVIDENCE_LEVELS[label]
                                >= _WORKER_MINIMUM_RANK
                            )
                        )
                        if include:
                            rows.append(result)
                            if len(rows) >= WRITE_BATCH_SIZE:
                                writer.write_rows(rows)
                                rows.clear()

                        next_position = position + 1
                        if (
                            next_position < len(geometry_indices)
                            and next_position
                            % CHECKPOINT_ORF_INTERVAL
                            == 0
                        ):
                            writer.write_rows(rows)
                            rows.clear()
                            writer.checkpoint_progress(
                                unit,
                                next_position,
                            )
                            message_print(
                                "Checkpoint sample={sample}, "
                                "strand={strand}, chrom={chrom}, "
                                "ORFs={position:,}/{total:,}, "
                                "written={written:,}.".format(
                                    sample=sample,
                                    strand=prepared_track.strand,
                                    chrom=chrom,
                                    position=next_position,
                                    total=len(geometry_indices),
                                    written=writer.counters.written,
                                )
                            )

                    writer.write_rows(rows)
                    writer.checkpoint(unit)
                    message_print(
                        "Checkpoint sample={sample}, strand={strand}, "
                        "chrom={chrom}, processed={processed:,}, "
                        "written={written:,}.".format(
                            sample=sample,
                            strand=prepared_track.strand,
                            chrom=chrom,
                            processed=writer.counters.processed,
                            written=writer.counters.written,
                        )
                    )

                # Chromosomes absent from the density track are NoEvidence.
                for chrom, geometry_indices in groups.items():
                    if chrom in seen_chromosomes:
                        continue
                    unit = f"{track_index}\t{chrom}"
                    if unit in writer.completed_units:
                        continue
                    missing_count = len(geometry_indices)
                    writer.counters.processed += missing_count
                    writer.counters.no_evidence += missing_count
                    if _WORKER_KEEP_NO_EVIDENCE:
                        rows = []
                        for geometry_index in geometry_indices:
                            geometry = (
                                _WORKER_GEOMETRY_INDEX
                                .geometries[geometry_index]
                            )
                            result = _base_result(
                                geometry,
                                prepared_track,
                            )
                            result.update(
                                {
                                    "profile_status": "NoDensityChrom",
                                    "profile_reason": (
                                        "chromosome_absent_from_density"
                                    ),
                                    "rpf_sum": 0.0,
                                    "rpf_mean": 0.0,
                                    "rpf_per_codon": 0.0,
                                    "covered_nt": 0,
                                    "coverage_ratio": 0.0,
                                    "covered_codon": 0,
                                    "covered_codon_ratio": 0.0,
                                    "max_density": 0.0,
                                    "frame0_density": 0.0,
                                    "frame1_density": 0.0,
                                    "frame2_density": 0.0,
                                    "frame0_ratio": 0.0,
                                    "frame1_ratio": 0.0,
                                    "frame2_ratio": 0.0,
                                    "frame0_vs_alt_ratio": 0.0,
                                    "periodicity_score": 0.0,
                                    "periodicity_label": "NA",
                                    "periodicity_evaluable": False,
                                    "evidence_components": "none",
                                    "evidence_reason": "zero_rpf",
                                    "translation_evidence": "NoEvidence",
                                }
                            )
                            result.update(_empty_pausing())
                            result.update(_empty_release("no_density"))
                            result.update(_empty_shape())
                            rows.append(result)
                            if len(rows) >= WRITE_BATCH_SIZE:
                                writer.write_rows(rows)
                                rows.clear()
                        writer.write_rows(rows)
                    writer.checkpoint(unit)
            finally:
                prepared.cleanup()

        writer.finalize()
        summary_path = _summary_output_path(output)
        _write_summary(
            path=summary_path,
            sample=sample,
            counters=writer.counters,
            sorted_track_count=sorted_track_count,
            resumed=writer.resumed,
            partial_path=writer.partial_path,
        )

        result = SampleEvidenceResult(
            sample=sample,
            output=str(output),
            written_rows=writer.counters.written,
            total_orfs=writer.counters.processed,
            no_evidence=writer.counters.no_evidence,
            low_confidence=writer.counters.low_confidence,
            medium_confidence=writer.counters.medium_confidence,
            high_confidence=writer.counters.high_confidence,
            sorted_track_count=sorted_track_count,
            resumed=writer.resumed,
            partial_output=str(writer.partial_path),
        )

    message_print(
        "Completed sample={sample}, processed={processed:,}, "
        "written={written:,}, NoEvidence={no_evidence:,}.".format(
            sample=sample,
            processed=result.total_orfs,
            written=result.written_rows,
            no_evidence=result.no_evidence,
        )
    )
    return result


def _multiprocessing_context(
    requested_threads: int,
    sample_count: int,
) -> tuple[multiprocessing.context.BaseContext | None, int]:
    """Resolve a safe multiprocessing context."""
    workers = max(1, min(int(requested_threads), int(sample_count)))
    if workers <= 1:
        return None, 1

    methods = multiprocessing.get_all_start_methods()
    if sys.platform.startswith("linux") and "fork" in methods:
        return multiprocessing.get_context("fork"), workers

    warning_print(
        "Multiprocessing requires fork for shared ORF geometry. "
        "Falling back to one process."
    )
    return None, 1


def run_riboseq_evidence(args: object) -> EvidenceRunResult:
    """Run sparse, resumable sample-separated evidence analysis."""
    message_print("Reading smorf_cluster family table.")
    orf_table = read_orf_table(
        args.orf_table,
        args.coord_mode,
    )
    message_print(f"Loaded clustered ORFs: {len(orf_table):,}")

    tracks = read_density_list(args)
    sample_groups = _sample_groups(tracks)
    _validate_sample_output_names(args.output, sample_groups)

    genepred_records: dict[str, GenePredRecord] = {}
    genepred_path = getattr(args, "genepred", None)
    if genepred_path is not None:
        required_names = set(
            orf_table["transcript_id"].astype(str)
        )
        if not {
            "exon_starts",
            "exon_ends",
        }.issubset(orf_table.columns):
            required_names.update(orf_table["orf_id"].astype(str))
        message_print(
            "Reading selected genePred records: "
            f"required={len(required_names):,}."
        )
        genepred_records = read_genepred(
            genepred_path,
            args.coord_mode,
            required_names=required_names,
        )
        message_print(
            f"Retained genePred records: {len(genepred_records):,}"
        )
        if not genepred_records:
            warning_print(
                "No genePred records matched transcript_id. "
                "Transcript-aware release analysis is disabled."
            )

    message_print("Compiling ORF geometry once.")
    geometry_index = compile_orf_geometry(
        orf_table=orf_table,
        genepred_records=genepred_records,
        post_stop_codons=args.post_stop_codons,
    )
    message_print(
        "Compiled ORF geometry: ORFs={orfs:,}, chromosomes={chroms:,}.".format(
            orfs=len(geometry_index.geometries),
            chroms=len(geometry_index.all_by_chrom),
        )
    )
    del orf_table
    del genepred_records

    thresholds = build_thresholds(args)
    keep_no_evidence = bool(
        getattr(args, "keep_no_evidence", False)
    )
    minimum_rank = _minimum_evidence_rank(
        getattr(args, "min_evidence", "low")
    )
    output_detail = str(
        getattr(args, "output_detail", "compact")
    ).lower()
    if output_detail not in {"compact", "full"}:
        raise ValueError(
            "--output-detail must be compact or full."
        )

    requested_threads = max(
        1,
        int(getattr(args, "threads", 1)),
    )
    context, workers = _multiprocessing_context(
        requested_threads=requested_threads,
        sample_count=len(sample_groups),
    )
    orf_signature = _input_signature(args.orf_table)
    _configure_worker(
        geometry_index=geometry_index,
        thresholds=thresholds,
        output_template=str(args.output),
        post_stop_codons=args.post_stop_codons,
        minimum_rank=minimum_rank,
        keep_no_evidence=keep_no_evidence,
        output_detail=output_detail,
        restart=bool(getattr(args, "restart", False)),
        orf_signature=orf_signature,
    )

    if workers == 1:
        results = tuple(
            _process_sample(task)
            for task in sample_groups
        )
    else:
        message_print(
            f"Parallel samples: workers={workers}, "
            f"samples={len(sample_groups)}."
        )
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
        ) as executor:
            results = tuple(
                executor.map(
                    _process_sample,
                    sample_groups,
                    chunksize=1,
                )
            )

    return EvidenceRunResult(
        samples=results,
        threads=workers,
    )
