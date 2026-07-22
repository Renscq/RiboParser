#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Extract ORF and transcript-downstream P-site density profiles.
# Input: Sparse chromosome density, ORF blocks, strand, and transcript annotation.
# Output: Transcript-oriented nucleotide and codon density profiles.

"""Profile extraction for smORF Ribo-seq evidence analysis."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from .smorf_riboseq_constants import STOP_CODONS
from .smorf_riboseq_density import ChromDensity
from .smorf_riboseq_io import GenePredRecord


def validate_blocks(
    starts: Sequence[int],
    ends: Sequence[int],
    *,
    orf_id: str = "unknown",
) -> tuple[list[int], list[int]]:
    """Validate and sort genomic exon blocks.

    Args:
        starts: Genomic starts.
        ends: Genomic ends.
        orf_id: Identifier used in error messages.

    Returns:
        Sorted start and end lists.

    Raises:
        ValueError: If blocks are empty, overlapping, or invalid.
    """
    if len(starts) != len(ends) or not starts:
        raise ValueError(
            f"ORF {orf_id} has mismatched or empty exon blocks."
        )

    blocks = sorted(
        (
            (int(start), int(end))
            for start, end in zip(starts, ends)
        ),
        key=lambda block: block[0],
    )

    previous_end = None
    for start, end in blocks:
        if start < 0 or end <= start:
            raise ValueError(
                f"ORF {orf_id} has invalid block [{start}, {end})."
            )
        if previous_end is not None and start < previous_end:
            raise ValueError(
                f"ORF {orf_id} has overlapping exon blocks."
            )
        previous_end = end

    return (
        [start for start, _ in blocks],
        [end for _, end in blocks],
    )


def extract_transcript_profile(
    density: ChromDensity | None,
    starts: Sequence[int],
    ends: Sequence[int],
    strand: str,
    *,
    orf_id: str = "unknown",
) -> np.ndarray:
    """Extract an ORF profile in transcript 5-prime-to-3-prime order.

    Args:
        density: Sparse chromosome density. ``None`` represents an entirely
            uncovered chromosome.
        starts: ORF genomic block starts.
        ends: ORF genomic block ends.
        strand: ORF genomic strand.
        orf_id: Identifier used in errors.

    Returns:
        Float32 nucleotide profile.

    Raises:
        ValueError: If the strand or blocks are invalid.
    """
    if strand not in {"+", "-"}:
        raise ValueError(f"Invalid ORF strand: {strand}")

    starts, ends = validate_blocks(
        starts,
        ends,
        orf_id=orf_id,
    )

    block_iterator = (
        zip(starts, ends)
        if strand == "+"
        else zip(reversed(starts), reversed(ends))
    )
    pieces: list[np.ndarray] = []

    for start, end in block_iterator:
        piece = (
            np.zeros(end - start, dtype=np.float32)
            if density is None
            else density.query(start, end)
        )
        if strand == "-":
            piece = piece[::-1]
        pieces.append(piece)

    return (
        np.concatenate(pieces).astype(np.float32, copy=False)
        if pieces
        else np.zeros(0, dtype=np.float32)
    )


def get_coding_nt_length(
    row: pd.Series,
    profile_length: int,
) -> int:
    """Return coding length excluding a complete terminal stop codon.

    Args:
        row: ORF metadata row.
        profile_length: Extracted ORF profile length.

    Returns:
        Non-negative coding length divisible by three.
    """
    declared_length = min(int(row["nt_length"]), int(profile_length))
    stop_codon = (
        str(row.get("stop_codon", ""))
        .strip()
        .upper()
        .replace("U", "T")
    )
    completeness = str(
        row.get("completeness", "")
    ).strip().lower()

    complete_labels = {"complete", "cmpl", "full", "true", "yes"}
    has_complete_stop = (
        stop_codon in STOP_CODONS
        and (
            completeness in complete_labels
            or completeness == ""
        )
    )

    if has_complete_stop and declared_length >= 6:
        declared_length -= 3

    return declared_length - (declared_length % 3)


def make_codon_profile(
    nt_profile: np.ndarray,
    coding_nt_length: int,
) -> np.ndarray:
    """Aggregate coding P-site density into codons.

    Args:
        nt_profile: Transcript-oriented ORF nucleotide profile.
        coding_nt_length: Coding nucleotide length excluding stop codon.

    Returns:
        Codon-summed float64 density.
    """
    usable_length = max(
        0,
        min(int(coding_nt_length), len(nt_profile)),
    )
    usable_length -= usable_length % 3
    if usable_length == 0:
        return np.zeros(0, dtype=np.float64)

    return (
        nt_profile[:usable_length]
        .reshape(-1, 3)
        .sum(axis=1, dtype=np.float64)
    )


def _plus_downstream_intervals(
    transcript: GenePredRecord,
    boundary: int,
    nt_window: int,
) -> list[tuple[int, int]]:
    """Return plus-strand transcript intervals after an ORF boundary."""
    intervals: list[tuple[int, int]] = []
    remaining = nt_window
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
    remaining = nt_window
    started = False

    blocks = list(
        zip(transcript.exon_starts, transcript.exon_ends)
    )

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


def extract_transcript_downstream_profile(
    density: ChromDensity | None,
    orf_starts: Sequence[int],
    orf_ends: Sequence[int],
    strand: str,
    transcript: GenePredRecord | None,
    nt_window: int,
    *,
    orf_id: str = "unknown",
) -> tuple[np.ndarray, str]:
    """Extract spliced downstream transcript density after the ORF stop.

    Args:
        density: Sparse chromosome density.
        orf_starts: ORF genomic block starts.
        orf_ends: ORF genomic block ends.
        strand: ORF strand.
        transcript: Source transcript genePred record.
        nt_window: Requested downstream nucleotides.
        orf_id: Identifier used in errors.

    Returns:
        Downstream profile and context label. The context is
        ``transcript`` when a matching transcript model is available and
        ``unavailable`` otherwise.

    Raises:
        ValueError: If transcript metadata conflicts with the ORF.
    """
    if nt_window <= 0:
        return np.zeros(0, dtype=np.float32), "disabled"
    if transcript is None:
        return np.zeros(0, dtype=np.float32), "unavailable"
    if transcript.strand != strand:
        raise ValueError(
            f"Transcript/ORF strand mismatch for {orf_id}."
        )

    starts, ends = validate_blocks(
        orf_starts,
        orf_ends,
        orf_id=orf_id,
    )
    boundary = max(ends) if strand == "+" else min(starts)

    intervals = (
        _plus_downstream_intervals(
            transcript,
            boundary,
            nt_window,
        )
        if strand == "+"
        else _minus_downstream_intervals(
            transcript,
            boundary,
            nt_window,
        )
    )

    pieces: list[np.ndarray] = []
    for start, end in intervals:
        piece = (
            np.zeros(end - start, dtype=np.float32)
            if density is None
            else density.query(start, end)
        )
        if strand == "-":
            piece = piece[::-1]
        pieces.append(piece)

    profile = (
        np.concatenate(pieces).astype(np.float32, copy=False)
        if pieces
        else np.zeros(0, dtype=np.float32)
    )
    return profile, "transcript"
