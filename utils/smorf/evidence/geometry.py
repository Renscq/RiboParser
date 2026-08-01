#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Compare spliced ORF geometries and coding frames.
# Input: ORF exon blocks and strands.
# Output: Containment, overlap, and phase relationships.

"""Compare spliced ORF geometries and coding frames."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from .input import (
    _oriented_blocks,
    _split_ints,
)


def _geometry_blocks(row: Mapping[str, Any]) -> tuple[tuple[int, int], ...]:
    """Return sorted genomic blocks from one geometry row."""
    starts = _split_ints(row["exon_starts"])
    ends = _split_ints(row["exon_ends"])
    return tuple(sorted(zip(starts, ends), key=lambda item: item[0]))


def _blocks_contain(
    parent: Sequence[tuple[int, int]],
    child: Sequence[tuple[int, int]],
) -> bool:
    """Return whether every child coding block is covered by parent blocks."""
    parent_index = 0
    for child_start, child_end in child:
        while parent_index < len(parent) and parent[parent_index][1] <= child_start:
            parent_index += 1
        if parent_index >= len(parent):
            return False
        parent_start, parent_end = parent[parent_index]
        if parent_start > child_start or parent_end < child_end:
            return False
    return True


def _block_length(blocks: Sequence[tuple[int, int]]) -> int:
    """Return the total spliced length of genomic blocks."""
    return sum(max(0, end - start) for start, end in blocks)


def _block_overlap_length(
    first: Sequence[tuple[int, int]],
    second: Sequence[tuple[int, int]],
) -> int:
    """Return the genomic overlap length between two sorted block sets."""
    first_index = 0
    second_index = 0
    overlap = 0
    while first_index < len(first) and second_index < len(second):
        first_start, first_end = first[first_index]
        second_start, second_end = second[second_index]
        overlap += max(
            0,
            min(first_end, second_end) - max(first_start, second_start),
        )
        if first_end <= second_end:
            first_index += 1
        else:
            second_index += 1
    return overlap


def _shorter_overlap_fraction(
    first: Sequence[tuple[int, int]],
    second: Sequence[tuple[int, int]],
) -> float:
    """Return overlap as a fraction of the shorter spliced ORF."""
    denominator = min(_block_length(first), _block_length(second))
    if denominator <= 0:
        return 0.0
    return _block_overlap_length(first, second) / denominator


def _first_overlap_position(
    first: Sequence[tuple[int, int]],
    second: Sequence[tuple[int, int]],
    strand: str,
) -> int | None:
    """Return the first shared genomic base in transcript orientation."""
    overlaps: list[tuple[int, int]] = []
    first_index = 0
    second_index = 0
    while first_index < len(first) and second_index < len(second):
        start = max(first[first_index][0], second[second_index][0])
        end = min(first[first_index][1], second[second_index][1])
        if end > start:
            overlaps.append((start, end))
        if first[first_index][1] <= second[second_index][1]:
            first_index += 1
        else:
            second_index += 1
    if not overlaps:
        return None
    return overlaps[0][0] if strand == "+" else overlaps[-1][1] - 1


def _same_overlap_frame(
    first_blocks: Sequence[tuple[int, int]],
    second_blocks: Sequence[tuple[int, int]],
    strand: str,
) -> bool:
    """Return whether two overlapping ORFs use the same coding frame."""
    position = _first_overlap_position(first_blocks, second_blocks, strand)
    if position is None:
        return False
    first_offset = _transcript_offset(first_blocks, strand, position)
    second_offset = _transcript_offset(second_blocks, strand, position)
    return (
        first_offset is not None
        and second_offset is not None
        and first_offset % 3 == second_offset % 3
    )


def _transcript_offset(
    blocks: Sequence[tuple[int, int]],
    strand: str,
    genomic_position: int,
) -> int | None:
    """Map one genomic base to a spliced offset from the parent ORF start."""
    offset = 0
    for start, end in _oriented_blocks(
        [item[0] for item in blocks],
        [item[1] for item in blocks],
        strand,
    ):
        if start <= genomic_position < end:
            return (
                offset + genomic_position - start
                if strand == "+"
                else offset + end - 1 - genomic_position
            )
        offset += end - start
    return None


def _same_coding_frame(
    parent_blocks: Sequence[tuple[int, int]],
    child_blocks: Sequence[tuple[int, int]],
    strand: str,
) -> bool:
    """Return whether the child start is frame-compatible with the parent."""
    oriented_child = _oriented_blocks(
        [item[0] for item in child_blocks],
        [item[1] for item in child_blocks],
        strand,
    )
    if not oriented_child:
        return False
    child_start = oriented_child[0][0] if strand == "+" else oriented_child[0][1] - 1
    offset = _transcript_offset(parent_blocks, strand, child_start)
    return offset is not None and offset % 3 == 0
