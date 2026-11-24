#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-11-24
# Version: dev001
# Function: Compare strand-aware spliced ORF coordinates.
# Input: Half-open genomic exon blocks in transcript order.
# Output: Ordered base coordinates and splice-suffix compatibility.

"""Compare spliced ORF coordinates in transcript order."""

from collections.abc import Iterator, Sequence
from itertools import islice


def spliced_positions(
    blocks: Sequence[tuple[int, int]], reverse: bool
) -> Iterator[int]:
    """Yield zero-based genomic positions in transcript order."""
    for start, end in blocks:
        yield from range(end - 1, start - 1, -1) if reverse else range(start, end)


def terminal_positions(blocks: Sequence[tuple[int, int]], reverse: bool) -> tuple[int, ...]:
    """Return the ordered triplet at the transcript's 3-prime end."""
    return tuple(islice(spliced_positions(blocks[::-1], not reverse), 3))[::-1]


def is_spliced_suffix(
    longer: Sequence[tuple[int, int]],
    shorter: Sequence[tuple[int, int]],
    reverse: bool,
) -> bool:
    """Check that two ordered block paths differ only at the upstream extent."""
    if not shorter or len(shorter) > len(longer):
        return False
    tail = longer[-len(shorter):]
    if tail[1:] != shorter[1:]:
        return False
    long_start, long_end = tail[0]
    short_start, short_end = shorter[0]
    if reverse:
        return short_start == long_start and short_end <= long_end
    return short_end == long_end and short_start >= long_start
