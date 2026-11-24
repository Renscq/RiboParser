#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-11-24
# Version: dev001
# Function: Measure P-site support on individual smORF coding paths.
# Input: Indexed density tracks, candidate coding blocks, and splice diagnostics.
# Output: Per-sample profiles, priority tiers, and splice-path support.

"""Measure candidate-specific P-site support from the existing density cache."""

from __future__ import annotations

import sqlite3
from collections import defaultdict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from itertools import islice
from pathlib import Path

import numpy as np

from ..coordinates import spliced_positions
from ..density import ChromDensity
from .features import _track_for_family
from .geometry import _geometry_blocks
from .input import _coding_nt_length, _load_cached_density, _split_ints
from .models import DensityTrack, EngineConfig, EvidenceEngineError


MIN_PRIORITY_RPF = 10.0
MIN_PRIORITY_FRAME0_RATIO = 0.60
CANDIDATE_COLUMNS = (
    "orf_id", "family_id", "gene_id", "transcript_id", "category", "chrom", "strand",
    "start0", "aa_length", "exon_starts", "exon_ends", "sample", "group",
    "total_rpf", "frame0_rpf", "frame0_ratio", "frame0_codons", "frame0_segments",
    "priority_sample_pass", "priority_supporting_samples", "candidate_tier",
    "splice_path_role", "splice_scaffold_ids", "splice_diagnostic_codons",
    "splice_diagnostic_frame0_rpf", "splice_diagnostic_frame0_codons", "splice_sample_pass",
    "splice_supporting_samples",
)


@dataclass(frozen=True, slots=True)
class CountSummary:
    """Store depth, periodicity, and distribution for stop-excluded codons."""

    total: float
    frame0: float
    frame0_codons: int
    frame0_segments: int

    @property
    def frame0_ratio(self) -> float:
        """Return the fraction of P-site signal in the candidate reading frame."""
        return self.frame0 / self.total if self.total else 0.0


def coding_positions(
    blocks: Sequence[tuple[int, int]], strand: str, coding_nt_length: int
) -> np.ndarray:
    """Return genomic bases of complete coding codons in translation order."""
    reverse = strand == "-"
    ordered = blocks[::-1] if reverse else blocks
    positions = islice(spliced_positions(ordered, reverse), coding_nt_length)
    return np.fromiter(positions, dtype=np.int64, count=coding_nt_length)


def summarize_counts(counts: np.ndarray) -> CountSummary:
    """Describe a codon-by-frame matrix using three consecutive coding segments."""
    frame0 = counts[:, 0]
    return CountSummary(
        total=float(counts.sum()),
        frame0=float(frame0.sum()),
        frame0_codons=int(np.count_nonzero(frame0)),
        frame0_segments=sum(bool(segment.any()) for segment in np.array_split(frame0, 3)),
    )


def priority_sample_pass(summary: CountSummary, aa_length: int) -> bool:
    """Test the short-ORF-aware priority rule within one sample."""
    minimum_codons = 3 if 10 <= aa_length <= 19 else 4
    return (
        summary.total >= MIN_PRIORITY_RPF
        and summary.frame0_ratio >= MIN_PRIORITY_FRAME0_RATIO
        and summary.frame0_codons >= minimum_codons
        and summary.frame0_segments >= 2
    )


def splice_sample_pass(
    counts: np.ndarray,
    positions: np.ndarray,
    diagnostic_positions: Sequence[int],
    min_rpf: float,
) -> tuple[float, int, bool]:
    """Test support at diagnostic frame-zero positions within the selected ORF."""
    diagnostic = np.isin(positions[::3], diagnostic_positions)
    minimum_codons = min(2, int(diagnostic.sum()))
    frame0 = counts[diagnostic, 0]
    rpf = float(frame0.sum())
    occupied = int(np.count_nonzero(frame0))
    return rpf, occupied, bool(minimum_codons and occupied >= minimum_codons and rpf >= min_rpf)


class DensityProfiles:
    """Extract codon-by-frame P-site matrices and cache the active chromosome."""

    def __init__(self, work_directory: Path, tracks: Sequence[DensityTrack]) -> None:
        self.work_directory = work_directory
        self.tracks: dict[str, dict[str, DensityTrack]] = defaultdict(dict)
        self.groups: dict[str, str] = {}
        for track in tracks:
            self.tracks[track.sample][track.strand] = track
            self.groups[track.sample] = track.group
        self.chromosome = ""
        self.densities: dict[DensityTrack, ChromDensity | None] = {}

    def counts(
        self, chrom: str, strand: str, positions: np.ndarray
    ) -> Iterator[tuple[str, np.ndarray | None]]:
        """Yield one matrix per sample, with missing strand tracks marked as unavailable."""
        if self.chromosome != chrom:
            self.densities.clear()
            self.chromosome = chrom
        for sample in sorted(self.tracks):
            track = _track_for_family(self.tracks, sample, strand)
            if track is None:
                yield sample, None
                continue
            if track not in self.densities:
                self.densities[track] = _load_cached_density(self.work_directory, track, chrom)
            yield sample, self._counts_at(self.densities[track], positions)

    @staticmethod
    def _counts_at(density: ChromDensity | None, positions: np.ndarray) -> np.ndarray:
        """Read half-open density intervals at the requested genomic bases."""
        counts = np.zeros(positions.size, dtype=np.float64)
        if density is not None and density.starts.size:
            indices = np.searchsorted(density.ends, positions, side="right")
            candidates = np.flatnonzero(indices < density.starts.size)
            hits = candidates[density.starts[indices[candidates]] <= positions[candidates]]
            values = np.abs(density.values[indices[hits]])
            counts[hits] = np.where(np.isfinite(values), values, 0.0)
        return counts.reshape(-1, 3)


@dataclass(frozen=True, slots=True)
class SplicePath:
    """Describe a branch or an ORF shared by multiple coding paths."""

    role: str = ""
    scaffolds: str = ""
    diagnostic_positions: tuple[int, ...] = ()


class SplicePathEvidence:
    """Resolve path-specific evidence before dependent ORFs enter competition."""

    def __init__(
        self, connection: sqlite3.Connection, profiles: DensityProfiles, config: EngineConfig
    ) -> None:
        self.connection = connection
        self.profiles = profiles
        self.config = config
        self.paths = {
            row["family_id"]: SplicePath(
                row["splice_path_role"], row["splice_scaffold_ids"],
                _split_ints(row["splice_unique_positions"]),
            )
            for row in connection.execute(
                "SELECT family_id, splice_path_role, splice_scaffold_ids, splice_unique_positions "
                "FROM families WHERE splice_path_role != ''"
            )
        }

    def unresolved(self, family_id: str, orf_id: str) -> bool:
        """Check whether the selected ORF has enough branch-specific support."""
        path = self.paths.get(family_id)
        if path is None:
            return False
        if path.role == "SharedBody" or not path.diagnostic_positions:
            return True
        geometry = self.connection.execute(
            "SELECT * FROM geometry WHERE orf_id=?", (orf_id,)
        ).fetchone()
        if geometry is None:
            raise EvidenceEngineError("Splice-path ORF lacks geometry: {}".format(orf_id))
        return not self.supports(path, geometry)

    def supports(self, path: SplicePath, geometry: sqlite3.Row) -> bool:
        """Count diagnostic support in independent samples for one ORF."""
        chrom, strand = geometry["chrom"], geometry["strand"]
        blocks = _geometry_blocks(geometry)
        coding_nt = _coding_nt_length(dict(geometry))
        positions = coding_positions(blocks, strand, coding_nt)
        passes = sum(
            splice_sample_pass(counts, positions, path.diagnostic_positions, self.config.min_window_rpf)[2]
            for _, counts in self.profiles.counts(chrom, strand, positions)
            if counts is not None
        )
        return passes >= self.config.reliable_sample
