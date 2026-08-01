#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Define smORF quantification data models.
# Input: Quantification settings and ORF geometry.
# Output: Typed quantification records.

"""Define smORF quantification data models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np

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
