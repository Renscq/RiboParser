#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Define family-evidence constants and defaults.
# Input: None.
# Output: Validated evidence configuration constants.

"""Define family-evidence constants and defaults."""

from __future__ import annotations

from typing import Final

SCHEMA_VERSION: Final[int] = 9

STOP_CODONS: Final[frozenset[str]] = frozenset({"TAA", "TAG", "TGA"})

ANNOTATED_CATEGORIES: Final[frozenset[str]] = frozenset({"annotated_orf", "annotated_morf"})

FAMILY_BATCH_SIZE: Final[int] = 10_000

SQL_BATCH_SIZE: Final[int] = 50_000

MAX_WORKERS: Final[int] = 8

MEMORY_PER_WORKER_BYTES: Final[int] = 4 * 1024**3

LEVEL_RANK: Final[dict[str, int]] = {
    "NoEvidence": 0,
    "LowConfidence": 1,
    "MediumConfidence": 2,
    "HighConfidence": 3,
}

MANUAL_THRESHOLD_DEFAULTS: Final[dict[str, float | int]] = {
    "min_rpf_sum": 5.0,
    "min_rpf_per_codon": 0.10,
    "min_covered_codon": 3,
    "min_coverage_ratio": 0.10,
    "moderate_periodicity": 0.50,
    "strong_periodicity": 0.60,
    "min_window_rpf": 3.0,
    "min_window_covered": 3,
}

ADVANCED_DEFAULTS: Final[dict[str, float | int]] = {
    "short_max_codons": 30,
    "long_min_codons": 100,
    "window_codons": 20,
    "window_step_codons": 5,
    "min_supported_windows": 2,
    "min_window_gap_codons": 15,
    "min_signal_span": 0.35,
    "localized_span_max": 0.20,
    "localized_top_window_fraction": 0.70,
    "boundary_codons": 5,
    "extent_bins": 5,
    "min_extent_bins": 3,
    "start_resolution_codons": 5,
    "min_frame_margin": 0.10,
    "leading_window_codons": 20,
    "min_exclusion_codons": 10,
    "silent_extension_density_ratio": 0.20,
    "noncanonical_min_exclusive_codons": 8,
    "noncanonical_override_density_ratio": 0.75,
    "noncanonical_override_min_coverage_ratio": 0.50,
    "noncanonical_override_min_frame_margin": 0.25,
    "nested_min_frame_margin": 0.15,
    "nested_min_phase_rpf": 10.0,
    "high_overlap_fraction": 0.70,
    "overlap_min_frame_margin": 0.25,
    "overlap_min_phase_rpf": 15.0,
}

SHORT_MODEL_FLOOR: Final[int] = 20

SHORT_MODEL_CEILING: Final[int] = 30

LONG_MODEL_FLOOR: Final[int] = 90

LONG_MODEL_CEILING: Final[int] = 100

TARGET_LENGTH_CATEGORIES: Final[frozenset[str]] = frozenset({"uorf", "dorf", "lncorf"})
