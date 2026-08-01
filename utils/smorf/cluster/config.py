#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Define smORF cluster constants and output contracts.
# Input: None.
# Output: Shared filtering, family, and I/O constants.

"""Define smORF cluster constants and output contracts."""

from __future__ import annotations

from typing import Final

STOP_CODONS: Final[frozenset[bytes]] = frozenset({b"TAA", b"TAG", b"TGA"})
ANNOTATED_CATEGORIES: Final[frozenset[bytes]] = frozenset({b"annotated_ORF", b"annotated_mORF"})
ANNOTATION_BIN_SIZE: Final[int] = 32_768
ANNOTATED_OVERLAP_COLUMNS: Final[tuple[bytes, ...]] = (
    b"matched_annotated_orf_id",
    b"annotated_overlap_relation",
    b"annotated_overlap_nt",
    b"annotated_overlap_codon",
)
CLUSTER_REQUIRED_COLUMNS: Final[tuple[bytes, ...]] = (
    b"orf_id",
    b"gene_id",
    b"transcript_id",
    b"chrom",
    b"strand",
    b"source_strand",
    b"category",
    b"priority",
    b"start_codon",
    b"stop_codon",
    b"nt_length",
    b"aa_length",
    b"completeness",
    b"exon_count",
    b"exon_starts",
    b"exon_ends",
)
FILTER_COLUMNS: Final[tuple[bytes, ...]] = (
    b"kozak_pwm_name",
    b"kozak_pwm_score",
    b"kozak_pwm_level",
    b"kozak_valid_ratio",
    b"structure_status",
    b"filter_status",
    b"filter_reason",
)
CLUSTER_COLUMNS: Final[tuple[bytes, ...]] = (
    b"family_id",
    b"family_role",
    b"family_type",
    b"family_size",
    b"family_representative_count",
    b"family_transcript_count",
    b"family_exact_duplicate_count",
    b"family_alt_start_count",
    b"family_collapsed_count",
    b"family_categories",
    b"family_start_codons",
    b"primary_selection_rule",
)
MEMBER_COLUMNS: Final[tuple[bytes, ...]] = (
    b"family_id",
    b"primary_orf_id",
    b"representative_orf_id",
    b"member_orf_id",
    b"member_transcript_id",
    b"member_category",
    b"member_start_codon",
    b"member_aa_length",
    b"member_transcript_length",
    b"family_role",
    b"collapse_reason",
)
PRIMARY_SELECTION_RULE: Final[bytes] = (
    b"annotated>start_codon>aa_length>kozak_score>transcript_length>input_order"
)
_BUFFER_BYTES: Final[int] = 8 * 1024 * 1024
_FLUSH_BYTES: Final[int] = 8 * 1024 * 1024
_MIN_PARALLEL_BYTES: Final[int] = 128 * 1024 * 1024
_KOZAK_CACHE_SIZE: Final[int] = 131_072
_DEFAULT_START_INDEX: Final[int] = 6
_BASE_BYTES: Final[tuple[int, ...]] = (65, 67, 71, 84)
