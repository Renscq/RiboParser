#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Define stable cross-stage smORF file contracts.
# Input: None.
# Output: Shared stage names and output suffixes.

"""Stable names used to connect the four smORF workflow stages."""

from typing import Final

STAGE_ORDER: Final[tuple[str, ...]] = (
    "scanner",
    "cluster",
    "evidence",
    "quant",
)

SCANNER_MESSAGE_SUFFIX: Final[str] = ".message.txt"
SCANNER_GENEPRED_SUFFIX: Final[str] = ".genepred"
CLUSTER_FAMILY_SUFFIX: Final[str] = ".family.message.txt"
CLUSTER_MEMBER_SUFFIX: Final[str] = ".family.members.txt"
EVIDENCE_SUFFIX: Final[str] = ".smorf_evidence.txt"
QUANT_SUFFIX: Final[str] = ".density_quant.txt"
