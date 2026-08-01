#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Expose family-aware translation-evidence analysis.
# Input: None.
# Output: Evidence pipeline API.

"""Family-aware Ribo-seq evidence analysis for clustered smORFs."""

from .config import ADVANCED_DEFAULTS, MANUAL_THRESHOLD_DEFAULTS
from .pipeline import run_family_evidence_engine

__all__ = [
    "ADVANCED_DEFAULTS",
    "MANUAL_THRESHOLD_DEFAULTS",
    "run_family_evidence_engine",
]
