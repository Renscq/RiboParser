#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Expose frame-aware smORF density quantification.
# Input: None.
# Output: Quantification pipeline API.

"""Frame-aware density quantification for reliable smORFs."""

from .pipeline import run_smorf_quant

__all__ = ["run_smorf_quant"]
