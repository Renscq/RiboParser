#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Expose the transcript-centric smORF scanner pipeline.
# Input: None.
# Output: Scanner pipeline API.

"""Transcript-centric smORF candidate scanning."""

from .pipeline import SmORFPipeline

__all__ = ["SmORFPipeline"]
