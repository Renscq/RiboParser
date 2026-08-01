#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Write smORF density count matrices.
# Input: Quantification records and sample count vectors.
# Output: Atomic wide count matrix.

"""Write smORF density count matrices."""

from __future__ import annotations

import math
import os
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from .models import (
    _BUFFER_BYTES,
    METADATA_COLUMNS,
    OUTPUT_SUFFIX,
    QuantORF,
)


def _format_density(value: float) -> str:
    """Format one non-negative finite density value compactly."""
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Invalid quantified density value: {value}")
    rounded = round(value)
    if math.isclose(value, rounded, rel_tol=0.0, abs_tol=1e-10):
        return str(int(rounded))
    return format(value, ".10g")


def output_path_from_prefix(prefix: str | Path) -> Path:
    """Derive the standard count-matrix path from an output prefix."""
    return Path(f"{prefix}{OUTPUT_SUFFIX}")


def _write_matrix(
    output_path: Path,
    records: Sequence[QuantORF],
    sample_order: Sequence[str],
    counts_by_sample: Mapping[str, np.ndarray],
) -> None:
    """Write a complete wide count matrix atomically."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + ".tmp")
    try:
        with temporary.open(
            "w",
            encoding="utf-8",
            newline="",
            buffering=_BUFFER_BYTES,
        ) as handle:
            handle.write("\t".join((*METADATA_COLUMNS, *sample_order)) + "\n")
            for record in records:
                row = [
                    record.orf_id,
                    record.gene_id,
                    record.chrom,
                    record.strand,
                    str(record.tx_start),
                    str(record.tx_end),
                    str(record.cds_start),
                    str(record.cds_end),
                    str(record.exon_count),
                    str(record.coding_nt_length),
                    str(record.coding_codon_count),
                ]
                row.extend(
                    _format_density(counts_by_sample[sample][record.order])
                    for sample in sample_order
                )
                handle.write("\t".join(row) + "\n")
        os.replace(temporary, output_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
