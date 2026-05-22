#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Rensc
date: 2026-05-23

Input and output functions for smORF Ribo-seq evidence analysis.
"""

import gzip
import os
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .smorf_riboseq_constants import DensityTrack, SUPPORTED_DENSITY_FORMATS, VALID_STRANDS


def eprint(message: str) -> None:
    """Print progress information to stderr."""
    print(message, file=sys.stderr, flush=True)


def smart_open(path: str, mode: str = "rt"):
    """Open plain or gzip-compressed text files."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_comma_ints(value) -> List[int]:
    """Parse comma-separated integer fields."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if text == "" or text == ".":
        return []
    text = text.rstrip(",")
    if text == "":
        return []
    return [int(x) for x in text.split(",") if x != ""]


def infer_density_format(path: str, user_format: Optional[str] = None) -> str:
    """Infer density file format from suffix or user-provided option."""
    if user_format and user_format != "auto":
        fmt = user_format.lower()
        if fmt not in SUPPORTED_DENSITY_FORMATS:
            raise ValueError(f"Unsupported density format: {fmt}")
        return fmt

    lower = path.lower()
    if lower.endswith((".bedgraph", ".bdg", ".bedgraph.gz", ".bdg.gz")):
        return "bedgraph"
    if lower.endswith((".wig", ".wiggle", ".wig.gz", ".wiggle.gz")):
        return "wig"

    return "wig"


def read_chrom_sizes(path: Optional[str]) -> Dict[str, int]:
    """Read chromosome sizes from a two-column file."""
    if path is None:
        return {}

    chrom_sizes = {}
    with smart_open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 2:
                continue
            chrom_sizes[fields[0]] = int(fields[1])

    return chrom_sizes


def read_density_list(args) -> List[DensityTrack]:
    """Read density tracks from command-line arguments or a list file."""
    tracks = []

    if args.density_list:
        table = pd.read_csv(args.density_list, sep="\t", comment="#")
        required = {"sample", "strand", "path"}
        missing = required - set(table.columns)
        if missing:
            raise ValueError(f"--density-list requires columns: sample, strand, path. Missing: {missing}")

        for _, row in table.iterrows():
            fmt = "auto"
            if "format" in table.columns and not pd.isna(row["format"]):
                fmt = str(row["format"])

            tracks.append(
                DensityTrack(
                    sample=str(row["sample"]),
                    strand=str(row["strand"]),
                    path=str(row["path"]),
                    file_format=fmt,
                )
            )

    if args.density_plus:
        tracks.append(DensityTrack(args.sample, "+", args.density_plus, args.density_format))

    if args.density_minus:
        tracks.append(DensityTrack(args.sample, "-", args.density_minus, args.density_format))

    if args.density:
        tracks.append(DensityTrack(args.sample, ".", args.density, args.density_format))

    if not tracks:
        raise ValueError("No density file was provided.")

    for track in tracks:
        if track.strand not in VALID_STRANDS:
            raise ValueError(f"Invalid strand in density track: {track.strand}")
        if not os.path.exists(track.path):
            raise FileNotFoundError(track.path)

    return tracks


def read_orf_table(path: str, coord_mode: str) -> pd.DataFrame:
    """Read filtered smORF table."""
    try:
        table = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        table = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")

    required = {
        "orf_id", "gene_id", "transcript_id", "chrom", "strand",
        "genomic_start", "genomic_end", "nt_length"
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"ORF table is missing required columns: {missing}")

    for column in ["genomic_start", "genomic_end", "nt_length"]:
        table[column] = pd.to_numeric(table[column], errors="coerce").astype("Int64")

    if "aa_length" in table.columns:
        table["aa_length"] = pd.to_numeric(table["aa_length"], errors="coerce").astype("Int64")

    if coord_mode == "1based-closed":
        table["genomic_start"] = table["genomic_start"] - 1
        table["genomic_end"] = table["genomic_end"]

    table = table.dropna(subset=["genomic_start", "genomic_end", "nt_length"]).copy()
    table["genomic_start"] = table["genomic_start"].astype(int)
    table["genomic_end"] = table["genomic_end"].astype(int)
    table["nt_length"] = table["nt_length"].astype(int)

    table = table[table["genomic_end"] > table["genomic_start"]].copy()

    if "filter_status" in table.columns:
        table = table[table["filter_status"].astype(str).str.upper().eq("PASS")].copy()

    table["strand"] = table["strand"].astype(str)
    table = table[table["strand"].isin(["+", "-"])].copy()

    return table.reset_index(drop=True)


def read_genepred(path: Optional[str], coord_mode: str) -> Dict[str, Tuple[List[int], List[int]]]:
    """Read genePred blocks keyed by ORF ID."""
    if path is None:
        return {}

    columns = [
        "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
        "exonCount", "exonStarts", "exonEnds", "score", "name2",
        "cdsStartStat", "cdsEndStat", "exonFrames"
    ]
    table = pd.read_csv(path, sep="\t", header=None, names=columns, dtype=str)
    block_map = {}

    for _, row in table.iterrows():
        starts = parse_comma_ints(row["exonStarts"])
        ends = parse_comma_ints(row["exonEnds"])

        if coord_mode == "1based-closed":
            starts = [x - 1 for x in starts]

        if len(starts) != len(ends) or len(starts) == 0:
            continue

        block_map[str(row["name"])] = (starts, ends)

    return block_map


def get_orf_blocks(row: pd.Series, genepred_blocks: Dict[str, Tuple[List[int], List[int]]]) -> Tuple[List[int], List[int]]:
    """Get ORF exon blocks from ORF table, genePred, or genomic interval."""
    orf_id = str(row["orf_id"])

    if "exon_starts" in row.index and "exon_ends" in row.index:
        starts = parse_comma_ints(row["exon_starts"])
        ends = parse_comma_ints(row["exon_ends"])
        if len(starts) == len(ends) and len(starts) > 0:
            return starts, ends

    if orf_id in genepred_blocks:
        return genepred_blocks[orf_id]

    return [int(row["genomic_start"])], [int(row["genomic_end"])]


def write_output(table: pd.DataFrame, output: str) -> None:
    """Write evidence table."""
    if output.endswith(".gz"):
        table.to_csv(output, sep="\t", index=False, compression="gzip")
    else:
        table.to_csv(output, sep="\t", index=False)
