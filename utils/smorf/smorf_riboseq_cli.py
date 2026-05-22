#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Rensc
date: 2026-05-23

Command-line interface for smORF Ribo-seq evidence analysis.
"""

import argparse

from .smorf_riboseq_io import eprint, write_output
from .smorf_riboseq_pipeline import run_riboseq_evidence


def get_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Evaluate smORF translation evidence using Ribo-seq P-site density."
    )

    # Required inputs
    parser.add_argument("-i", "--orf-table", required=True, help="Filtered smORF table from smorf_filter.")
    parser.add_argument("-o", "--output", required=True, help="Output evidence table in TSV format.")

    # Optional annotation
    parser.add_argument("--genepred", default=None, help="Optional genePred file for ORF exon blocks.")
    parser.add_argument("--chrom-sizes", default=None, help="Two-column chromosome size file.")

    # Density input
    parser.add_argument("--density-list", default=None, help="TSV with columns: sample, strand, path, optional format.")
    parser.add_argument("--density-plus", default=None, help="Plus-strand P-site density file.")
    parser.add_argument("--density-minus", default=None, help="Minus-strand P-site density file.")
    parser.add_argument("--density", default=None, help="Unstranded P-site density file.")
    parser.add_argument("--sample", default="sample1", help="Sample name for direct density input.")
    parser.add_argument("--density-format", default="auto", choices=["auto", "wig", "bedgraph"], help="Density file format.")

    # Coordinate settings
    parser.add_argument(
        "--coord-mode",
        default="0based-half-open",
        choices=["0based-half-open", "1based-closed"],
        help="Coordinate mode for ORF table and genePred-like blocks.",
    )

    # Windows
    parser.add_argument("--post-stop-codons", type=int, default=10, help="Number of codons after stop codon used for release signal.")
    parser.add_argument("--pseudocount", type=float, default=0.1, help="Pseudocount for ratio calculation.")

    # Basic filters
    parser.add_argument("--min-rpf-sum", type=float, default=3.0, help="Minimum ORF-level RPF sum for evidence scoring.")
    parser.add_argument("--min-covered-codon", type=int, default=2, help="Minimum covered codon count.")
    parser.add_argument("--min-coverage-ratio", type=float, default=0.10, help="Minimum nucleotide-level coverage ratio.")

    # Periodicity thresholds
    parser.add_argument("--strong-periodicity", type=float, default=0.70, help="Frame-0 ratio threshold for strong periodicity.")
    parser.add_argument("--moderate-periodicity", type=float, default=0.55, help="Frame-0 ratio threshold for moderate periodicity.")

    # Pausing thresholds
    parser.add_argument("--strong-start-pause", type=float, default=1.50, help="Start pausing ratio threshold for strong signal.")
    parser.add_argument("--moderate-start-pause", type=float, default=1.20, help="Start pausing ratio threshold for moderate signal.")
    parser.add_argument("--strong-stop-pause", type=float, default=1.50, help="Pre-stop pausing ratio threshold for strong signal.")
    parser.add_argument("--moderate-stop-pause", type=float, default=1.20, help="Pre-stop pausing ratio threshold for moderate signal.")

    # Release thresholds
    parser.add_argument("--strong-release", type=float, default=3.0, help="Release ratio threshold for strong signal.")
    parser.add_argument("--moderate-release", type=float, default=1.5, help="Release ratio threshold for moderate signal.")

    # Coverage shape thresholds
    parser.add_argument("--uniform-coverage-ratio", type=float, default=0.40, help="Coverage ratio threshold for Uniform shape.")
    parser.add_argument("--uniform-gini", type=float, default=0.50, help="Gini threshold for Uniform shape.")
    parser.add_argument("--uniform-max-to-mean", type=float, default=5.0, help="Max/mean threshold for Uniform shape.")
    parser.add_argument("--skewed-max-to-mean", type=float, default=10.0, help="Max/mean threshold for Skewed shape.")
    parser.add_argument("--skewed-top-fraction", type=float, default=0.70, help="Top 10 percent density fraction threshold for Skewed shape.")
    parser.add_argument("--disperse-coverage-ratio", type=float, default=0.20, help="Coverage ratio threshold below which coverage is Disperse.")

    # Runtime
    parser.add_argument("--progress-every", type=int, default=10000, help="Print progress every N ORFs per chromosome.")
    parser.add_argument("--keep-no-evidence", action="store_true", help="Keep ORFs with no RPF evidence in output.")

    return parser


def main() -> None:
    """Run command-line workflow."""
    parser = get_parser()
    args = parser.parse_args()

    output_table = run_riboseq_evidence(args)

    eprint(f"[Info] Writing output: {args.output}")
    write_output(output_table, args.output)
    eprint("[Info] Done.")
