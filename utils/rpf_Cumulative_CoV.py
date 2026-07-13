#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.002
# Function: Calculate cumulative coefficient of variation along CDS regions.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Cumulative CoV tables, summary JSON, outlier table, and figures.

"""Command-line entry point for rpf_Cumulative_CoV."""

from __future__ import annotations

import argparse
import os

from utils.ribo import Cumulative_CoV
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Calculate cumulative coefficient of variation along CDS regions.",
    )

    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-r", "--rpf", required=True, type=str,
        help="Input RPF density file in JSONL, JSONL.GZ, TXT, or TSV format.",
    )
    required.add_argument(
        "-o", "--output", required=True, type=str,
        help="Output file prefix.",
    )

    filtering = parser.add_argument_group("Input filtering arguments")
    filtering.add_argument(
        "-l", "--list", type=str, default=None,
        help="Optional transcript ID list.",
    )
    filtering.add_argument(
        "-s", "--site", choices=["E", "P", "A"], default="P",
        help="Ribosome site used for coordinate assignment.",
    )
    filtering.add_argument(
        "-f", "--frame", choices=["0", "1", "2", "all"], default="all",
        help="Reading frame used for cumulative CoV calculation.",
    )
    filtering.add_argument(
        "-m", "--min", type=float, default=0,
        help="Minimum sample-specific CDS RPF count required per transcript.",
    )
    filtering.add_argument(
        "--tis", type=int, default=0,
        help="Number of codons discarded after TIS.",
    )
    filtering.add_argument(
        "--tts", type=int, default=0,
        help="Number of codons discarded before TTS.",
    )
    filtering.add_argument(
        "--min-positions", type=int, default=10,
        help="Minimum number of retained positions required per transcript.",
    )

    calculation = parser.add_argument_group("Cumulative CoV arguments")
    calculation.add_argument(
        "-t", "--trim", type=int, default=150,
        help="Maximum CDS-relative position included in meta analysis and plots.",
    )
    calculation.add_argument(
        "--resolution", choices=["nucleotide", "codon"], default="nucleotide",
        help="Position resolution used for cumulative CoV calculation.",
    )
    calculation.add_argument(
        "--ddof", choices=[0, 1], type=int, default=1,
        help="Delta degrees of freedom used for cumulative standard deviation.",
    )
    calculation.add_argument(
        "-n", "--normal", action="store_true",
        help="Convert density to RPM before reporting mean and SD; CoV is unchanged.",
    )
    calculation.add_argument(
        "--thread", type=int, default=1,
        help="Number of sample-level worker threads.",
    )

    outlier = parser.add_argument_group("Outlier arguments")
    outlier.add_argument(
        "--remove-outlier", action="store_true",
        help="Remove extreme local RPF pileups before cumulative CoV calculation.",
    )
    outlier.add_argument(
        "--outlier-iqr", type=float, default=8.0,
        help="IQR multiplier for the global log1p density cutoff.",
    )
    outlier.add_argument(
        "--outlier-window", type=int, default=5,
        help="Number of neighboring positions used on each side of a candidate outlier.",
    )
    outlier.add_argument(
        "--outlier-local-fold", type=float, default=10.0,
        help="Minimum fold above local background required to remove a candidate.",
    )

    plotting = parser.add_argument_group("Output and plotting arguments")
    plotting.add_argument(
        "--plot-stat", choices=["median", "mean"], default="median",
        help="Center statistic used for the cumulative CoV meta curve.",
    )
    plotting.add_argument(
        "--plot-transform", choices=["none", "sqrt", "log1p", "log2", "log10"],
        default="none", help="Transformation applied only to plotted CoV values.",
    )
    plotting.add_argument(
        "--ci-low", type=float, default=0.25,
        help="Lower transcript quantile used for the curve ribbon.",
    )
    plotting.add_argument(
        "--ci-high", type=float, default=0.75,
        help="Upper transcript quantile used for the curve ribbon.",
    )
    plotting.add_argument(
        "--gene-fig", choices=["png", "pdf", "both"], default="png",
        help="Output format for one cumulative CoV figure per transcript.",
    )
    plotting.add_argument(
        "--all", action="store_true",
        help="Output transcript-position cumulative CoV details.",
    )
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)
    if args.min < 0:
        raise ValueError("--min must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if args.trim < 2:
        raise ValueError("--trim must be >= 2.")
    if args.min_positions < 2:
        raise ValueError("--min-positions must be >= 2.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if not 0 <= args.ci_low < args.ci_high <= 1:
        raise ValueError("Require 0 <= --ci-low < --ci-high <= 1.")


def _parse_args() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args()
    _validate_args(args)
    args_print(args)
    return args


def _print_step(number: int, message: str) -> None:
    """Print a formatted pipeline step."""
    print(f"\nStep{number}: {message}", flush=True)


def _run_pipeline(args: argparse.Namespace) -> None:
    """Run cumulative CoV analysis."""
    analysis = Cumulative_CoV.CumulativeCoV(args)

    _print_step(2, "Import the RPF density file.")
    analysis.import_rpf()

    _print_step(3, "Calculate sample-specific cumulative CoV.")
    analysis.calculate_cumulative_cov()

    _print_step(4, "Draw cumulative CoV summary and per-transcript figures.")
    analysis.draw_cumulative_cov()
    analysis.draw_transcript_count()
    analysis.draw_gene_cumulative_cov()

    _print_step(5, "Output cumulative CoV results.")
    analysis.output_results()


def main() -> None:
    """Run the command-line program."""
    now_time()
    print("\nCalculate cumulative coefficient of variation.", flush=True)
    _print_step(1, "Check input arguments.")
    args = _parse_args()
    _run_pipeline(args)
    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
