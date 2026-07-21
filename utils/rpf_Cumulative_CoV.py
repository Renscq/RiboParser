#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate cumulative coefficient of variation along CDS regions.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Cumulative CoV tables, summary JSON, outlier table, and figures.

"""Command-line entry point for cumulative CoV analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Cumulative_CoV
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Calculate cumulative coefficient of variation along CDS regions.",
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density file in JSONL, JSONL.GZ, TXT, or TSV format.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output file prefix.",
    )

    filtering_group = parser.add_argument_group("Input filtering arguments")
    filtering_group.add_argument(
        "-l",
        "--list",
        dest="list",
        default=None,
        type=str,
        help="Optional transcript ID list.",
    )
    filtering_group.add_argument(
        "-s",
        "--site",
        dest="site",
        choices=["E", "P", "A"],
        default="P",
        help="Ribosome site used for coordinate assignment.",
    )
    filtering_group.add_argument(
        "-f",
        "--frame",
        dest="frame",
        choices=["0", "1", "2", "all"],
        default="all",
        help="Reading frame used for cumulative CoV calculation.",
    )
    filtering_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=0,
        type=float,
        help="Minimum sample-specific CDS RPF count per transcript.",
    )
    filtering_group.add_argument(
        "--tis",
        dest="tis",
        default=0,
        type=int,
        help="Number of codons discarded after TIS.",
    )
    filtering_group.add_argument(
        "--tts",
        dest="tts",
        default=0,
        type=int,
        help="Number of codons discarded before TTS.",
    )
    filtering_group.add_argument(
        "--min-positions",
        dest="min_positions",
        default=10,
        type=int,
        help="Minimum retained positions required per transcript.",
    )

    calculation_group = parser.add_argument_group("Cumulative CoV arguments")
    calculation_group.add_argument(
        "-t",
        "--trim",
        dest="trim",
        default=150,
        type=int,
        help="Maximum CDS-relative position included in meta analysis.",
    )
    calculation_group.add_argument(
        "--resolution",
        dest="resolution",
        choices=["nucleotide", "codon"],
        default="nucleotide",
        help="Position resolution used for cumulative CoV calculation.",
    )
    calculation_group.add_argument(
        "--ddof",
        dest="ddof",
        choices=[0, 1],
        default=1,
        type=int,
        help="Delta degrees of freedom used for cumulative standard deviation.",
    )
    calculation_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        default=False,
        help="Convert density to RPM before reporting mean and SD.",
    )
    calculation_group.add_argument(
        "--thread",
        dest="thread",
        default=1,
        type=int,
        help="Number of sample-level worker threads.",
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier",
        dest="remove_outlier",
        action="store_true",
        default=False,
        help="Remove extreme local RPF pileups.",
    )
    outlier_group.add_argument(
        "--outlier-iqr",
        dest="outlier_iqr",
        default=8.0,
        type=float,
        help="IQR multiplier for the global log1p-density cutoff.",
    )
    outlier_group.add_argument(
        "--outlier-window",
        dest="outlier_window",
        default=5,
        type=int,
        help="Neighboring positions used on each side of a candidate outlier.",
    )
    outlier_group.add_argument(
        "--outlier-local-fold",
        dest="outlier_local_fold",
        default=10.0,
        type=float,
        help="Minimum fold above local background required for removal.",
    )

    plotting_group = parser.add_argument_group("Output and plotting arguments")
    plotting_group.add_argument(
        "--plot-stat",
        dest="plot_stat",
        choices=["median", "mean"],
        default="median",
        help="Center statistic used for the cumulative CoV meta curve.",
    )
    plotting_group.add_argument(
        "--plot-transform",
        dest="plot_transform",
        choices=["none", "sqrt", "log1p", "log2", "log10"],
        default="none",
        help="Transformation applied only to plotted CoV values.",
    )
    plotting_group.add_argument(
        "--ci-low",
        dest="ci_low",
        default=0.25,
        type=float,
        help="Lower transcript quantile used for the curve ribbon.",
    )
    plotting_group.add_argument(
        "--ci-high",
        dest="ci_high",
        default=0.75,
        type=float,
        help="Upper transcript quantile used for the curve ribbon.",
    )
    plotting_group.add_argument(
        "--gene-fig",
        dest="gene_fig",
        choices=["png", "pdf", "both"],
        default="png",
        help="Output format for per-transcript cumulative CoV figures.",
    )
    plotting_group.add_argument(
        "--all",
        dest="all",
        action="store_true",
        default=False,
        help="Output transcript-position cumulative CoV details.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
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
    if args.outlier_iqr < 0:
        raise ValueError("--outlier-iqr must be >= 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 0:
        raise ValueError("--outlier-local-fold must be > 0.")
    if not 0 <= args.ci_low < args.ci_high <= 1:
        raise ValueError("Require 0 <= --ci-low < --ci-high <= 1.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print a standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_cov_pipeline(args: Namespace) -> None:
    """Run cumulative CoV analysis."""
    analysis = Cumulative_CoV.CumulativeCoV(args)

    _print_step(2, "Import the RPF density file.")
    analysis.import_rpf()

    _print_step(3, "Calculate sample-specific cumulative CoV.")
    analysis.calculate_cumulative_cov()

    _print_step(4, "Draw cumulative CoV figures.")
    analysis.draw_cumulative_cov()
    analysis.draw_transcript_count()
    analysis.draw_gene_cumulative_cov()

    _print_step(5, "Output cumulative CoV results.")
    analysis.output_results()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Cumulative_CoV."""
    now_time()
    print("\nCalculate cumulative coefficient of variation.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_cov_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
