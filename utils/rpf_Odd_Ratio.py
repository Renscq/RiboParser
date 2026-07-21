#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate codon-level RPF odds ratios between treatment groups.
# Input: Frame-resolved RPF density file and control/treatment sample names.
# Output: Codon odds-ratio tables, statistical summaries, and figures.

"""Command-line entry point for codon-level odds-ratio analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Odd_Ratio
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate codon-level RPF odds ratios."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="Input frame-resolved RPF density file.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output file prefix.",
    )
    required_group.add_argument(
        "-c",
        dest="control",
        required=True,
        type=str,
        help="Comma-separated control sample names.",
    )
    required_group.add_argument(
        "-t",
        dest="treat",
        required=True,
        type=str,
        help="Comma-separated treatment sample names.",
    )

    filtering_group = parser.add_argument_group("Input filtering arguments")
    filtering_group.add_argument(
        "-l",
        dest="list",
        default=None,
        type=str,
        help="Optional transcript ID list.",
    )
    filtering_group.add_argument(
        "-s",
        dest="site",
        choices=["E", "P", "A"],
        default="P",
        type=str,
        help="Ribosome site used for coordinate assignment (default: %(default)s).",
    )
    filtering_group.add_argument(
        "-f",
        dest="frame",
        choices=["0", "1", "2", "all"],
        default="all",
        type=str,
        help="Reading frame used for calculation (default: %(default)s).",
    )
    filtering_group.add_argument(
        "-m",
        dest="min",
        default=50,
        type=int,
        help="Minimum RPF count required for retained transcripts (default: %(default)s).",
    )
    filtering_group.add_argument(
        "--tis",
        dest="tis",
        default=0,
        type=int,
        help="Number of codons removed after TIS (default: %(default)s).",
    )
    filtering_group.add_argument(
        "--tts",
        dest="tts",
        default=0,
        type=int,
        help="Number of codons removed before TTS (default: %(default)s).",
    )
    filtering_group.add_argument(
        "--stop",
        dest="stop",
        action="store_true",
        default=False,
        help="Include stop codons in the analysis (default: %(default)s).",
    )

    calculation_group = parser.add_argument_group("Calculation arguments")
    calculation_group.add_argument(
        "--thread",
        dest="thread",
        default=1,
        type=int,
        help="Number of worker processes (default: %(default)s).",
    )
    calculation_group.add_argument(
        "-n",
        dest="normal",
        action="store_true",
        default=False,
        help="Normalize RPF counts to RPM (default: %(default)s).",
    )
    calculation_group.add_argument(
        "-z",
        dest="zero",
        action="store_true",
        default=False,
        help="Retain zero-density positions (default: %(default)s).",
    )
    calculation_group.add_argument(
        "--fdr",
        dest="fdr",
        choices=["bhfdr", "pvalue"],
        default="bhfdr",
        type=str,
        help="Significance field used for filtering (default: %(default)s).",
    )
    calculation_group.add_argument(
        "-v",
        dest="value",
        default=0.05,
        type=float,
        help="Significance threshold (default: %(default)s).",
    )
    calculation_group.add_argument(
        "--scale",
        dest="scale",
        choices=["zscore", "minmax"],
        default="minmax",
        type=str,
        help="Scaling method used for visualization (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.min < 0:
        raise ValueError("-m must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if not 0 < args.value <= 1:
        raise ValueError("-v must be in the interval (0, 1].")


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


def _run_odd_ratio_pipeline(args: Namespace) -> None:
    """Run the codon-level odds-ratio workflow."""
    odd_ratio = Odd_Ratio.OddRatio(args)

    _print_step(2, "Import the RPF density file.")
    odd_ratio.read_rpf()

    _print_step(3, "Build the codon contingency tables.")
    odd_ratio.make_two_dimensional_table()

    _print_step(4, "Calculate codon odds ratios.")
    odd_ratio.calc_odd_ratio()

    _print_step(5, "Perform codon-level statistical tests.")
    odd_ratio.calc_chi2_test2()

    _print_step(6, "Output codon odds-ratio results.")
    odd_ratio.output_odd_ratio()
    odd_ratio.summarize_odd_ratio()

    _print_step(7, "Draw codon odds-ratio figures.")
    odd_ratio.draw_odd_ratio_line()
    odd_ratio.draw_odd_ratio_scatter()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Odd_Ratio."""
    now_time()
    print("\nCalculate codon-level RPF odds ratios.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_odd_ratio_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
