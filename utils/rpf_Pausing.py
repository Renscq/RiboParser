#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.005
# Function: Calculate sample-specific codon pausing scores.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Pausing score tables, outlier records, summary JSON, and codon-level figures.

"""Command-line entry point for RPF codon pausing analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Pausing
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate sample-specific codon pausing scores from JSONL or TXT "
            "RPF density files."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    input_group = parser.add_argument_group("Input filtering arguments")
    input_group.add_argument(
        "-l",
        "--list",
        dest="list",
        required=False,
        type=str,
        default=None,
        help=(
            "Optional transcript filter table in TXT format. The transcript_id "
            "column is used when available; otherwise the first column is used."
        ),
    )
    input_group.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=50,
        help=(
            "Minimum sample-specific CDS RPF count required for a transcript. "
            "Default: %(default)s."
        ),
    )
    input_group.add_argument(
        "--tis",
        dest="tis",
        required=False,
        type=int,
        default=10,
        help="Discard this number of CDS codons after the TIS. Default: %(default)s AA.",
    )
    input_group.add_argument(
        "--tts",
        dest="tts",
        required=False,
        type=int,
        default=5,
        help="Discard this number of CDS codons before the TTS. Default: %(default)s AA.",
    )
    calculation_group = parser.add_argument_group("Pausing calculation arguments")
    calculation_group.add_argument(
        "-s",
        dest="site",
        choices=["E", "P", "A"],
        required=False,
        type=str,
        default="P",
        help="Ribosomal site used for pausing calculation. Default: %(default)s.",
    )
    calculation_group.add_argument(
        "-f",
        dest="frame",
        choices=["0", "1", "2", "all"],
        required=False,
        type=str,
        default="all",
        help="Reading frame used for pausing calculation. Default: %(default)s.",
    )
    calculation_group.add_argument(
        "-b",
        dest="background",
        required=False,
        type=int,
        default=0,
        help=(
            "Number of neighboring codons on each side used as local background. "
            "The focal codon is excluded. Use 0 for the mean CDS density of each "
            "gene and sample. Default: %(default)s."
        ),
    )
    calculation_group.add_argument(
        "--thread",
        dest="thread",
        required=False,
        type=int,
        default=1,
        help=(
            "Number of sample-level worker threads. Values larger than the sample "
            "count are capped automatically. Default: %(default)s."
        ),
    )
    calculation_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Normalize density to RPM before pausing calculation. Ratios are "
            "mathematically unchanged within a sample but normalized density is "
            "retained in exported intermediate tables. Default: %(default)s."
        ),
    )
    calculation_group.add_argument(
        "--ind",
        dest="individual",
        action="store_true",
        required=False,
        default=False,
        help=(
            "For valid-codon summaries, require positive RPF density independently "
            "within each sample. Default: %(default)s."
        ),
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier",
        dest="remove_outlier",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Remove isolated extreme raw-density pileups before background and "
            "pausing-score calculation. Default: %(default)s."
        ),
    )
    outlier_group.add_argument(
        "--outlier-iqr",
        dest="outlier_iqr",
        type=float,
        required=False,
        default=8.0,
        help="IQR multiplier for global extreme-pileup detection. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-window",
        dest="outlier_window",
        type=int,
        required=False,
        default=5,
        help=(
            "Neighboring codons on each side used to confirm an isolated outlier. "
            "Default: %(default)s."
        ),
    )
    outlier_group.add_argument(
        "--outlier-local-fold",
        dest="outlier_local_fold",
        type=float,
        required=False,
        default=10.0,
        help="Minimum fold over local background for outlier removal. Default: %(default)s.",
    )

    output_group = parser.add_argument_group("Output and plotting arguments")
    output_group.add_argument(
        "--scale",
        dest="scale",
        choices=["none", "minmax", "zscore"],
        required=False,
        type=str,
        default="minmax",
        help="Sample-wise scaling for relative codon pausing scores. Default: %(default)s.",
    )
    output_group.add_argument(
        "--plot-transform",
        dest="plot_transform",
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        required=False,
        type=str,
        default="none",
        help=(
            "Transform absolute pausing scores for plotting only; output tables are "
            "not changed. 'log' is an alias of log1p. Default: %(default)s."
        ),
    )
    output_group.add_argument(
        "--all",
        dest="all",
        action="store_true",
        required=False,
        default=False,
        help=("Output detailed transcript-position and transcript-codon pausing "
              "tables. Default: %(default)s."),
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)

    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.tis < 0:
        raise ValueError("--tis must be >= 0.")
    if args.tts < 0:
        raise ValueError("--tts must be >= 0.")
    if args.background < 0:
        raise ValueError("-b/--background must be >= 0.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.outlier_iqr <= 0:
        raise ValueError("--outlier-iqr must be > 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 1:
        raise ValueError("--outlier-local-fold must be > 1.")


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


def _run_pausing_pipeline(args: Namespace) -> None:
    """Run the codon pausing analysis workflow."""
    pausing = Pausing.Pausing(args)

    _print_step(2, "Import the RPF density file.")
    pausing.import_rpf()

    _print_step(3, "Calculate sample-specific pausing scores.")
    pausing.calculate_pausing()

    _print_step(4, "Output pausing score tables.")
    pausing.output_cds_pausing()
    pausing.output_cds_codon_pausing()
    pausing.output_sum_codon_pausing()
    pausing.output_all_pausing()

    _print_step(5, "Draw codon pausing figures.")
    pausing.draw_pausing_corr()
    pausing.draw_codon_pausing_heatmap()
    pausing.draw_codon_rank_plot()

    _print_step(6, "Write pausing summary.")
    pausing.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Pausing."""
    now_time()
    print("\nCalculate sample-specific codon pausing scores.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_pausing_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
