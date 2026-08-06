#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.002
# Function: Calculate sample-specific codon occupancy.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Codon occupancy tables, outlier records, summary JSON, and figures.

"""Command-line entry point for RPF codon occupancy analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Occupancy

from utils.ribo.ArgsParser import (
    args_print,
    complete_print,
    file_check,
    now_time,
    step_print,
    title_print,
)


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate sample-specific codon occupancy from JSONL or TXT RPF "
            "density files."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    input_group = parser.add_argument_group("Filtering arguments")
    input_group.add_argument(
        "-l",
        "--list",
        dest="list",
        type=str,
        default=None,
        help=(
            "Optional transcript filter table in TXT format. The transcript_id "
            "column is used when available; otherwise the first column is used."
        ),
    )
    input_group.add_argument(
        "-m",
        "--min",
        dest="min",
        type=int,
        default=30,
        help=(
            "Minimum sample-specific CDS RPF count required for a transcript. "
            "Default: %(default)s."
        ),
    )
    input_group.add_argument(
        "--tis",
        dest="tis",
        type=int,
        default=15,
        help="Discard this number of CDS codons after the TIS. Default: %(default)s AA.",
    )
    input_group.add_argument(
        "--tts",
        dest="tts",
        type=int,
        default=5,
        help="Discard this number of CDS codons before the TTS. Default: %(default)s AA.",
    )

    calculation_group = parser.add_argument_group("Occupancy calculation arguments")
    calculation_group.add_argument(
        "-s",
        "--site",
        dest="site",
        choices=["E", "P", "A"],
        type=str,
        default="P",
        help="Ribosomal site used for occupancy calculation. Default: %(default)s.",
    )
    calculation_group.add_argument(
        "-f",
        "--frame",
        dest="frame",
        choices=["0", "1", "2", "all"],
        type=str,
        default="all",
        help="Reading frame used for occupancy calculation. Default: %(default)s.",
    )
    calculation_group.add_argument(
        "--thread",
        dest="thread",
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
        default=False,
        help=(
            "Normalize density to RPM before occupancy calculation. The occupancy "
            "ratio is unchanged within a sample, but normalized raw density is "
            "retained in detailed outputs. Default: %(default)s."
        ),
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier",
        dest="remove_outlier",
        action="store_true",
        default=False,
        help=(
            "Remove isolated extreme raw-density pileups before occupancy "
            "calculation. Default: %(default)s."
        ),
    )
    outlier_group.add_argument(
        "--outlier-iqr",
        dest="outlier_iqr",
        type=float,
        default=8.0,
        help="IQR multiplier for global extreme-pileup detection. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-window",
        dest="outlier_window",
        type=int,
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
        default=10.0,
        help="Minimum fold over local background for outlier removal. Default: %(default)s.",
    )

    output_group = parser.add_argument_group("Plotting arguments")
    output_group.add_argument(
        "--scale",
        dest="scale",
        choices=["none", "minmax", "zscore"],
        type=str,
        default="minmax",
        help="Sample-wise scaling for relative codon occupancy. Default: %(default)s.",
    )
    output_group.add_argument(
        "--plot-transform",
        dest="plot_transform",
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        type=str,
        default="none",
        help=(
            "Transform absolute occupancy for plotting only; output tables are not "
            "changed. 'log' is an alias of log1p. Default: %(default)s."
        ),
    )
    output_group.add_argument(
        "--rankplot-ncol",
        dest="rankplot_ncol",
        required=False,
        type=int,
        default=1,
        help=(
            "Number of sample panels per row in the codon rank plot. Every codon is "
            "labeled on the x-axis (rotated 90 degrees), so 1 gives one wide panel "
            "per row and 2 gives two narrower panels per row. Default: %(default)s."
        ),
    )
    output_group.add_argument(
        "--all",
        dest="all",
        action="store_true",
        default=False,
        help=(
            "Output detailed transcript-position and transcript-codon occupancy "
            "tables. Default: %(default)s."
        ),
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
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.outlier_iqr <= 0:
        raise ValueError("--outlier-iqr must be > 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 1:
        raise ValueError("--outlier-local-fold must be > 1.")
    if args.rankplot_ncol < 1:
        raise ValueError("--rankplot-ncol must be >= 1.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_occupancy_pipeline(args: Namespace) -> None:
    """Run the codon occupancy workflow."""
    occupancy = Occupancy.Occupancy(args)

    step_print(2, "Import the RPF density file.")
    occupancy.import_rpf()

    step_print(3, "Calculate sample-specific codon occupancy.")
    occupancy.calculate_occupancy()

    step_print(4, "Output codon occupancy tables.")
    occupancy.output_tables()

    step_print(5, "Draw codon occupancy figures.")
    occupancy.draw_occupancy_corr()
    occupancy.draw_occupancy_heatmap()
    occupancy.draw_occupancy_rankplot()

    step_print(6, "Write occupancy summary.")
    occupancy.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Occupancy."""
    now_time()
    title_print('Calculate codon occupancy.')
    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)
    _run_occupancy_pipeline(args)
    complete_print()
    now_time()


if __name__ == "__main__":
    main()
