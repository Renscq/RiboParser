#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate transcript-level CDS RPF coverage percentages.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Gene coverage percentage tables, summary tables, and figures.

"""Command-line entry point for transcript-level RPF coverage analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence
from pathlib import Path

from utils.ribo import Percentage
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate the percentage of valid CDS codons covered by RPFs "
            "for each transcript and sample."
        )
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

    input_group = parser.add_argument_group("Input filtering arguments")
    input_group.add_argument(
        "-t",
        "--transcript",
        dest="transcript",
        default=None,
        type=str,
        help=(
            "Optional transcript annotation or ID list. The transcript_id "
            "column is used when present; otherwise the first column is used."
        ),
    )
    input_group.add_argument(
        "-f",
        "--frame",
        dest="frame",
        choices=["0", "1", "2", "all"],
        default="all",
        type=str,
        help="Reading frame used for coverage calculation (default: %(default)s).",
    )
    input_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=50,
        type=int,
        help=(
            "Minimum sample-specific CDS RPF count required for a transcript "
            "(default: %(default)s)."
        ),
    )
    input_group.add_argument(
        "--tis",
        dest="tis",
        default=0,
        type=int,
        help="Number of CDS codons removed after TIS (default: %(default)s).",
    )
    input_group.add_argument(
        "--tts",
        dest="tts",
        default=0,
        type=int,
        help="Number of CDS codons removed before TTS (default: %(default)s).",
    )

    calculation_group = parser.add_argument_group("Calculation arguments")
    calculation_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        default=False,
        help=(
            "Use RPM instead of raw RPF count in abundance plots and summaries "
            "(default: %(default)s)."
        ),
    )

    plotting_group = parser.add_argument_group("Plotting arguments")
    plotting_group.add_argument(
        "--fig-format",
        dest="fig_format",
        choices=["png", "pdf", "both"],
        default="png",
        type=str,
        help="Figure output format (default: %(default)s).",
    )
    plotting_group.add_argument(
        "--dpi",
        dest="dpi",
        default=300,
        type=int,
        help="PNG resolution (default: %(default)s).",
    )
    plotting_group.add_argument(
        "--font-size",
        dest="font_size",
        default=11.0,
        type=float,
        help="Base figure font size (default: %(default)s).",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        type=str,
        help="Output prefix. Default: input RPF filename without density suffix.",
    )

    return parser


def _default_output_prefix(rpf_file: str) -> str:
    """Return a stable output prefix derived from the RPF filename."""
    file_name = Path(rpf_file).name
    lower_name = file_name.lower()

    for suffix in (
        ".jsonl.gz",
        ".json.gz",
        ".jsonl",
        ".json",
        ".tsv",
        ".txt",
        ".tab",
    ):
        if lower_name.endswith(suffix):
            return file_name[: -len(suffix)]
    return str(Path(file_name).stem)


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)

    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if args.dpi < 72:
        raise ValueError("--dpi must be >= 72.")
    if args.font_size <= 0:
        raise ValueError("--font-size must be > 0.")

    if args.output is None:
        args.output = _default_output_prefix(args.rpf)


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        _validate_args(args)
    except ValueError as error:
        parser.error(str(error))

    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print a standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_percentage_pipeline(args: Namespace) -> None:
    """Run the transcript-level CDS RPF coverage workflow."""
    percentage = Percentage.Percentage(args)

    _print_step(2, "Import the RPF density file.")
    percentage.read_rpf()
    percentage.import_gene()

    _print_step(3, "Calculate transcript-level CDS coverage.")
    percentage.calc_density_percent()

    _print_step(4, "Draw RPF coverage and abundance figures.")
    percentage.draw_rpf_histogram()
    percentage.draw_rpf_boxplot()

    _print_step(5, "Output RPF coverage tables.")
    percentage.output_density_percent()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Percent."""
    now_time()
    print("\nCalculate transcript-level CDS RPF coverage.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_percentage_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
