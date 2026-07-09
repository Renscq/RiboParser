#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.7.7
# Function: Draw 3-nt periodicity plots from RPF density files.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Periodicity summary table and periodicity figures.

"""Command-line entry point for 3-nt periodicity analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Periodicity
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(description="Draw 3-nt periodicity plots from RPF density files.")

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

    parser.add_argument(
        "-t",
        dest="transcript",
        required=False,
        type=str,
        default=None,
        help="Optional transcript filter table in TXT format.",
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=50,
        help="Retain transcripts with more than this minimum RPF count. Default: %(default)s.",
    )
    parser.add_argument(
        "--tis",
        dest="tis",
        required=False,
        type=int,
        default=0,
        help="Number of codons after TIS to discard. Default: %(default)s AA.",
    )
    parser.add_argument(
        "--tts",
        dest="tts",
        required=False,
        type=int,
        default=0,
        help="Number of codons before TTS to discard. Default: %(default)s AA.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)


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


def _run_periodicity_pipeline(args: Namespace) -> None:
    """Run the 3-nt periodicity analysis workflow."""
    rpfs = Periodicity.Periodicity(args)

    _print_step(2, "Import the RPFs file.")
    rpfs.import_rpf()

    _print_step(3, "Calculate the 3-nt periodicity.")
    rpfs.calc_3nt_period()

    _print_step(4, "Output the 3-nt periodicity.")
    rpfs.output_meta()

    _print_step(5, "Draw the 3-nt periodicity plots.")
    rpfs.draw_3nt_period_count()
    rpfs.draw_3nt_period_ratio()
    rpfs.draw_3nt_period_stacked()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Periodicity."""
    now_time()
    print("\nDraw the periodicity plot.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_periodicity_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
