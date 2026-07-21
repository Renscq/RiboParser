#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate and visualize meta-codon RPF density profiles.
# Input: Frame-resolved RPF density file and optional codon/transcript lists.
# Output: Meta-codon density tables, sequence tables, and figures.

"""Command-line entry point for meta-codon analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import args_print, file_check, now_time
from utils.ribo.MetaCodon import MetaCodon


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate meta-codon RPF density profiles."
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

    filtering_group = parser.add_argument_group("Input filtering arguments")
    filtering_group.add_argument(
        "-l",
        dest="list",
        default=None,
        type=str,
        help="Optional transcript ID list.",
    )
    filtering_group.add_argument(
        "-c",
        dest="codon",
        default=None,
        type=str,
        help="Optional codon list.",
    )
    filtering_group.add_argument(
        "-f",
        dest="frame",
        choices=["0", "1", "2", "all"],
        default="all",
        type=str,
        help="Reading frame used for analysis (default: %(default)s).",
    )
    filtering_group.add_argument(
        "-a",
        dest="around",
        default=20,
        type=int,
        help="Number of codons retained around the target codon (default: %(default)s).",
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
        "-u",
        dest="unique",
        action="store_true",
        default=False,
        help="Only retain uniquely occurring target codons (default: %(default)s).",
    )

    calculation_group = parser.add_argument_group("Calculation arguments")
    calculation_group.add_argument(
        "-n",
        dest="normal",
        action="store_true",
        default=False,
        help="Normalize RPF counts to RPM (default: %(default)s).",
    )
    calculation_group.add_argument(
        "-s",
        dest="scale",
        action="store_true",
        default=False,
        help="Scale meta-codon density profiles (default: %(default)s).",
    )
    calculation_group.add_argument(
        "--smooth",
        dest="smooth",
        default=None,
        type=str,
        help="Optional smoothing method.",
    )
    calculation_group.add_argument(
        "--thread",
        dest="thread",
        default=1,
        type=int,
        help="Number of worker processes (default: %(default)s).",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--fig",
        dest="fig",
        action="store_true",
        default=False,
        help="Draw meta-codon figures (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)
    if args.codon:
        file_check(args.codon)
    if args.around < 0:
        raise ValueError("-a must be >= 0.")
    if args.min < 0:
        raise ValueError("-m must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")


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


def _run_meta_codon_pipeline(args: Namespace) -> None:
    """Run the meta-codon density workflow."""
    meta_codon = MetaCodon(args)

    _print_step(2, "Import the codon list.")
    meta_codon.import_codon()

    _print_step(3, "Import the RPF density file.")
    meta_codon.import_rpf()

    _print_step(4, "Smooth the RPF density profiles.")
    meta_codon.smooth_rpf_density()

    _print_step(5, "Retrieve meta-codon density profiles.")
    meta_codon.reterieve_codon_density()

    _print_step(6, "Output meta-codon results.")
    meta_codon.output_meta_codon_density()
    meta_codon.output_meta_codon_seq()

    if args.fig:
        _print_step(7, "Draw meta-codon figures.")
        meta_codon.draw_meta_codon()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Meta_Codon."""
    now_time()
    print("\nCalculate meta-codon RPF density profiles.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_meta_codon_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
