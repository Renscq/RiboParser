#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Convert BAM alignments to strand-resolved bedGraph or WIG tracks.
# Input: BAM alignment file and RiboParser P-site offset table.
# Output: Strand-resolved bedGraph or WIG density tracks.

"""Command-line entry point for converting BAM alignments to density tracks."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Bam2Wig
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Convert BAM alignments to bedGraph or WIG density tracks."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="Input BAM alignment file.",
    )
    required_group.add_argument(
        "-p",
        dest="psite",
        required=True,
        type=str,
        help="Input P-site offset table.",
    )

    alignment_group = parser.add_argument_group("Alignment arguments")
    alignment_group.add_argument(
        "-t",
        dest="times",
        default=3,
        type=int,
        help="Maximum number of reported alignments retained per read (default: %(default)s).",
    )
    alignment_group.add_argument(
        "--second",
        dest="secondary",
        action="store_true",
        default=False,
        help="Retain secondary alignments (default: %(default)s).",
    )
    alignment_group.add_argument(
        "--supply",
        dest="supplementary",
        action="store_true",
        default=False,
        help="Retain supplementary alignments (default: %(default)s).",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-f",
        dest="format",
        choices=["bedgraph", "wig"],
        default="bedgraph",
        type=str,
        help="Output track format (default: %(default)s).",
    )
    output_group.add_argument(
        "-n",
        dest="norm",
        action="store_true",
        default=False,
        help="Normalize density values to RPM (default: %(default)s).",
    )
    output_group.add_argument(
        "-m",
        dest="merge",
        action="store_true",
        default=False,
        help="Merge strand-resolved output tracks (default: %(default)s).",
    )
    output_group.add_argument(
        "-o",
        dest="output",
        default=None,
        type=str,
        help="Output prefix. Default: BAM basename.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.bam, args.psite)
    if args.times < 1:
        raise ValueError("-t must be >= 1.")


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


def _run_bam_to_track_pipeline(args: Namespace) -> None:
    """Run the BAM-to-track conversion workflow."""
    converter = Bam2Wig.Bam2Wig(args)

    _print_step(2, "Import the P-site offset table.")
    converter.read_offset()

    _print_step(3, "Configure alignment multiplicity tags.")
    converter.set_tag_num()

    _print_step(4, "Import BAM alignments.")
    converter.import_bam()

    _print_step(5, "Convert alignments to density tables.")
    converter.convert_dict_to_dataframe()
    converter.norm_rpm()

    _print_step(6, "Output strand-resolved density tracks.")
    converter.output_bed()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Bam2bw."""
    now_time()
    print("\nConvert BAM alignments to density tracks.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_bam_to_track_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
