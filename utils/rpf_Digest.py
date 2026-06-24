#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Rensc
# Date: 2026-06-24
# Version: 0.2.7.2
# Function: This script is used to detect the digestion sites.
# Input: Command-line arguments and input files specified by the user.
# Output: Digestion-site tables, figures, and sequence-logo files.

"""Command-line entry point for rpf_Digest."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="This script is used to detect the digestion sites."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-t",
        dest="transcript",
        required=True,
        type=str,
        help="the name of input transcript file in TXT format.",
    )
    required_group.add_argument(
        "-s",
        dest="sequence",
        required=True,
        type=str,
        help="the name of input transcript sequence file in FA format.",
    )
    required_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="the name of mapping file in BAM format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="the name of output file. (prefix + _digestion_sites.txt)",
    )

    parser.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help=(
            "only retain the transcript with longest CDS of each gene "
            "(default: %(default)s). Recommended: True"
        ),
    )
    parser.add_argument(
        "--scale",
        dest="scale",
        action="store_true",
        required=False,
        default=False,
        help="scale the motif matrix (default: %(default)s).",
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=20,
        help="the minimum reads length to keep (default: %(default)s nt).",
    )
    parser.add_argument(
        "-M",
        dest="max",
        required=False,
        type=int,
        default=100,
        help="the maximum reads length to keep (default: %(default)s nt).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.transcript, args.sequence, args.bam)

    if args.min <= 0:
        raise ValueError("--min/-m must be greater than 0.")

    if args.max < args.min:
        raise ValueError("--max/-M must be greater than or equal to --min/-m.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print a standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_digest_pipeline(args: Namespace) -> None:
    """Run digestion-site detection pipeline."""
    from utils.ribo.Digestion import Ribo

    ribo_attr = Ribo(args)

    _print_step(1, "Import the annotation of transcripts.")
    ribo_attr.read_transcript()

    _print_step(2, "Detect the digestion sites.")
    ribo_attr.get_digest_sites()

    _print_step(3, "Output the digestion sites.")
    ribo_attr.output_digest_sites()

    _print_step(4, "Draw the heatmap of digestion sites.")
    ribo_attr.digestion_plot()

    _print_step(5, "Output digestion-site motif counts.")
    ribo_attr.output_counts()

    _print_step(6, "Draw the sequence logo of digestion sites.")
    ribo_attr.seq_logo_plot2()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Digest."""
    args = _parse_args(argv)

    now_time()
    print("\nDetect the digestion sites.", flush=True)

    _run_digest_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
