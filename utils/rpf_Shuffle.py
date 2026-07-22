#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.001
# Function: Shuffle RPF density profiles within each transcript to generate randomized control data.
# Input: RPF density file in JSONL or TXT format and an optional transcript filter.
# Output: Shuffled density file in JSONL and/or TXT format and a summary JSON file.

"""Command-line entry point for RPF density shuffling."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Shuffle

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
            "Shuffle codon-frame RPF density within each transcript to generate "
            "reproducible randomized control data."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density file in JSONL or legacy TXT format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix. Shuffled files use prefix + '_shuffle'.",
    )

    filter_group = parser.add_argument_group("Input filtering arguments")
    filter_group.add_argument(
        "-l",
        "--list",
        dest="list",
        required=False,
        type=str,
        default=None,
        help=(
            "Optional transcript filter table. The transcript_id column is used "
            "when available; otherwise the first column is used."
        ),
    )

    shuffle_group = parser.add_argument_group("Shuffle arguments")
    shuffle_group.add_argument(
        "-s",
        "--seed",
        dest="seed",
        required=False,
        type=int,
        default=0,
        help="Global random seed. Default: %(default)s.",
    )
    shuffle_group.add_argument(
        "-i",
        "--individual",
        dest="individual",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Shuffle each sample independently within each transcript. By default, "
            "one shared position permutation is applied to all samples, preserving "
            "cross-sample site-level covariance. Default: %(default)s."
        ),
    )
    shuffle_group.add_argument(
        "--thread",
        dest="thread",
        required=False,
        type=int,
        default=1,
        help=(
            "Number of transcript-level worker threads. Results are deterministic "
            "across worker counts. Default: %(default)s."
        ),
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["auto", "json", "txt", "both"],
        default="auto",
        help=(
            "Output density format. 'auto' preserves the input format. "
            "Default: %(default)s."
        ),
    )
    output_group.add_argument(
        "--density-encoding",
        dest="density_encoding",
        choices=["auto", "sparse", "dense"],
        default="auto",
        help=(
            "JSON density encoding. 'auto' preserves each source record encoding. "
            "Default: %(default)s."
        ),
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_shuffle_pipeline(args: Namespace) -> None:
    """Run the RPF density shuffle workflow."""
    shuffler = Shuffle.Shuffle(args)

    step_print(2, "Inspect the RPF density file.")
    shuffler.import_rpf()

    step_print(3, "Shuffle transcript-level RPF density profiles.")
    shuffler.shuffle_rpfs()

    step_print(4, "Write the shuffle summary.")
    shuffler.output_rpfs()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Shuffle."""
    now_time()
    title_print('Shuffle RPF density profiles.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_shuffle_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
