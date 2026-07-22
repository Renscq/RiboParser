#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Retrieve selected transcript RPF density profiles.
# Input: Frame-resolved RPF density file and optional transcript ID list.
# Output: Retrieved RPF density table and optional per-transcript files.

"""Command-line entry point for retrieving selected RPF density profiles."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence


from utils.ribo.Retrieve import Retrieve

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
        description="Retrieve selected transcript RPF density profiles."
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
        required=False,
        default=None,
        type=str,
        help="Output file prefix. Default: determined from the input filename.",
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
        "-m",
        dest="min",
        default=0,
        type=int,
        help="Minimum RPF count required for retained transcripts (default: %(default)s).",
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
        "-f",
        dest="format",
        action="store_true",
        default=False,
        help="Melt frame-resolved sample columns into long format (default: %(default)s).",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-s",
        dest="split",
        action="store_true",
        default=False,
        help="Write one RPF density file per transcript (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)
    if args.min < 0:
        raise ValueError("-m must be >= 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_retrieve_pipeline(args: Namespace) -> None:
    """Run the RPF density retrieval workflow."""
    rpfs = Retrieve(args)

    step_print(2, "Retrieve the selected transcript RPF density.")
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    step_print(3, "Format the RPF density table.")
    rpfs.melt_rpf_table()

    step_print(4, "Output the retrieved RPF density table.")
    rpfs.output_rpf_table()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Retrieve."""
    now_time()
    title_print('Retrieve selected transcript RPF density profiles.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_retrieve_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
