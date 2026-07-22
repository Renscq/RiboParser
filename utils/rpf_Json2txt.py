#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Convert RiboParser RPF density data from JSONL to TXT format.
# Input: RPF density file in JSONL or JSONL.GZ format.
# Output: Frame-resolved RPF density table in TXT format.

"""Command-line entry point for converting RPF density JSONL to TXT."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Ribo

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
        description="Convert RPF density data from JSONL to TXT format."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-i",
        "--input",
        dest="json_file",
        required=True,
        type=str,
        help="Input RPF density file in JSONL or JSONL.GZ format.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output_txt",
        required=True,
        type=str,
        help="Output RPF density file in TXT format.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.json_file)


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_conversion_pipeline(args: Namespace) -> None:
    """Convert the input JSONL density file to TXT format."""
    step_print(2, "Convert RPF density JSONL to TXT.")
    Ribo.Ribo.json_to_txt(
        json_file=args.json_file,
        output_txt=args.output_txt,
    )


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Json2txt."""
    now_time()
    title_print('Convert RPF density JSONL to TXT.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_conversion_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
