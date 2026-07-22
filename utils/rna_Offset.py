#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Generate a uniform RNA-seq offset table for a read-length range.
# Input: Minimum/maximum read lengths, expected offset, and output prefix.
# Output: RNA-seq offset table in TXT format.

"""Command-line entry point for generating RNA-seq offset tables."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

import numpy as np
import pandas as pd

from utils.ribo.ArgsParser import args_print, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate a uniform RNA-seq offset table."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output file prefix.",
    )

    offset_group = parser.add_argument_group("Offset arguments")
    offset_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=25,
        type=int,
        help="Minimum read length (default: %(default)s nt).",
    )
    offset_group.add_argument(
        "-M",
        "--max",
        dest="max",
        default=151,
        type=int,
        help="Maximum read length (default: %(default)s nt).",
    )
    offset_group.add_argument(
        "-e",
        "--exp_offset",
        dest="exp_offset",
        default=12,
        type=int,
        help="Expected RNA-seq read offset (default: %(default)s nt).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    if args.min <= 0:
        raise ValueError("-m must be > 0.")
    if args.max < args.min:
        raise ValueError("-M must be >= -m.")
    if args.exp_offset < 0:
        raise ValueError("-e must be >= 0.")


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


def _build_offset_table(args: Namespace) -> pd.DataFrame:
    """Build the RNA-seq offset table."""
    read_lengths = np.arange(args.min, args.max + 1, dtype=int)
    offset_table = pd.DataFrame(index=read_lengths)
    offset_table["length"] = read_lengths
    offset_table["frame0"] = args.exp_offset
    offset_table["rpfs0"] = 0
    offset_table["frame1"] = args.exp_offset + 1
    offset_table["rpfs1"] = 0
    offset_table["frame2"] = args.exp_offset + 2
    offset_table["rpfs2"] = 0
    offset_table["rpfs"] = 0
    offset_table["p_site"] = args.exp_offset + 1
    offset_table["periodicity"] = 100
    offset_table["ribo"] = "first"
    return offset_table


def _run_rna_offset_pipeline(args: Namespace) -> None:
    """Generate and write the RNA-seq offset table."""
    _print_step(2, "Build the RNA-seq offset table.")
    offset_table = _build_offset_table(args)

    _print_step(3, "Output the RNA-seq offset table.")
    output_file = f"{args.output}_offset.txt"
    offset_table.to_csv(output_file, sep="\t", index=False)


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rna_Offset."""
    now_time()
    print("\nGenerate a uniform RNA-seq offset table.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_rna_offset_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
