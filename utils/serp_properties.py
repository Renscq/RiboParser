#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate codon-usage and protein physicochemical properties.
# Input: Coding sequence file in FASTA format.
# Output: Codon-usage and protein-property tables.

"""Command-line entry point for SeRP sequence-property analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence
import textwrap

from utils.ribo.ArgsParser import args_print, file_check, now_time
from utils.serp import Properties


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Calculate codon-usage and protein properties.",
        epilog=textwrap.dedent(
            """\
            Calculated properties include:
              - CAI, RSCU, and codon occurrence frequency
              - GRAVY
              - Secondary-structure fractions
              - Flexibility
              - Instability index
              - Isoelectric point
            """
        ),
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-f",
        dest="fasta",
        required=True,
        type=str,
        help="Input coding sequence file in FASTA format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output file prefix.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.fasta)


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


def _run_properties_pipeline(args: Namespace) -> None:
    """Run the sequence-property workflow."""
    _print_step(2, "Import coding sequences.")
    sequence = Properties.Sequence(args)

    _print_step(3, "Calculate gene-level codon usage.")
    sequence.create_codon_table()
    sequence.calc_gene_codon_usage()

    _print_step(4, "Calculate whole-dataset codon usage.")
    sequence.calc_whole_codon_usage()

    _print_step(5, "Calculate protein physicochemical properties.")
    sequence.protein_analysis()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for serp_properties."""
    now_time()
    print("\nCalculate sequence and protein properties.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_properties_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
