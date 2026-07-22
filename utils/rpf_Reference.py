#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Build normalized RiboParser reference files with enhanced transcript annotation.
# Input: Genome FASTA and GTF/GFF annotation files.
# Output: Normalized genePred, GTF, norm TXT, mRNA FASTA, and CDS FASTA files.

"""Command-line entry point for building RiboParser reference files."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

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
        description="Build normalized RiboParser reference files."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-g",
        dest="genome",
        required=True,
        type=str,
        help="Input genome sequence file in FASTA format.",
    )
    required_group.add_argument(
        "-t",
        dest="gtf",
        required=True,
        type=str,
        help="Input transcript annotation file in GTF or GFF/GFF3 format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix for normalized reference files.",
    )

    annotation_group = parser.add_argument_group("Annotation arguments")
    annotation_group.add_argument(
        "-u",
        dest="utr",
        default=0,
        type=int,
        help="Pseudo-UTR length added to leaderless transcripts (default: %(default)s nt).",
    )
    annotation_group.add_argument(
        "-c",
        dest="coding",
        action="store_true",
        default=False,
        help="Only retain protein-coding transcripts (default: %(default)s).",
    )
    annotation_group.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        default=False,
        help=(
            "Only retain the transcript with the longest CDS per gene "
            "(default: %(default)s)."
        ),
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-w",
        dest="whole",
        action="store_true",
        default=False,
        help="Output the complete intermediate annotation table (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.genome, args.gtf)
    if args.utr < 0:
        raise ValueError("-u must be >= 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_reference_pipeline(args: Namespace) -> None:
    """Run the RiboParser reference-building workflow."""
    from utils.ribo.GenePred import GenePred

    reference = GenePred(args)

    step_print(2, "Import the genome sequence.")
    reference.read_genome()

    step_print(3, "Format the GTF/GFF transcript annotation.")
    reference.gtf2gp()
    reference.read_genepred()
    reference.get_rep_transcript()
    reference.add_utr()

    step_print(4, "Output normalized transcript annotation.")
    reference.write_txt()
    reference.gp2gtf()

    step_print(5, "Retrieve mRNA and CDS sequences.")
    reference.get_seq()
    reference.write_seq()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Reference."""
    now_time()
    title_print('Build the reference files for RiboParser.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_reference_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
