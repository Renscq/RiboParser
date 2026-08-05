#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-04
# Version: 0.2.8.19
# Function: Convert RNA-seq BAM/SAM alignments to multi-sample-compatible self-contained read density JSONL.
# Input: RiboParser norm TXT or genePred annotation, transcript FASTA, BAM/SAM alignment, and RNA-seq read offset table.
# Output: Gzip-compressed compact RNA-seq density JSONL, summary JSON, and optional legacy TXT table.

"""Command-line entry point for compact RNA-seq read density generation."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import RNA

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
            "Convert RNA-seq BAM/SAM alignments to multi-sample-compatible "
            "self-contained read density JSONL."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-t",
        "--transcript",
        dest="transcript",
        required=True,
        type=str,
        help="Input transcript annotation file in RiboParser norm TXT or genePred format.",
    )
    required_group.add_argument(
        "-s",
        "--sequence",
        dest="sequence",
        required=True,
        type=str,
        help="Input transcript sequence file in FASTA format.",
    )
    required_group.add_argument(
        "-b",
        "--bam",
        dest="bam",
        required=True,
        type=str,
        help="Input transcriptome-aligned mapping file in BAM or SAM format.",
    )
    required_group.add_argument(
        "-p",
        "--psite",
        dest="psite",
        required=True,
        type=str,
        help="Input RNA-seq read offset file in TXT format.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help=(
            "Output prefix. Default JSON output is prefix + '_rna.jsonl.gz'. "
            "The JSON uses a samples field to support downstream multi-sample merging."
        ),
    )

    filtering_group = parser.add_argument_group("Filtering arguments")
    filtering_group.add_argument(
        "-l",
        "--longest",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Only retain the transcript with longest CDS of each gene. "
            "Recommended: True. Default: %(default)s."
        ),
    )
    filtering_group.add_argument(
        "-m",
        "--min",
        dest="min",
        required=False,
        type=int,
        default=25,
        help="Minimum read length to keep. Default: %(default)s nt.",
    )
    filtering_group.add_argument(
        "-M",
        "--max",
        dest="max",
        required=False,
        type=int,
        default=150,
        help="Maximum read length to keep. Default: %(default)s nt.",
    )
    filtering_group.add_argument(
        "-r",
        "--rolling",
        dest="rolling",
        required=False,
        action="store_true",
        default=False,
        help="Assign reads with a rolling window over each aligned block. Default: %(default)s.",
    )
    filtering_group.add_argument(
        "--pe",
        dest="pair_end",
        required=False,
        action="store_true",
        default=False,
        help="Process paired-end RNA-seq alignments, keeping reverse reads. Default: %(default)s.",
    )
    filtering_group.add_argument(
        "--drop-warning",
        dest="drop_warning",
        required=False,
        action="store_true",
        default=False,
        help=(
            "Drop offset rows whose warning column is not PASS when the offset "
            "file contains a warning column. Default: %(default)s."
        ),
    )
    filtering_group.add_argument(
        "--silence",
        dest="silence",
        required=False,
        action="store_true",
        default=False,
        help="Discard warning information. Default: %(default)s.",
    )

    file_group = parser.add_argument_group("File arguments")
    file_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["json", "txt", "both"],
        default="json",
        help=(
            "Output format. 'json' writes gzip-compressed compact JSONL; "
            "'txt' writes legacy TXT; 'both' writes both. Default: %(default)s."
        ),
    )
    file_group.add_argument(
        "--density-encoding",
        dest="density_encoding",
        choices=["sparse", "dense"],
        default="sparse",
        help="Density encoding in JSONL. Sparse encoding is recommended. Default: %(default)s.",
    )

    runtime_group = parser.add_argument_group("Runtime arguments")
    runtime_group.add_argument(
        "--thread",
        dest="thread",
        type=int,
        required=False,
        default=1,
        help=(
            "Number of workers for indexed BAM parallel scanning. If BAM is not "
            "indexed, stream mode is used. Default: %(default)s."
        ),
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.transcript, args.sequence, args.bam, args.psite)
    if args.min <= 0:
        raise ValueError("-m must be > 0.")
    if args.max < args.min:
        raise ValueError("-M must be >= -m.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_density_pipeline(args: Namespace) -> None:
    """Run the RNA-seq read density generation pipeline."""
    rna_attr = RNA.RNA(args)

    step_print(2, "Import the RNA-seq read offset.")
    rna_attr.read_offset()

    step_print(3, "Import the transcripts annotation.")
    rna_attr.read_transcript()
    rna_attr.check_transcript()

    step_print(4, "Import the BAM/SAM file and count read density.")
    rna_attr.read_bam()

    step_print(5, "Output the RNA-seq read density.")
    rna_attr.output_density()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rna_Density."""
    now_time()
    title_print('Convert reads to multi-sample-compatible compact RNA-seq read density.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_density_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
