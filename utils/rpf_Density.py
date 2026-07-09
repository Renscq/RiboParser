#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.009
# Function: Convert Ribo-seq BAM/SAM alignments to multi-sample-compatible self-contained P-site density JSONL.
# Input: RiboParser norm TXT or genePred annotation, transcript FASTA, BAM/SAM alignment, and P-site offset table.
# Output: Gzip-compressed compact RPF density JSONL, summary JSON, and optional legacy TXT table.

"""Command-line entry point for compact RPF density generation."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Ribo
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Convert Ribo-seq BAM/SAM alignments to multi-sample-compatible "
            "self-contained P-site density JSONL."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-t",
        dest="transcript",
        required=True,
        type=str,
        help="Input transcript annotation file in RiboParser norm TXT or genePred format.",
    )
    required_group.add_argument(
        "-s",
        dest="sequence",
        required=True,
        type=str,
        help="Input transcript sequence file in FASTA format.",
    )
    required_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="Input transcriptome-aligned mapping file in BAM or SAM format.",
    )
    required_group.add_argument(
        "-p",
        dest="psite",
        required=True,
        type=str,
        help="Input P-site offset file in TXT format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help=(
            "Output prefix. Default JSON output is prefix + '_rpf.jsonl.gz'. "
            "The JSON uses a samples field to support downstream multi-sample merging."
        ),
    )

    parser.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Only retain the transcript with longest CDS of each gene. "
            "Recommended: True. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=27,
        help="Minimum read length to keep. Default: %(default)s nt.",
    )
    parser.add_argument(
        "-M",
        dest="max",
        required=False,
        type=int,
        default=33,
        help="Maximum read length to keep. Default: %(default)s nt.",
    )
    parser.add_argument(
        "--period",
        dest="periodicity",
        required=False,
        type=float,
        default=40,
        help="Minimum 3-nt periodicity to keep. Default: %(default)s.",
    )
    parser.add_argument(
        "--min-confidence",
        dest="min_confidence",
        required=False,
        type=float,
        default=50.0,
        help=(
            "Minimum offset confidence to keep when the offset file contains "
            "a confidence column. Default: %(default)s."
        ),
    )
    parser.add_argument(
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
    parser.add_argument(
        "--silence",
        dest="silence",
        required=False,
        action="store_true",
        default=False,
        help="Discard warning information. Default: %(default)s.",
    )
    parser.add_argument(
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

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["json", "txt", "both"],
        default="json",
        help=(
            "Output format. 'json' writes gzip-compressed compact JSONL; "
            "'txt' writes legacy TXT; 'both' writes both. Default: %(default)s."
        ),
    )
    output_group.add_argument(
        "--density-encoding",
        dest="density_encoding",
        choices=["sparse", "dense"],
        default="sparse",
        help="Density encoding in JSONL. Sparse encoding is recommended. Default: %(default)s.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.transcript, args.sequence, args.bam, args.psite)


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


def _run_density_pipeline(args: Namespace) -> None:
    """Run the RPF density generation pipeline."""
    ribo_attr = Ribo.Ribo(args)

    _print_step(2, "Import the P-site offset.")
    ribo_attr.read_offset()

    _print_step(3, "Import the transcripts annotation.")
    ribo_attr.read_transcript()
    ribo_attr.check_transcript()

    _print_step(4, "Import the BAM/SAM file and count P-site density.")
    ribo_attr.read_bam()

    _print_step(5, "Output the RPF density.")
    ribo_attr.output_density()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Density."""
    now_time()
    print("\nConvert reads to multi-sample-compatible compact RPF density.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_density_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
