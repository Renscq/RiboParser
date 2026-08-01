#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev003
# Function: Scan complete transcript-centric smORFs and discard invalid candidates.
# Input: Genome FASTA and genePred transcript annotation.
# Output: smORF genePred, metadata, nucleotide FASTA, and peptide FASTA files.

"""Command-line entry point for transcript-centric smORF scanning."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import (
    args_print,
    file_check,
    now_time,
    step_print,
    title_print,
)
from utils.smorf.scanner import SmORFPipeline


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser.

    Returns:
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Scan transcript-centric small open reading frames."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-g",
        "--genome",
        dest="genome",
        required=True,
        type=str,
        help="Input genome sequence file in FASTA or FASTA.GZ format.",
    )
    required_group.add_argument(
        "-a",
        "--annotation",
        dest="annotation",
        required=True,
        type=str,
        help="Input transcript annotation file in genePred format.",
    )

    scanning_group = parser.add_argument_group("ORF scanning arguments")
    scanning_group.add_argument(
        "-o",
        "--out-prefix",
        dest="out_prefix",
        default="ORF",
        type=str,
        help="Output file prefix (default: %(default)s).",
    )
    scanning_group.add_argument(
        "-p",
        "--orf-prefix",
        dest="orf_prefix",
        default="ORF",
        type=str,
        help="Generated ORF identifier prefix (default: %(default)s).",
    )
    scanning_group.add_argument(
        "-s",
        "--start-codons",
        dest="start_codons",
        default="ATG",
        type=str,
        help="Comma-separated start codons (default: %(default)s).",
    )
    scanning_group.add_argument(
        "-m",
        "--min-aa",
        dest="min_aa",
        default=8,
        type=int,
        help="Minimum peptide length (default: %(default)s aa).",
    )
    scanning_group.add_argument(
        "-M",
        "--max-aa",
        dest="max_aa",
        default=10000,
        type=int,
        help="Maximum peptide length (default: %(default)s aa).",
    )
    scanning_group.add_argument(
        "-x",
        "--scan-strand",
        dest="scan_strand",
        choices=["sense", "antisense", "both"],
        default="sense",
        type=str,
        help="Transcript orientation scanned for ORFs (default: %(default)s).",
    )
    scanning_group.add_argument(
        "-u",
        "--kozak-up",
        dest="kozak_up",
        default=6,
        type=int,
        help="Upstream Kozak-context length (default: %(default)s nt).",
    )
    scanning_group.add_argument(
        "-d",
        "--kozak-down",
        dest="kozak_down",
        default=6,
        type=int,
        help=("Downstream Kozak-context length after the start codon (default: %(default)s nt)."),
    )
    scanning_group.add_argument(
        "-I",
        "--include-stop",
        dest="include_stop",
        action="store_true",
        default=False,
        help=("Retain the terminal stop symbol in peptide sequences (default: %(default)s)."),
    )

    overlap_group = parser.add_argument_group("Overlap arguments")
    overlap_group.add_argument(
        "-O",
        "--mark-overlap",
        dest="mark_overlap",
        action="store_true",
        default=False,
        help="Annotate nested and overlapping ORFs (default: %(default)s).",
    )
    overlap_group.add_argument(
        "-R",
        "--remove-discarded",
        dest="remove_discarded",
        action="store_true",
        default=False,
        help=("Remove discarded same-frame internal ORFs (default: %(default)s)."),
    )

    runtime_group = parser.add_argument_group("Runtime arguments")
    runtime_group.add_argument(
        "-t",
        "--thread",
        dest="threads",
        default=1,
        type=int,
        help="Number of worker processes (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments.

    Args:
        args: Parsed command-line arguments.

    Raises:
        ValueError: If a parameter value is invalid.
    """
    file_check(args.genome, args.annotation)

    if args.min_aa < 1:
        raise ValueError("--min-aa must be >= 1.")
    if args.max_aa < args.min_aa:
        raise ValueError("--max-aa must be >= --min-aa.")
    if args.kozak_up < 0 or args.kozak_down < 0:
        raise ValueError("Kozak-context lengths must be >= 0.")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1.")
    if not args.out_prefix:
        raise ValueError("--out-prefix must not be empty.")
    if not args.orf_prefix:
        raise ValueError("--orf-prefix must not be empty.")

    start_codons = [
        codon.strip().upper().replace("U", "T")
        for codon in args.start_codons.split(",")
        if codon.strip()
    ]
    if not start_codons:
        raise ValueError("--start-codons must contain at least one codon.")

    invalid_codons = [
        codon for codon in start_codons if len(codon) != 3 or set(codon).difference("ACGT")
    ]
    if invalid_codons:
        raise ValueError("Invalid start codon(s): " + ", ".join(invalid_codons))
    stop_as_start = sorted(set(start_codons).intersection({"TAA", "TAG", "TGA"}))
    if stop_as_start:
        raise ValueError("Stop codons cannot be used as start codons: " + ", ".join(stop_as_start))


def _parse_args(
    argv: Sequence[str] | None = None,
) -> Namespace:
    """Parse, validate, and print command-line arguments.

    Args:
        argv: Optional argument sequence used by tests or embedded callers.

    Returns:
        Validated command-line arguments.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        _validate_args(args)
    except ValueError as error:
        parser.error(str(error))

    args_print(args)
    return args


def _run_scanner_pipeline(args: Namespace) -> None:
    """Run the complete smORF scanning workflow.

    Args:
        args: Validated command-line arguments.
    """
    pipeline = SmORFPipeline(
        genome=args.genome,
        annotation=args.annotation,
        out_prefix=args.out_prefix,
        orf_prefix=args.orf_prefix,
        start_codons=args.start_codons,
        min_aa=args.min_aa,
        max_aa=args.max_aa,
        scan_strand=args.scan_strand,
        kozak_up=args.kozak_up,
        kozak_down=args.kozak_down,
        mark_overlap=args.mark_overlap,
        remove_discarded=args.remove_discarded,
        include_stop=args.include_stop,
        keep_partial=False,
        allow_ambiguous=False,
        threads=args.threads,
        retain_records=False,
    )
    pipeline.run()


def main(argv: Sequence[str] | None = None) -> None:
    """Run the smorf_scanner command.

    Args:
        argv: Optional argument sequence used by tests or embedded callers.
    """
    now_time()
    title_print("Scan transcript-centric small open reading frames.")
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)

    step_print(2, "Scan ORFs from transcript sequences.")
    _run_scanner_pipeline(args)

    title_print("All done.")
    now_time()


if __name__ == "__main__":
    main()
