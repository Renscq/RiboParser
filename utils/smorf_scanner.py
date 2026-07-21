#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Scan transcript-centric small open reading frames.
# Input: Genome FASTA and genePred transcript annotation.
# Output: Candidate smORF annotation, sequence, and summary files.

"""Command-line entry point for transcript-centric smORF scanning."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import args_print, file_check, now_time
from utils.smorf import SmORFPipeline


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
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
        help="Input genome sequence file in FASTA format.",
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
        help="Minimum ORF length (default: %(default)s aa).",
    )
    scanning_group.add_argument(
        "-M",
        "--max-aa",
        dest="max_aa",
        default=10000,
        type=int,
        help="Maximum ORF length (default: %(default)s aa).",
    )
    scanning_group.add_argument(
        "-x",
        "--scan-strand",
        dest="scan_strand",
        choices=["sense", "antisense", "both"],
        default="sense",
        type=str,
        help="Transcript strand scanned for ORFs (default: %(default)s).",
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
        help="Downstream Kozak-context length (default: %(default)s nt).",
    )
    scanning_group.add_argument(
        "-O",
        "--mark-overlap",
        dest="mark_overlap",
        action="store_true",
        default=False,
        help="Annotate nested and overlapping ORFs (default: %(default)s).",
    )
    scanning_group.add_argument(
        "-R",
        "--remove-discarded",
        dest="remove_discarded",
        action="store_true",
        default=False,
        help="Remove discarded same-frame internal ORFs (default: %(default)s).",
    )
    scanning_group.add_argument(
        "-I",
        "--include-stop",
        dest="include_stop",
        action="store_true",
        default=False,
        help="Retain the stop-codon symbol in peptide sequences (default: %(default)s).",
    )

    runtime_group = parser.add_argument_group("Runtime arguments")
    runtime_group.add_argument(
        "-t",
        "--threads",
        dest="threads",
        default=1,
        type=int,
        help="Number of worker processes (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.genome, args.annotation)
    if args.min_aa < 1:
        raise ValueError("--min-aa must be >= 1.")
    if args.max_aa < args.min_aa:
        raise ValueError("--max-aa must be >= --min-aa.")
    if args.kozak_up < 0 or args.kozak_down < 0:
        raise ValueError("Kozak-context lengths must be >= 0.")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1.")

    start_codons = [
        codon.strip().upper()
        for codon in args.start_codons.split(",")
        if codon.strip()
    ]
    if not start_codons:
        raise ValueError("--start-codons must contain at least one codon.")
    invalid_codons = [
        codon
        for codon in start_codons
        if len(codon) != 3 or set(codon).difference("ACGT")
    ]
    if invalid_codons:
        raise ValueError(
            "Invalid start codon(s): " + ", ".join(invalid_codons)
        )


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


def _run_scanner_pipeline(args: Namespace) -> None:
    """Run transcript-centric smORF scanning."""
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
        threads=args.threads,
    )

    _print_step(2, "Scan ORFs from transcript sequences.")
    pipeline.run()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for smorf_scanner."""
    now_time()
    print("\nScan transcript-centric small open reading frames.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_scanner_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
