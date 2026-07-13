#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.001
# Function: Calculate sample-specific codon decoding time from paired RPF and RNA density data.
# Input: RPF and RNA density files in JSONL or TXT format and an optional transcript filter.
# Output: CDT tables, outlier records, summary JSON, correlation, heatmap, and rank plots.

"""Command-line entry point for codon decoding time analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import CDT
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate sample-specific codon decoding time from paired RPF and "
            "RNA density files in JSONL or TXT format."
        )
    )

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--rpf", required=True, type=str,
                          help="Input RPF density file in JSONL or TXT format.")
    required.add_argument("--rna", required=True, type=str,
                          help="Input RNA density file in JSONL or TXT format.")
    required.add_argument("-o", dest="output", required=True, type=str,
                          help="Output prefix.")

    filtering = parser.add_argument_group("Input filtering arguments")
    filtering.add_argument(
        "-l", "--list", dest="list", default=None, type=str,
        help=(
            "Optional transcript filter table. The transcript_id column is used "
            "when available; otherwise the first column is used."
        ),
    )
    filtering.add_argument(
        "-m", dest="min", default=30, type=int,
        help="Minimum sample-specific CDS RPF count per transcript. Default: %(default)s.",
    )
    filtering.add_argument(
        "--min-rna", dest="min_rna", default=1, type=int,
        help="Minimum sample-specific CDS RNA count per transcript. Default: %(default)s.",
    )
    filtering.add_argument(
        "--tis", default=15, type=int,
        help="Discard this number of CDS codons after TIS. Default: %(default)s AA.",
    )
    filtering.add_argument(
        "--tts", default=5, type=int,
        help="Discard this number of CDS codons before TTS. Default: %(default)s AA.",
    )

    calculation = parser.add_argument_group("CDT calculation arguments")
    calculation.add_argument(
        "-s", dest="site", choices=["E", "P", "A"], default="P",
        help="Ribosomal site used for RPF CDT calculation. Default: %(default)s.",
    )
    calculation.add_argument(
        "-f", dest="frame", choices=["0", "1", "2", "all"], default="all",
        help="Reading frame used for RPF and RNA density import. Default: %(default)s.",
    )
    calculation.add_argument(
        "--thread", type=int, default=1,
        help="Number of sample-pair worker threads. Default: %(default)s.",
    )

    outlier = parser.add_argument_group("Outlier arguments")
    outlier.add_argument(
        "--remove-outlier", action="store_true", default=False,
        help="Remove isolated extreme RPF pileups before CDT calculation. Default: %(default)s.",
    )
    outlier.add_argument(
        "--outlier-iqr", type=float, default=8.0,
        help="IQR multiplier for global extreme-pileup detection. Default: %(default)s.",
    )
    outlier.add_argument(
        "--outlier-window", type=int, default=5,
        help="Neighboring codons on each side used to confirm outliers. Default: %(default)s.",
    )
    outlier.add_argument(
        "--outlier-local-fold", type=float, default=10.0,
        help="Minimum fold over local background for outlier removal. Default: %(default)s.",
    )

    output = parser.add_argument_group("Output and plotting arguments")
    output.add_argument(
        "--scale", choices=["none", "minmax", "zscore"], default="minmax",
        help="Sample-wise scaling for relative CDT values. Default: %(default)s.",
    )
    output.add_argument(
        "--plot-transform",
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        default="none",
        help="Transform absolute CDT values for plotting only. Default: %(default)s.",
    )
    output.add_argument(
        "--all", action="store_true", default=False,
        help="Output detailed position-level and gene-codon CDT tables. Default: %(default)s.",
    )
    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf, args.rna)
    if args.list:
        file_check(args.list)
    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.min_rna < 0:
        raise ValueError("--min-rna must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.outlier_iqr <= 0:
        raise ValueError("--outlier-iqr must be > 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 1:
        raise ValueError("--outlier-local-fold must be > 1.")


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


def _run_pipeline(args: Namespace) -> None:
    """Run the codon decoding time workflow."""
    cdt = CDT.CodonDecodingTime(args)

    _print_step(2, "Import paired RPF and RNA density files.")
    cdt.import_density()

    _print_step(3, "Calculate sample-specific codon decoding time.")
    cdt.calculate_cdt()

    _print_step(4, "Output codon decoding time tables.")
    cdt.output_tables()

    _print_step(5, "Draw codon decoding time figures.")
    cdt.draw_cdt_corr()
    cdt.draw_cdt_heat()
    cdt.draw_cdt_rank()

    _print_step(6, "Write codon decoding time summary.")
    cdt.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_CDT."""
    now_time()
    print("\nCalculate codon decoding time.", flush=True)
    _print_step(1, "Checking the input arguments.")
    args = _parse_args(argv)
    _run_pipeline(args)
    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
