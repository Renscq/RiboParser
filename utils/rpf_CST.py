#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.001
# Function: Calculate sample-specific codon selection time from paired RPF and RNA density data.
# Input: RPF and RNA density files in JSONL or TXT format and an optional transcript filter.
# Output: CST tables, outlier records, summary JSON, correlation, heatmap, rank, and convergence plots.

"""Command-line entry point for codon selection time analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import CST

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
            "Calculate sample-specific codon selection time from paired RPF "
            "and RNA density files in JSONL or TXT format."
        )
    )

    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "--rpf", 
        dest="rpf",
        required=True, type=str,
        help="Input RPF density file in JSONL or TXT format.",
    )
    required.add_argument(
        "--rna", 
        dest="rna",
        required=True, type=str,
        help="Input RNA density file in JSONL or TXT format.",
    )
    required.add_argument(
        "-o", 
        "--output",
        dest="output", required=True, type=str,
        help="Output prefix.",
    )

    filtering = parser.add_argument_group("Filtering arguments")
    filtering.add_argument(
        "-l", 
        "--list", 
        dest="list", default=None, type=str,
        help=(
            "Optional transcript filter table. The transcript_id column is used "
            "when available; otherwise the first column is used."
        ),
    )
    filtering.add_argument(
        "--min",
        dest="min", default=30, type=int,
        help="Minimum sample-specific CDS RPF count per transcript. Default: %(default)s.",
    )
    filtering.add_argument(
        "--min-rna", 
        dest="min_rna", default=1, type=int,
        help="Minimum sample-specific CDS RNA count per transcript. Default: %(default)s.",
    )
    filtering.add_argument(
        "--tis", default=0, type=int,
        help="Discard this number of CDS codons after TIS. Default: %(default)s AA.",
    )
    filtering.add_argument(
        "--tts", default=0, type=int,
        help="Discard this number of CDS codons before TTS. Default: %(default)s AA.",
    )

    calculation = parser.add_argument_group("CST calculation arguments")
    calculation.add_argument(
        "-s", "--site", dest="site", choices=["E", "P", "A"], default="P",
        help="Ribosomal site used for RPF CST calculation. Default: %(default)s.",
    )
    calculation.add_argument(
        "-f", "--frame", dest="frame", choices=["0", "1", "2", "all"], default="all",
        help="Reading frame used for RPF and RNA density import. Default: %(default)s.",
    )
    calculation.add_argument(
        "-t", "--times", dest="times", type=int, default=10,
        help="Maximum CST iteration count. Default: %(default)s.",
    )
    calculation.add_argument(
        "--tolerance", type=float, default=1e-6,
        help="Stop iterations when maximum absolute CST change is below this value. Default: %(default)s.",
    )
    calculation.add_argument(
        "--thread", type=int, default=1,
        help="Number of sample-pair worker threads. Default: %(default)s.",
    )

    outlier = parser.add_argument_group("Outlier arguments")
    outlier.add_argument(
        "--remove-outlier", action="store_true", default=False,
        help="Remove isolated extreme RPF pileups before CST calculation. Default: %(default)s.",
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

    output = parser.add_argument_group("Plotting arguments")
    output.add_argument(
        "--scale", choices=["none", "minmax", "zscore"], default="minmax",
        help="Sample-wise scaling for relative CST values. Default: %(default)s.",
    )
    output.add_argument(
        "--plot-transform",
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        default="none",
        help="Transform absolute CST values for plotting only. Default: %(default)s.",
    )
    output.add_argument(
        "--rankplot-ncol",
        dest="rankplot_ncol", type=int, default=1,
        help=(
            "Number of panels per row in the rank plot. Set to 2 for a more "
            "compact two-column layout. Default: %(default)s."
        ),
    )
    output.add_argument(
        "--all", action="store_true", default=False,
        help="Output detailed gene-level CST metrics. Default: %(default)s.",
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
    if args.times < 0:
        raise ValueError("-t/--times must be >= 0.")
    if args.tolerance <= 0:
        raise ValueError("--tolerance must be > 0.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.rankplot_ncol < 1:
        raise ValueError("--rankplot-ncol must be >= 1.")
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


def _run_pipeline(args: Namespace) -> None:
    """Run the codon selection time workflow."""
    cst = CST.CodonSelectiveTime(args)

    step_print(2, "Import paired RPF and RNA density files.")
    cst.import_density()

    step_print(3, "Calculate sample-specific codon selection time.")
    cst.calculate_cst()

    step_print(4, "Output codon selection time tables.")
    cst.output_tables()

    step_print(5, "Draw codon selection time figures.")
    cst.draw_cst_corr()
    cst.draw_cst_heat()
    cst.draw_cst_rank()
    cst.draw_cst_convergence()

    step_print(6, "Write codon selection time summary.")
    cst.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_CST."""
    now_time()
    title_print('Calculate codon selection time.')
    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)
    _run_pipeline(args)
    complete_print()
    now_time()


if __name__ == "__main__":
    main()
