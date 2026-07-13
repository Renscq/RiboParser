#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.003
# Function: Calculate transcript-level CDS coefficient of variation from RPF density.
# Input: RPF density file in JSONL or TXT format, optional transcript list, and optional sample-group table.
# Output: Gene-level CoV tables, outlier records, group statistics, fitted mean-CoV curves, and summary JSON.

"""Command-line entry point for RPF coefficient-of-variation analysis."""

from __future__ import annotations

import argparse
import textwrap
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Coefficient_of_Variation
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Calculate transcript-level CDS coefficient of variation from JSONL or TXT RPF density.",
        epilog=textwrap.dedent(
            """\
            Group table format:
            Name\tGroup
            WT_1\tWT
            WT_2\tWT
            Treat_1\tTreat
            Treat_2\tTreat
            """
        ),
    )

    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-r", dest="rpf", required=True, type=str,
        help="Input RPF density file in JSONL or TXT format.",
    )
    required.add_argument(
        "-o", dest="output", required=True, type=str,
        help="Output prefix.",
    )

    input_group = parser.add_argument_group("Input filtering arguments")
    input_group.add_argument(
        "-g", "--group", dest="group", default=None, type=str,
        help="Optional sample-group table containing Name and Group columns.",
    )
    input_group.add_argument(
        "-l", "--list", dest="list", default=None, type=str,
        help="Optional transcript filter table. The transcript_id column or first column is used.",
    )
    input_group.add_argument(
        "-s", "--site", dest="site", choices=["E", "P", "A"], default="P",
        help="Ribosomal site used for positional density. Default: %(default)s.",
    )
    input_group.add_argument(
        "-f", "--frame", dest="frame", choices=["0", "1", "2", "all"], default="all",
        help="Reading frame used for CoV calculation. Default: %(default)s.",
    )
    input_group.add_argument(
        "-m", "--min", dest="min", type=int, default=5,
        help="Minimum sample-specific CDS RPF count required for a transcript. Default: %(default)s.",
    )
    input_group.add_argument(
        "--tis", dest="tis", type=int, default=15,
        help="Discard this many codons after the start codon. Default: %(default)s.",
    )
    input_group.add_argument(
        "--tts", dest="tts", type=int, default=5,
        help="Discard this many codons before the stop codon. Default: %(default)s.",
    )

    calculation = parser.add_argument_group("CoV calculation arguments")
    calculation.add_argument(
        "-n", "--normal", dest="normal", action="store_true", default=False,
        help="Convert each sample to RPM before reporting sums and means. CoV itself is scale invariant. Default: %(default)s.",
    )
    calculation.add_argument(
        "--ddof", dest="ddof", type=int, choices=[0, 1], default=1,
        help="Delta degrees of freedom used for positional SD. Use 1 for sample SD and 0 for population SD. Default: %(default)s.",
    )
    calculation.add_argument(
        "--min-codons", dest="min_codons", type=int, default=10,
        help="Minimum number of non-outlier CDS codon positions required for CoV. Default: %(default)s.",
    )
    calculation.add_argument(
        "--thread", dest="thread", type=int, default=1,
        help="Number of sample-level worker threads. Default: %(default)s.",
    )

    outlier = parser.add_argument_group("Outlier arguments")
    outlier.add_argument(
        "--remove-outlier", dest="remove_outlier", action="store_true", default=False,
        help="Remove isolated extreme RPF pileups before CoV calculation. Default: %(default)s.",
    )
    outlier.add_argument(
        "--outlier-iqr", dest="outlier_iqr", type=float, default=8.0,
        help="Robust log1p cutoff multiplier for candidate pileups. Default: %(default)s.",
    )
    outlier.add_argument(
        "--outlier-window", dest="outlier_window", type=int, default=5,
        help="Neighboring codons on each side used for local background. Default: %(default)s.",
    )
    outlier.add_argument(
        "--outlier-local-fold", dest="outlier_local_fold", type=float, default=10.0,
        help="Minimum fold over local background required for removal. Default: %(default)s.",
    )

    plot_group = parser.add_argument_group("Curve fitting and plotting arguments")
    plot_group.add_argument(
        "--fit-model", dest="fit_model", choices=["nb"], default="nb",
        help="Mean-CoV model. 'nb' fits CV = sqrt(alpha + beta / mean). Default: %(default)s.",
    )
    plot_group.add_argument(
        "--fit-quantile", dest="fit_quantile", type=float, default=0.01,
        help="Symmetric tail fraction excluded before fitting; 0 disables trimming. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--plot-transform", dest="plot_transform", choices=["log2"], default="log2",
        help="Axis transformation for the mean-CoV scatter. Default: %(default)s.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.group:
        file_check(args.group)
    if args.list:
        file_check(args.list)
    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if args.min_codons < 2:
        raise ValueError("--min-codons must be >= 2.")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.outlier_iqr <= 0:
        raise ValueError("--outlier-iqr must be > 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 1:
        raise ValueError("--outlier-local-fold must be > 1.")
    if not 0 <= args.fit_quantile < 0.5:
        raise ValueError("--fit-quantile must be in [0, 0.5).")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print a standardized workflow step."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_pipeline(args: Namespace) -> None:
    """Run the complete coefficient-of-variation workflow."""
    cov = Coefficient_of_Variation.CoV(args)

    _print_step(2, "Import the RPF density file.")
    cov.import_rpf()

    _print_step(3, "Calculate sample-specific transcript CoV.")
    cov.calculate_cov()

    _print_step(4, "Read sample groups and compare CoV distributions.")
    cov.read_group()
    cov.compare_groups()

    _print_step(5, "Fit the corrected mean-CoV relationship.")
    cov.fit_mean_cov()

    _print_step(6, "Output CoV tables.")
    cov.output_tables()

    _print_step(7, "Draw CoV figures.")
    cov.draw_fit_plot()
    cov.draw_distribution_plot()

    _print_step(8, "Write the analysis summary.")
    cov.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_CoV."""
    now_time()
    print("\nCalculate positional coefficient of variation in CDS regions.", flush=True)
    _print_step(1, "Checking the input arguments.")
    args = _parse_args(argv)
    _run_pipeline(args)
    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
