#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.007
# Function: Calculate and visualize sample correlation from RPF density files.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Gene-level and codon-level correlation tables, heatmaps, and summary JSON.

"""Command-line entry point for RPF sample correlation analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Corr

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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Calculate sample correlations from RPF density files in JSONL or TXT format.",
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    input_group = parser.add_argument_group("Input and density arguments")
    input_group.add_argument(
        "-t",
        dest="transcript",
        required=False,
        default=None,
        type=str,
        help="Optional transcript filter table. If provided, transcript_id is used when available.",
    )
    input_group.add_argument(
        "--level",
        dest="level",
        choices=["gene", "rpf", "both"],
        default="both",
        type=str,
        help="Correlation level to calculate. Default: %(default)s.",
    )
    input_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        default=False,
        help="Normalize RPF counts to RPM before correlation. Default: %(default)s.",
    )
    input_group.add_argument(
        "--method",
        dest="method",
        choices=["pearson", "spearman", "kendall"],
        default="pearson",
        type=str,
        help="Correlation method. Default: %(default)s.",
    )

    filter_group = parser.add_argument_group("Filtering arguments")
    filter_group.add_argument(
        "--region",
        dest="region",
        choices=Corr.REGION_CHOICES,
        default="cds",
        type=str,
        help="Transcript region used for both gene-level and codon-level correlation. Default: %(default)s.",
    )
    filter_group.add_argument(
        "-m",
        "--min-gene-count",
        dest="min_gene_count",
        default=0.0,
        type=float,
        help="Retain genes with total raw RPF count >= this value. Default: %(default)s.",
    )

    plot_group = parser.add_argument_group("Plot arguments")
    plot_group.add_argument(
        "--cmap",
        dest="cmap",
        default="Blues",
        type=str,
        help="Heatmap colormap. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--figure-width",
        dest="figure_width",
        default=None,
        type=float,
        help="Optional figure width in inches.",
    )
    plot_group.add_argument(
        "--figure-height",
        dest="figure_height",
        default=None,
        type=float,
        help="Optional figure height in inches.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)

    if args.min_gene_count < 0:
        raise ValueError("--min-gene-count must be >= 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_corr_pipeline(args: Namespace) -> None:
    """Run the RPF correlation workflow."""
    corr = Corr.RPFCorrelation(args)

    step_print(2, "Import the RPF density file.")
    corr.import_rpf()

    step_print(3, "Build frame-specific density matrices.")
    corr.build_frame_tables()

    step_print(4, "Calculate sample correlations.")
    corr.calculate_correlations()

    step_print(5, "Output correlation tables.")
    corr.output_tables()

    step_print(6, "Draw correlation heatmaps.")
    if corr.gene_corr:
        corr.draw_corr_plot("gene", corr.gene_corr)
    if corr.rpf_corr:
        corr.draw_corr_plot("rpf", corr.rpf_corr)

    step_print(7, "Write summary JSON.")
    corr.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Corr."""
    now_time()
    title_print('Draw the correlation of samples.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_corr_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
