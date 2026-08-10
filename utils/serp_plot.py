#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev003
# Function: Draw compact SeRP peak profiles for selected genes after peak calling.
# Input: serp_peak peak table, per-codon enrichment profile, and target gene/transcript identifiers.
# Output: Selected transcript-level figures with enrichment and core peaks.

"""Command-line entry point for targeted SeRP peak plotting."""

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
from utils.serp.Plot import SeRPPlot


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Draw SeRP enrichment profiles and called peak regions only for "
            "selected genes or transcripts. Run serp_peak first."
        ),
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        type=str,
        help=(
            "Input prefix previously used by serp_peak. The plotter reads "
            "<prefix>_peaks.log and <prefix>_peaks_ratio.txt."
        ),
    )

    target_group = parser.add_argument_group("Target genes")
    target_selection = target_group.add_mutually_exclusive_group(required=True)
    target_selection.add_argument(
        "-g",
        "--gene",
        "--target",
        dest="target",
        type=str,
        help="One transcript ID or gene name to plot.",
    )
    target_selection.add_argument(
        "--target-list",
        dest="target_list",
        type=str,
        help="Target list; the first non-empty column is used.",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        type=str,
        help="Optional output prefix; defaults to the serp_peak input prefix.",
    )
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["pdf", "png", "both"],
        default="pdf",
        help="Figure output format.",
    )

    display_group = parser.add_argument_group("Peak display arguments")
    display_group.add_argument(
        "--threshold",
        dest="threshold",
        default=None,
        type=float,
        help=(
            "Optional peak-threshold reference line. Set this to the same "
            "--enrich value used by serp_peak; it is not inferred automatically."
        ),
    )
    display_group.add_argument(
        "--show-max-site",
        dest="show_max_site",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Mark the reported maximum-enrichment site of each called peak.",
    )
    display_group.add_argument(
        "--shade-utr",
        dest="shade_utr",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use a subtle background to distinguish 5-prime and 3-prime UTRs.",
    )

    figure_group = parser.add_argument_group("Figure arguments")
    figure_group.add_argument(
        "--y-max",
        dest="y_max",
        default=None,
        type=float,
        help="Optional fixed y-axis maximum.",
    )
    figure_group.add_argument(
        "--font-size",
        dest="font_size",
        default=9.0,
        type=float,
        help="Base figure font size.",
    )
    figure_group.add_argument(
        "--dpi",
        dest="dpi",
        default=300,
        type=int,
        help="PNG output resolution.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate plotting arguments."""
    peak_file = args.input + "_peaks.log"
    profile_file = args.input + "_peaks_ratio.txt"
    file_check(peak_file, profile_file)
    if args.target_list:
        file_check(args.target_list)
    if args.y_max is not None and args.y_max <= 0:
        raise ValueError("--y-max must be > 0 when provided.")
    if args.font_size <= 0:
        raise ValueError("--font-size must be > 0.")
    if args.dpi < 1:
        raise ValueError("--dpi must be >= 1.")
    if args.threshold is not None and args.threshold < 0:
        args.threshold = None


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for ``serp_plot``."""
    now_time()
    title_print("Draw selected selective ribosome profiling peak profiles.")
    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)

    workflow = SeRPPlot(args)
    step_print(2, "Read target genes and called peak regions.")
    workflow.read_targets()
    workflow.read_peak_table()

    step_print(3, "Retrieve target enrichment profiles.")
    workflow.read_profiles()

    step_print(4, "Draw selected SeRP enrichment and peak-region figures.")
    workflow.draw_targets()

    complete_print()
    now_time()


if __name__ == "__main__":
    main()