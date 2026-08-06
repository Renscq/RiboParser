#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-06
# Version: dev004
# Function: Visualize all detected ribosome pausing sites for selected genes.
# Input: rpf_Odd_Ratio site table and compact RPF density JSONL file.
# Output: Gene-level RPF pausing-site figures and optional metrics table.

"""Command-line entry point for ribosome pausing-site visualization."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import PSplot
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
            "Draw all pausing sites from selected genes on sample-resolved "
            "raw RPF profiles."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        type=str,
        help="Input rpf_Odd_Ratio site-level TXT table.",
    )
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input compact RPF density JSONL or JSONL.GZ file.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    target_group = parser.add_argument_group("Target genes")
    target_selection = target_group.add_mutually_exclusive_group(required=True)
    target_selection.add_argument(
        "-g",
        "--gene",
        "--target",
        dest="target",
        type=str,
        help="Gene ID, transcript ID, or result-table name to plot.",
    )
    target_selection.add_argument(
        "--target-list",
        dest="target_list",
        type=str,
        help="Target list; the first non-empty column is used.",
    )

    view_group = parser.add_argument_group("Plot view and scale")
    view_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        default=False,
        help="Plot RPM-normalized density; sample pause metrics still use raw counts.",
    )

    view_group.add_argument(
        "--view",
        dest="view",
        choices=["region", "gene"],
        default="gene",
        help="Pause-span regional or whole-gene view (default: %(default)s).",
    )
    view_group.add_argument(
        "--flank",
        dest="flank",
        default=30,
        type=int,
        help="Upstream/downstream codons in region view (default: %(default)s).",
    )
    view_group.add_argument(
        "--plot-transform",
        dest="plot_transform",
        choices=["none", "sqrt", "log1p", "log2", "log10"],
        default="none",
        help="Display-only density transformation (default: %(default)s).",
    )
    view_group.add_argument(
        "--y-scale",
        dest="y_scale",
        choices=["shared", "sample"],
        default="shared",
        help="Shared or per-sample y-axis range (default: %(default)s).",
    )
    view_group.add_argument(
        "--y-max",
        dest="y_max",
        default=None,
        type=float,
        help="Optional fixed y-axis maximum after transformation.",
    )
    view_group.add_argument(
        "--pause-region",
        dest="pause_region",
        action="store_true",
        default=False,
        help="Highlight each pause codon with a background region.",
    )

    output_group = parser.add_argument_group("Output and layout")
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["pdf", "png", "both"],
        default="pdf",
        help="Figure output format (default: %(default)s).",
    )
    output_group.add_argument(
        "--export-metrics",
        dest="export_metrics",
        action="store_true",
        default=False,
        help="Export per-event, per-sample pause metrics.",
    )
    output_group.add_argument(
        "--dpi",
        dest="dpi",
        default=300,
        type=int,
        help="PNG output resolution (default: %(default)s).",
    )
    output_group.add_argument(
        "--font-size",
        dest="font_size",
        default=9.0,
        type=float,
        help="Base figure font size (default: %(default)s).",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.input, args.rpf)
    if args.target_list:
        file_check(args.target_list)
    if args.flank < 1:
        raise ValueError("--flank must be >= 1.")
    if args.y_max is not None and args.y_max <= 0:
        raise ValueError("--y-max must be > 0 when provided.")
    if args.dpi < 1:
        raise ValueError("--dpi must be positive.")
    if args.font_size <= 0:
        raise ValueError("--font-size must be > 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_PSplot."""
    now_time()
    title_print("Draw sample-resolved ribosome pausing-site profiles.")
    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)

    workflow = PSplot.PSplot(args)
    step_print(2, "Import and filter pausing-site events.")
    workflow.read_events()
    step_print(3, "Stream the JSON density file and retrieve target transcripts.")
    workflow.read_rpf_records()
    step_print(4, "Draw pausing-site profiles across samples.")
    workflow.draw_events()
    if args.export_metrics:
        step_print(5, "Export sample-level pause metrics.")
        workflow.output_metrics()

    complete_print()
    now_time()


if __name__ == "__main__":
    main()