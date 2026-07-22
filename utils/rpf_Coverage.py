#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.013
# Function: Draw normalized metagene coverage profiles from RPF density files.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Metagene coverage tables, outlier tables, line plots, combined heatmap, and summary JSON.

"""Command-line entry point for RPF metagene coverage analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Coverage

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
        description="Draw normalized 5'UTR-CDS-3'UTR metagene coverage profiles from JSONL or TXT density files."
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

    parser.add_argument(
        "-t",
        dest="transcript",
        required=False,
        type=str,
        default=None,
        help="Optional transcript filter table in TXT format. If provided, transcript_id is used when available.",
    )
    parser.add_argument(
        "-f",
        dest="frame",
        choices=["0", "1", "2", "all"],
        required=False,
        type=str,
        default="all",
        help="Reading frame used for coverage calculation. Default: %(default)s.",
    )
    parser.add_argument(
        "--site",
        dest="site",
        choices=["E", "P", "A", "all"],
        required=False,
        type=str,
        default="P",
        help="Ribosomal site used for codon-level density shifting. Default: %(default)s.",
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=50,
        help="Retain transcripts with at least this sample-specific CDS RPF count. Default: %(default)s.",
    )
    parser.add_argument(
        "-b",
        "--bin",
        dest="bin",
        required=False,
        type=str,
        default="30,100,30",
        help="Normalize 5'UTR, CDS, and 3'UTR to these bin numbers. Default: %(default)s.",
    )
    parser.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        required=False,
        default=False,
        help="Normalize RPF counts to RPM before coverage aggregation. Default: %(default)s.",
    )
    parser.add_argument(
        "--set",
        dest="set",
        choices=["intersect", "union"],
        required=False,
        type=str,
        default="union",
        help=(
            "Transcript-region filtering strategy. 'union' treats missing region bins as zero; "
            "'intersect' keeps only transcripts with all requested regions. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--mode",
        dest="mode",
        required=False,
        choices=["line", "heatmap", "both", "all"],
        type=str,
        default="both",
        help="Coverage plot type. 'both' draws line plots and the combined sample heatmap. Default: %(default)s.",
    )
    parser.add_argument(
        "--plot-transform",
        dest="plot_transform",
        required=False,
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        type=str,
        default="none",
        help=(
            "Transform plotted density values without changing output tables. "
            "'log' is an alias of log1p. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--scale",
        dest="scale",
        required=False,
        choices=["none", "row"],
        type=str,
        default="none",
        help=(
            "Heatmap scaling method. 'row' applies row-wise min-max scaling "
            "to emphasize each sample profile shape. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--thread",
        dest="thread",
        type=int,
        required=False,
        default=1,
        help="Reserved for compatibility. Current coverage aggregation is vectorized. Default: %(default)s.",
    )

    output_group = parser.add_argument_group("Optional plot outputs")
    output_group.add_argument(
        "--bar",
        dest="barplot",
        action="store_true",
        required=False,
        default=False,
        help="Draw covered-transcript percentage barplots. Default: %(default)s.",
    )
    output_group.add_argument(
        "--gene-heatmap",
        "--heat",
        dest="gene_heatmap",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Deprecated compatibility option. Transcript-by-position heatmaps are no longer drawn. "
            "The combined sample heatmap is controlled by --mode heatmap/both/all."
        ),
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier",
        "--outlier",
        dest="remove_outlier",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Remove extreme transcript-bin RPF pileups before aggregation. "
            "This is useful for rRNA-like or mis-mapped fragments. Default: %(default)s."
        ),
    )
    outlier_group.add_argument(
        "--outlier-iqr",
        dest="outlier_iqr",
        type=float,
        required=False,
        default=8.0,
        help="IQR multiplier for extreme pileup detection. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-window",
        dest="outlier_window",
        type=int,
        required=False,
        default=5,
        help="Number of neighboring bins on each side used to estimate local background. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-local-fold",
        dest="outlier_local_fold",
        type=float,
        required=False,
        default=10.0,
        help="Minimum fold over local background for dynamic isolated-pileup filtering. Default: %(default)s.",
    )

    return parser


def _parse_bins(bin_arg: str) -> tuple[int, int, int]:
    """Parse and validate normalized metagene bin settings."""
    values = [item.strip() for item in str(bin_arg).split(",")]
    if len(values) != 3:
        raise ValueError("-b/--bin must contain three comma-separated integers, such as 30,100,30.")

    utr5_bin, cds_bin, utr3_bin = [int(value) for value in values]
    if utr5_bin < 0:
        raise ValueError("5'UTR bin number must be >= 0.")
    if cds_bin <= 0:
        raise ValueError("CDS bin number must be > 0.")
    if utr3_bin < 0:
        raise ValueError("3'UTR bin number must be >= 0.")
    return utr5_bin, cds_bin, utr3_bin


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)

    _parse_bins(args.bin)

    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.thread < 1:
        args.thread = 1
    if args.outlier_iqr <= 0:
        raise ValueError("--outlier-iqr must be > 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 1:
        raise ValueError("--outlier-local-fold must be > 1.")

    if getattr(args, "gene_heatmap", False):
        print(
            "Warning: --gene-heatmap/--heat is deprecated and will be ignored. "
            "Use --mode heatmap/both/all for the combined sample heatmap.",
            flush=True,
        )
        args.gene_heatmap = False


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_coverage_pipeline(args: Namespace) -> None:
    """Run the metagene coverage workflow."""
    coverage = Coverage.Coverage(args)

    step_print(2, "Import the RPFs file.")
    coverage.import_rpf()

    step_print(3, "Calculate metagene coverage profiles.")
    coverage.calc_coverage()

    step_print(4, "Output coverage tables.")
    coverage.output_coverage()

    step_print(5, "Draw coverage figures.")
    if args.mode in {"line", "both", "all"}:
        coverage.draw_meta_gene_line()
        coverage.draw_combined_line()
    if args.mode in {"heatmap", "both", "all"}:
        coverage.draw_heatmap()
    if args.barplot or args.mode == "all":
        coverage.draw_meta_gene_bar()

    step_print(6, "Write coverage summary.")
    coverage.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Coverage."""
    now_time()
    title_print('Draw the metagene coverage.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_coverage_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
