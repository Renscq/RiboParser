#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Detect transcript-level ribosomal frameshift candidates.
# Input: Frame-resolved RPF density file and optional transcript filter.
# Output: Frame-transition statistics, candidate tables, summaries, and figures.

"""Command-line entry point for ribosomal frameshift detection."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Shift

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
            "Scan transcript CDS regions for statistically supported changes "
            "in the three-frame RPF distribution."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input frame-resolved RPF density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    input_group = parser.add_argument_group("Filtering arguments")
    input_group.add_argument(
        "-t",
        "--transcript",
        dest="transcript",
        default=None,
        type=str,
        help="Optional transcript filter table.",
    )
    input_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=50,
        type=int,
        help="Minimum sample-specific CDS RPF count.",
    )
    input_group.add_argument(
        "--tis",
        dest="tis",
        default=10,
        type=int,
        help="CDS codons discarded after TIS.",
    )
    input_group.add_argument(
        "--tts",
        dest="tts",
        default=5,
        type=int,
        help="CDS codons discarded before TTS.",
    )

    scan_group = parser.add_argument_group("Frameshift scan arguments")
    scan_group.add_argument(
        "-s",
        "--site",
        dest="site",
        choices=["E", "P", "A"],
        default="P",
        help="Ribosomal site used for frame-density analysis.",
    )
    scan_group.add_argument(
        "-p",
        "--period",
        dest="period",
        default=45.0,
        type=float,
        help="Minimum downstream percentage in the shifted frame.",
    )
    scan_group.add_argument(
        "--pre-period",
        dest="pre_period",
        default=45.0,
        type=float,
        help="Minimum upstream percentage in frame 0.",
    )
    scan_group.add_argument(
        "--min-shift",
        dest="min_shift",
        default=0.15,
        type=float,
        help="Minimum increase in the shifted-frame proportion.",
    )
    scan_group.add_argument(
        "--min-segment",
        dest="min_segment",
        default=20,
        type=int,
        help="Minimum codon positions required on each side of a change point.",
    )
    scan_group.add_argument(
        "--min-segment-reads",
        dest="min_segment_reads",
        default=20.0,
        type=float,
        help="Minimum total RPF count required on each side of a change point.",
    )
    scan_group.add_argument(
        "--scan-step",
        dest="scan_step",
        default=1,
        type=int,
        help="Codon step between candidate change points.",
    )
    scan_group.add_argument(
        "--alpha",
        dest="alpha",
        default=0.05,
        type=float,
        help="Maximum adjusted significance value for reported candidates.",
    )
    scan_group.add_argument(
        "--thread",
        dest="thread",
        default=1,
        type=int,
        help="Number of sample-level worker threads.",
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier",
        dest="remove_outlier",
        action="store_true",
        default=False,
        help="Remove isolated extreme total-density pileups.",
    )
    outlier_group.add_argument(
        "--outlier-iqr",
        dest="outlier_iqr",
        default=8.0,
        type=float,
        help="IQR multiplier for global pileup detection.",
    )
    outlier_group.add_argument(
        "--outlier-window",
        dest="outlier_window",
        default=5,
        type=int,
        help="Neighboring codons used for local outlier confirmation.",
    )
    outlier_group.add_argument(
        "--outlier-local-fold",
        dest="outlier_local_fold",
        default=10.0,
        type=float,
        help="Minimum fold over local background for outlier removal.",
    )

    plotting_group = parser.add_argument_group("Plotting arguments")
    plotting_group.add_argument(
        "--smooth-window",
        dest="smooth_window",
        default=5,
        type=int,
        help="Centered codon window used for candidate profile plots.",
    )
    plotting_group.add_argument(
        "--plot-bin-size",
        dest="plot_bin_size",
        default=5,
        type=int,
        help="Codon bin size used for stacked frame-composition plots.",
    )
    plotting_group.add_argument(
        "--gene-plot-mode",
        dest="gene_plot_mode",
        choices=["bar", "heatmap", "both"],
        default="both",
        help="Candidate plot content selector retained for compatibility.",
    )
    plotting_group.add_argument(
        "--gene-fig",
        dest="gene_fig",
        choices=["none", "png", "pdf", "both"],
        default="png",
        help="Candidate transcript figure format.",
    )
    plotting_group.add_argument(
        "--max-gene-figures",
        dest="max_gene_figures",
        default=500,
        type=int,
        help="Maximum number of candidate figures; use 0 for all.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)
    if args.min < 0:
        raise ValueError("-m must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if not 0 <= args.period <= 100:
        raise ValueError("--period must be between 0 and 100.")
    if not 0 <= args.pre_period <= 100:
        raise ValueError("--pre-period must be between 0 and 100.")
    if not 0 <= args.min_shift <= 1:
        raise ValueError("--min-shift must be between 0 and 1.")
    if args.min_segment < 2:
        raise ValueError("--min-segment must be >= 2.")
    if args.min_segment_reads < 0:
        raise ValueError("--min-segment-reads must be >= 0.")
    if args.scan_step < 1:
        raise ValueError("--scan-step must be >= 1.")
    if not 0 < args.alpha <= 1:
        raise ValueError("--alpha must be in (0, 1].")
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.outlier_iqr < 0:
        raise ValueError("--outlier-iqr must be >= 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 0:
        raise ValueError("--outlier-local-fold must be > 0.")
    if args.smooth_window < 1:
        raise ValueError("--smooth-window must be >= 1.")
    if args.plot_bin_size < 1:
        raise ValueError("--plot-bin-size must be >= 1.")
    if args.max_gene_figures < 0:
        raise ValueError("--max-gene-figures must be >= 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_shift_pipeline(args: Namespace) -> None:
    """Run the complete frameshift-detection workflow."""
    shift = Shift.Shift(args)

    step_print(2, "Import the frame-resolved RPF density file.")
    shift.import_rpf()

    step_print(3, "Scan transcript change points and calculate p-values.")
    shift.scan_frame_shift()

    step_print(4, "Output frameshift tables and summary information.")
    shift.output_results()

    step_print(5, "Draw summary and candidate-level frameshift figures.")
    shift.draw_all()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Shift."""
    now_time()
    title_print('Detect transcript-level ribosomal frameshift candidates.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_shift_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
