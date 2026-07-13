#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.005
# Function: Detect transcript-level ribosomal frameshift candidates.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Frame-transition statistics, candidate tables, summaries, and figures.

"""Command-line entry point for ribosomal frameshift detection."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Shift
from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Scan transcript CDS regions for statistically supported changes in "
            "the three-frame RPF distribution."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r", dest="rpf", required=True, type=str,
        help="Input frame-resolved RPF density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-o", dest="output", required=True, type=str,
        help="Output prefix.",
    )

    input_group = parser.add_argument_group("Input filtering arguments")
    input_group.add_argument(
        "-t", "-l", "--transcript", dest="transcript", required=False,
        type=str, default=None,
        help=(
            "Optional transcript filter table. The transcript_id column is used "
            "when available; otherwise the first column is used."
        ),
    )
    input_group.add_argument(
        "-m", dest="min", required=False, type=int, default=50,
        help=(
            "Minimum sample-specific CDS RPF count required for scanning a "
            "transcript. Default: %(default)s."
        ),
    )
    input_group.add_argument(
        "--tis", dest="tis", required=False, type=int, default=10,
        help="Discard this number of CDS codons after the TIS. Default: %(default)s AA.",
    )
    input_group.add_argument(
        "--tts", dest="tts", required=False, type=int, default=5,
        help="Discard this number of CDS codons before the TTS. Default: %(default)s AA.",
    )

    scan_group = parser.add_argument_group("Frameshift scan arguments")
    scan_group.add_argument(
        "-s", dest="site", choices=["E", "P", "A"], default="P",
        help="Ribosomal site used for frame-density analysis. Default: %(default)s.",
    )
    scan_group.add_argument(
        "-p", "--period", dest="period", type=float, default=45.0,
        help=(
            "Minimum downstream percentage in the shifted frame. Range: 0-100. "
            "Default: %(default)s."
        ),
    )
    scan_group.add_argument(
        "--pre-period", dest="pre_period", type=float, default=45.0,
        help=(
            "Minimum upstream percentage in frame 0 before the candidate shift. "
            "Range: 0-100. Default: %(default)s."
        ),
    )
    scan_group.add_argument(
        "--min-shift", dest="min_shift", type=float, default=0.15,
        help=(
            "Minimum increase in the target shifted-frame proportion from upstream "
            "to downstream. Range: 0-1. Default: %(default)s."
        ),
    )
    scan_group.add_argument(
        "--min-segment", dest="min_segment", type=int, default=20,
        help=(
            "Minimum codon positions required on each side of a scanned change "
            "point. Default: %(default)s."
        ),
    )
    scan_group.add_argument(
        "--min-segment-reads", dest="min_segment_reads", type=float, default=20.0,
        help=(
            "Minimum total RPF count required on each side of a change point. "
            "Default: %(default)s."
        ),
    )
    scan_group.add_argument(
        "--scan-step", dest="scan_step", type=int, default=1,
        help="Codon step between candidate change points. Default: %(default)s.",
    )
    scan_group.add_argument(
        "--alpha", dest="alpha", type=float, default=0.05,
        help=(
            "Maximum BH-FDR for a reported candidate after within-transcript "
            "Bonferroni correction. Default: %(default)s."
        ),
    )
    scan_group.add_argument(
        "--thread", dest="thread", type=int, default=1,
        help="Number of sample-level worker threads. Default: %(default)s.",
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier", dest="remove_outlier", action="store_true", default=False,
        help=(
            "Remove isolated extreme total-density pileups before frame-shift "
            "scanning. Default: %(default)s."
        ),
    )
    outlier_group.add_argument(
        "--outlier-iqr", dest="outlier_iqr", type=float, default=8.0,
        help="IQR multiplier for global pileup detection. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-window", dest="outlier_window", type=int, default=5,
        help="Neighboring codons on each side for local confirmation. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-local-fold", dest="outlier_local_fold", type=float, default=10.0,
        help="Minimum fold over local background for removal. Default: %(default)s.",
    )

    plot_group = parser.add_argument_group("Output and plotting arguments")
    plot_group.add_argument(
        "--smooth-window", dest="smooth_window", type=int, default=5,
        help="Centered codon window used only for candidate profile plots. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--plot-bin-size", dest="plot_bin_size", type=int, default=5,
        help=(
            "Codon bin size used for stacked frame-composition bar plots in candidate "
            "figures. Default: %(default)s."
        ),
    )
    plot_group.add_argument(
        "--gene-plot-mode", dest="gene_plot_mode",
        choices=["bar", "heatmap", "both"], default="both",
        help=(
            "Candidate plot content selector retained for compatibility; the figure now uses a fixed "
            "three-row layout containing density, stacked frame proportions, and heatmap. "
            "Default: %(default)s."
        ),
    )
    plot_group.add_argument(
        "--gene-fig", dest="gene_fig", choices=["none", "png", "pdf", "both"],
        default="png",
        help="Candidate transcript figure format. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--max-gene-figures", dest="max_gene_figures", type=int, default=500,
        help=(
            "Maximum candidate figures to create, ranked by FDR and effect. Use 0 "
            "for all candidates. Default: %(default)s."
        ),
    )
    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)
    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
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
    """Parse and validate command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    return args


def _print_step(number: int, message: str) -> None:
    """Print a standardized pipeline step."""
    print(f"\nStep{number}: {message}", flush=True)


def _run_pipeline(args: Namespace) -> None:
    """Run the complete frameshift detection workflow."""
    shift = Shift.Shift(args)

    _print_step(2, "Import the frame-resolved RPF density file.")
    shift.import_rpf()

    _print_step(3, "Scan transcript change points and calculate p-values.")
    shift.scan_frame_shift()

    _print_step(4, "Output frameshift tables and summary information.")
    shift.output_results()

    _print_step(5, "Draw summary and candidate-level frameshift figures.")
    shift.draw_all()


def main(argv: Sequence[str] | None = None) -> None:
    """Run the command-line program."""
    now_time()
    print("\nDetect transcript-level ribosomal frameshift candidates.", flush=True)
    _print_step(1, "Check input arguments.")
    args = _parse_args(argv)
    args_print(args)
    _run_pipeline(args)
    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
