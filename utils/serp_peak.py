#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Detect selective ribosome profiling binding peaks.
# Input: SeRP RPF coverage table, sample definitions, and optional annotation.
# Output: Peak tables, optional ratios, RPM tables, and figures.

"""Command-line entry point for SeRP peak detection."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import args_print, file_check, now_time
from utils.serp.SeRP import SeRP


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Detect selective ribosome profiling binding peaks.",
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF coverage table.",
    )
    required_group.add_argument(
        "--ck",
        dest="control",
        required=True,
        type=str,
        help="Comma-separated control sample names.",
    )
    required_group.add_argument(
        "--ip",
        dest="ip",
        required=True,
        type=str,
        help="Comma-separated immunoprecipitation sample names.",
    )

    input_group = parser.add_argument_group("Input arguments")
    input_group.add_argument(
        "-n",
        dest="norm",
        default=None,
        type=str,
        help="Optional file containing sample-level RPF totals.",
    )
    input_group.add_argument(
        "--scale",
        dest="scale",
        default=1e6,
        type=float,
        help="Normalization scale.",
    )
    input_group.add_argument(
        "-a",
        dest="anno",
        default=None,
        type=str,
        help="Optional gene annotation table.",
    )

    filtering_group = parser.add_argument_group("Data filtering arguments")
    filtering_group.add_argument(
        "-m",
        dest="min",
        default=50,
        type=int,
        help="Minimum gene-level RPF coverage.",
    )
    filtering_group.add_argument(
        "--corr",
        dest="corr",
        default=0.3,
        type=float,
        help="Minimum replicate correlation.",
    )

    peak_group = parser.add_argument_group("Peak scanning arguments")
    peak_group.add_argument(
        "-f",
        dest="fill",
        default=30,
        type=int,
        help="Background fill mode or leading-codon window.",
    )
    peak_group.add_argument(
        "-s",
        dest="size",
        default=3,
        type=int,
        help="Savitzky-Golay smoothing window size.",
    )
    peak_group.add_argument(
        "-k",
        dest="k",
        default=1,
        type=int,
        help="Savitzky-Golay polynomial order.",
    )
    peak_group.add_argument(
        "-w",
        dest="width",
        default=5,
        type=int,
        help="Minimum binding-peak width in amino acids.",
    )
    peak_group.add_argument(
        "-e",
        dest="enrich",
        default=2.0,
        type=float,
        help="Peak-height enrichment threshold.",
    )
    peak_group.add_argument(
        "-c",
        dest="collision",
        default=1.5,
        type=float,
        help="Collision-region enrichment threshold.",
    )
    peak_group.add_argument(
        "-g",
        dest="gaps",
        default=1,
        type=int,
        help="Maximum consecutive gap length within a peak.",
    )
    peak_group.add_argument(
        "-p",
        dest="proportion",
        default=0.2,
        type=float,
        help="Maximum gap proportion within a peak.",
    )
    peak_group.add_argument(
        "--back",
        dest="background",
        default=0,
        type=int,
        help="Leading-codon background window.",
    )
    peak_group.add_argument(
        "--bf",
        dest="backFold",
        action="store_true",
        default=True,
        help="Apply background-fold filtering.",
    )
    peak_group.add_argument(
        "--up",
        dest="upstream",
        default=10,
        type=int,
        help="Upstream codons retained around each peak.",
    )
    peak_group.add_argument(
        "--down",
        dest="downstream",
        default=10,
        type=int,
        help="Downstream codons retained around each peak.",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-o",
        dest="output",
        default="results",
        type=str,
        help="Output file prefix.",
    )
    output_group.add_argument(
        "--all",
        dest="all",
        action="store_true",
        default=False,
        help="Retain all detected peak regions.",
    )
    output_group.add_argument(
        "--rpm",
        dest="rpm",
        action="store_true",
        default=False,
        help="Output peak-level RPM values.",
    )
    output_group.add_argument(
        "--ratio",
        dest="ratio",
        action="store_true",
        default=False,
        help="Output the original enrichment-ratio table.",
    )
    output_group.add_argument(
        "--fig",
        dest="fig",
        action="store_true",
        default=False,
        help="Draw peak-scanning demonstration figures.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.norm:
        file_check(args.norm)
    if args.anno:
        file_check(args.anno)

    if args.scale <= 0:
        raise ValueError("--scale must be > 0.")
    if args.min < 0:
        raise ValueError("-m must be >= 0.")
    if not -1 <= args.corr <= 1:
        raise ValueError("--corr must be in [-1, 1].")
    if args.size < 1 or args.size % 2 == 0:
        raise ValueError("-s must be a positive odd integer.")
    if args.k < 0 or args.k >= args.size:
        raise ValueError("-k must satisfy 0 <= k < -s.")
    if args.width < 1:
        raise ValueError("-w must be >= 1.")
    if args.gaps < 0:
        raise ValueError("-g must be >= 0.")
    if not 0 <= args.proportion <= 1:
        raise ValueError("-p must be in [0, 1].")
    if args.upstream < 0 or args.downstream < 0:
        raise ValueError("--up and --down must be >= 0.")


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


def _run_peak_pipeline(args: Namespace) -> None:
    """Run the selective ribosome profiling peak workflow."""
    serp = SeRP(args)
    serp.args_check()

    _print_step(2, "Import the RPF coverage data.")
    serp.rpf_txt_read()

    _print_step(3, "Import gene annotation.")
    serp.gene_anno()

    _print_step(4, "Detect selective ribosome profiling peaks.")
    serp.detect_binding_peaks()

    _print_step(5, "Output peak-detection results.")
    serp.output_peak()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for serp_peak."""
    now_time()
    print("\nDetect selective ribosome profiling peaks.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_peak_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
