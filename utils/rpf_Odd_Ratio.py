#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-06
# Version: dev003
# Function: Detect local ribosome pauses and differential codon enrichment.
# Input: Frame-resolved RPF density file and control/treatment sample names.
# Output: Local-pause, differential-enrichment, summary, and figure files.

"""Command-line entry point for codon-level odds-ratio analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Odd_Ratio

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
            "Detect local ribosome pauses and differential codon enrichment."
        )
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input frame-resolved RPF density file.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output file prefix.",
    )
    required_group.add_argument(
        "-c",
        "--control",
        dest="control",
        required=True,
        type=str,
        help="Comma-separated control sample names.",
    )
    required_group.add_argument(
        "-t",
        "--treat",
        dest="treat",
        required=True,
        type=str,
        help="Comma-separated treatment sample names.",
    )

    filtering_group = parser.add_argument_group("Input filtering arguments")
    filtering_group.add_argument(
        "-l",
        "--list",
        dest="list",
        default=None,
        type=str,
        help="Optional transcript ID list.",
    )
    filtering_group.add_argument(
        "-s",
        "--site",
        dest="site",
        choices=["E", "P", "A"],
        default="P",
        type=str,
        help="Ribosome site used for coordinate assignment (default: %(default)s).",
    )
    filtering_group.add_argument(
        "-f",
        "--frame",
        dest="frame",
        choices=["0", "1", "2", "all"],
        default="all",
        type=str,
        help="Reading frame used for calculation (default: %(default)s).",
    )
    filtering_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=50,
        type=int,
        help="Minimum RPF count required for retained transcripts (default: %(default)s).",
    )
    filtering_group.add_argument(
        "--tis",
        dest="tis",
        default=0,
        type=int,
        help="Number of codons removed after TIS (default: %(default)s).",
    )
    filtering_group.add_argument(
        "--tts",
        dest="tts",
        default=0,
        type=int,
        help="Number of codons removed before TTS (default: %(default)s).",
    )
    filtering_group.add_argument(
        "--stop",
        dest="stop",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Exclude stop codons from the analysis. Use --no-stop to include "
            "stop codons (default: %(default)s)."
        ),
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--detail",
        dest="detail",
        action="store_true",
        default=False,
        help=(
            "Output all tested positions to "
            "<prefix>_codon_odd_ratio_all.txt (default: %(default)s)."
        ),
    )

    calculation_group = parser.add_argument_group("Calculation arguments")
    calculation_group.add_argument(
        "--thread",
        dest="thread",
        default=1,
        type=int,
        help=(
            "Compatibility option; core statistics are vectorized and do not "
            "spawn per-position workers (default: %(default)s)."
        ),
    )
    calculation_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        default=False,
        help=(
            "Append per-sample RPM columns for reporting; statistical tests "
            "always use raw counts (default: %(default)s)."
        ),
    )
    calculation_group.add_argument(
        "-z",
        "--zero",
        dest="zero",
        action="store_true",
        default=False,
        help=(
            "Deprecated compatibility flag. Zero-density values are retained "
            "automatically and are never replaced in the raw count matrix."
        ),
    )
    calculation_group.add_argument(
        "--test",
        dest="test",
        choices=["auto", "beta-binomial", "pooled"],
        default="auto",
        type=str,
        help=(
            "Differential test. 'auto' uses a replicate-aware beta-binomial "
            "Wald test when both groups contain >=2 samples, otherwise a pooled "
            "Wald test (default: %(default)s)."
        ),
    )
    calculation_group.add_argument(
        "--pseudocount",
        dest="pseudocount",
        default=0.5,
        type=float,
        help=(
            "Haldane-Anscombe pseudocount used only for log2 odds-ratio "
            "estimation; raw RPF counts are unchanged (default: %(default)s)."
        ),
    )
    calculation_group.add_argument(
        "--min-log2-or",
        dest="min_log2_or",
        default=1.0,
        type=float,
        help=(
            "Minimum absolute log2 odds ratio for differential calls "
            "(default: %(default)s)."
        ),
    )
    calculation_group.add_argument(
        "--fdr",
        dest="fdr",
        choices=["bhfdr", "pvalue"],
        default="bhfdr",
        type=str,
        help="Significance field used for filtering (default: %(default)s).",
    )
    calculation_group.add_argument(
        "-v",
        "--value",
        dest="value",
        default=0.05,
        type=float,
        help="Significance threshold (default: %(default)s).",
    )
    calculation_group.add_argument(
        "--scale",
        dest="scale",
        choices=["zscore", "minmax"],
        default="minmax",
        type=str,
        help="Scaling method used for visualization (default: %(default)s).",
    )

    pause_group = parser.add_argument_group("Local pause arguments")
    pause_group.add_argument(
        "--local-window",
        dest="local_window",
        default=101,
        type=int,
        help=(
            "Centered local-background window in codons; must be odd "
            "(default: %(default)s)."
        ),
    )
    pause_group.add_argument(
        "--pause-score",
        dest="pause_score",
        default=20.0,
        type=float,
        help=(
            "Minimum site-to-local-mean ratio for a local pause "
            "(default: %(default)s)."
        ),
    )
    pause_group.add_argument(
        "--min-site-rpf",
        dest="min_site_rpf",
        default=3,
        type=int,
        help=(
            "Minimum site RPF count required in the pooled group and in each "
            "supporting replicate (default: %(default)s)."
        ),
    )
    pause_group.add_argument(
        "--min-local-coverage",
        dest="min_local_coverage",
        default=0.10,
        type=float,
        help=(
            "Minimum fraction of covered codons in the local window "
            "(default: %(default)s)."
        ),
    )
    pause_group.add_argument(
        "--min-pause-replicates",
        dest="min_pause_replicates",
        default=2,
        type=int,
        help=(
            "Minimum supporting replicates for a group-local pause; capped at "
            "the available group size (default: %(default)s)."
        ),
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.list:
        file_check(args.list)
    if args.thread < 1:
        raise ValueError("--thread must be >= 1.")
    if args.min < 0:
        raise ValueError("-m must be >= 0.")
    if args.tis < 0 or args.tts < 0:
        raise ValueError("--tis and --tts must be >= 0.")
    if not 0 < args.value <= 1:
        raise ValueError("-v must be in the interval (0, 1].")
    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be > 0.")
    if args.min_log2_or < 0:
        raise ValueError("--min-log2-or must be >= 0.")
    if args.local_window < 3 or args.local_window % 2 == 0:
        raise ValueError("--local-window must be an odd integer >= 3.")
    if args.pause_score <= 1:
        raise ValueError("--pause-score must be > 1.")
    if args.min_site_rpf < 1:
        raise ValueError("--min-site-rpf must be >= 1.")
    if not 0 <= args.min_local_coverage <= 1:
        raise ValueError("--min-local-coverage must be in the interval [0, 1].")
    if args.min_pause_replicates < 1:
        raise ValueError("--min-pause-replicates must be >= 1.")
    control_samples = [item.strip() for item in args.control.split(",") if item.strip()]
    treat_samples = [item.strip() for item in args.treat.split(",") if item.strip()]
    if not control_samples or not treat_samples:
        raise ValueError("-c/--control and -t/--treat must each contain a sample.")
    if len(control_samples) != len(set(control_samples)):
        raise ValueError("Duplicate sample names were found in -c/--control.")
    if len(treat_samples) != len(set(treat_samples)):
        raise ValueError("Duplicate sample names were found in -t/--treat.")
    overlap = sorted(set(control_samples).intersection(treat_samples))
    if overlap:
        raise ValueError(
            "Samples cannot occur in both groups: {samples}.".format(
                samples=", ".join(overlap)
            )
        )
    if args.test == "beta-binomial":
        if min(len(control_samples), len(treat_samples)) < 2:
            raise ValueError(
                "--test beta-binomial requires at least two samples in each group."
            )


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_odd_ratio_pipeline(args: Namespace) -> None:
    """Run the codon-level odds-ratio workflow."""
    odd_ratio = Odd_Ratio.OddRatio(args)

    step_print(2, "Import the RPF density file.")
    odd_ratio.read_rpf()

    step_print(3, "Build site-versus-gene count tables.")
    odd_ratio.make_two_dimensional_table()

    step_print(4, "Calculate local pause evidence.")
    odd_ratio.calc_local_pause()

    step_print(5, "Calculate differential codon enrichment.")
    odd_ratio.calc_odd_ratio()
    odd_ratio.calc_differential_test()
    odd_ratio.classify_pause()

    step_print(6, "Output differential and local-pause results.")
    odd_ratio.output_odd_ratio()
    odd_ratio.summarize_odd_ratio()

    step_print(7, "Draw codon differential-enrichment figures.")
    odd_ratio.draw_odd_ratio_line()
    odd_ratio.draw_odd_ratio_scatter()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Odd_Ratio."""
    now_time()
    title_print('Detect local pauses and differential codon enrichment.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_odd_ratio_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()