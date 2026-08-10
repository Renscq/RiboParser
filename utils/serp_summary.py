#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-09
# Version: dev002
# Function: Summarize SeRP peak-calling or pairwise peak-overlap results.
# Input: serp_peak or serp_overlap output prefix.
# Output: Peak statistics, shared/specific summaries, significance tables, and figures.

"""Command-line entry point for SeRP peak and overlap-result summarization."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence
from pathlib import Path

from utils.ribo.ArgsParser import (
    args_print,
    complete_print,
    file_check,
    now_time,
    result_print,
    step_print,
    title_print,
)


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Summarize either one serp_peak result or one pairwise serp_overlap "
            "comparison. Auto mode detects the result type from the input prefix."
        ),
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        type=str,
        help="Input prefix previously used by serp_peak or serp_overlap.",
    )

    mode_group = parser.add_argument_group("Analysis mode")
    mode_group.add_argument(
        "--mode",
        dest="mode",
        choices=["auto", "peak", "overlap"],
        default="auto",
        help="Summary mode. Auto detects serp_peak versus serp_overlap outputs.",
    )

    significance_group = parser.add_argument_group("Overlap significance arguments")
    significance_group.add_argument(
        "--significance-metric",
        dest="significance_metric",
        choices=["bhfdr", "pvalue", "none"],
        default="bhfdr",
        help=(
            "Metric used only to annotate statistical significance in overlap mode. "
            "Missing values are classified as untested."
        ),
    )
    significance_group.add_argument(
        "--significance-cutoff",
        dest="significance_cutoff",
        default=0.05,
        type=float,
        help="Maximum BHFDR/P-value annotated as significant in overlap mode.",
    )
    significance_group.add_argument(
        "--top-specific-genes",
        dest="top_specific_genes",
        default=20,
        type=int,
        help="Maximum genes shown in the condition-specific peak-burden figure.",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        type=str,
        help="Output prefix; defaults to <input>_summary.",
    )
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["pdf", "png", "both"],
        default="pdf",
        help="Summary figure output format.",
    )
    output_group.add_argument(
        "--plot",
        dest="plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Generate summary figures in addition to summary tables.",
    )

    peak_group = parser.add_argument_group("Peak-mode positional map arguments")
    peak_group.add_argument(
        "--bins",
        dest="bins",
        default=100,
        type=int,
        help="Number of normalized transcript bins used in peak mode.",
    )
    peak_group.add_argument(
        "--max-heatmap-cells",
        dest="max_heatmap_cells",
        default=100_000_000,
        type=int,
        help="Maximum dense full-length heatmap cells used in peak mode.",
    )
    peak_group.add_argument(
        "--chunksize",
        dest="chunksize",
        default=500_000,
        type=int,
        help="Rows read per chunk from *_peaks_ratio.txt in peak mode.",
    )

    figure_group = parser.add_argument_group("Figure arguments")
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


def _peak_inputs_exist(prefix: str) -> bool:
    """Return whether the required serp_peak summary inputs exist."""
    return Path(prefix + "_peaks.txt").is_file() and Path(
        prefix + "_peaks_ratio.txt"
    ).is_file()


def _overlap_inputs_exist(prefix: str) -> bool:
    """Return whether the core serp_overlap summary inputs exist."""
    required = [
        prefix + ".summary.txt",
        prefix + ".shared.peaks.txt",
        prefix + ".shared.clusters.txt",
        prefix + ".relationships.txt",
    ]
    return all(Path(path).is_file() for path in required)


def _resolve_mode(args: Namespace) -> str:
    """Resolve peak versus overlap summary mode from explicit or detected inputs."""
    peak_exists = _peak_inputs_exist(args.input)
    overlap_exists = _overlap_inputs_exist(args.input)

    if args.mode == "peak":
        if not peak_exists:
            file_check(args.input + "_peaks.txt", args.input + "_peaks_ratio.txt")
        return "peak"
    if args.mode == "overlap":
        if not overlap_exists:
            file_check(
                args.input + ".summary.txt",
                args.input + ".shared.peaks.txt",
                args.input + ".shared.clusters.txt",
                args.input + ".relationships.txt",
            )
        return "overlap"

    if peak_exists and overlap_exists:
        raise ValueError(
            "Both serp_peak and serp_overlap outputs match this prefix. "
            "Specify --mode peak or --mode overlap explicitly."
        )
    if overlap_exists:
        return "overlap"
    if peak_exists:
        return "peak"
    raise FileNotFoundError(
        "Cannot detect SeRP summary mode from prefix: {0}. Expected either "
        "*_peaks.txt + *_peaks_ratio.txt or .summary.txt + .shared.peaks.txt + "
        ".shared.clusters.txt + .relationships.txt.".format(args.input)
    )


def _validate_args(args: Namespace) -> None:
    """Validate shared and mode-specific numeric arguments."""
    if args.output is None:
        args.output = args.input + "_summary"

    output_parent = Path(args.output).expanduser().resolve().parent
    if not output_parent.exists():
        raise FileNotFoundError(
            "Output directory does not exist: {0}".format(output_parent)
        )
    if not 0 <= args.significance_cutoff <= 1:
        raise ValueError("--significance-cutoff must be within [0, 1].")
    if args.top_specific_genes < 0:
        raise ValueError("--top-specific-genes must be >= 0.")
    if args.bins < 10:
        raise ValueError("--bins must be >= 10.")
    if args.chunksize < 1:
        raise ValueError("--chunksize must be >= 1.")
    if args.max_heatmap_cells < 1:
        raise ValueError("--max-heatmap-cells must be >= 1.")
    if args.font_size <= 0:
        raise ValueError("--font-size must be > 0.")
    if args.dpi < 1:
        raise ValueError("--dpi must be >= 1.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, resolve, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        args.mode = _resolve_mode(args)
        _validate_args(args)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    args_print(args)
    return args


def _run_peak_summary(args: Namespace) -> None:
    """Run the existing single-serp_peak summary workflow."""
    from utils.serp.Summary import SeRPSummary

    workflow = SeRPSummary(args)

    step_print(2, "Read called peaks and transcript spans.")
    workflow.read_peaks()
    workflow.read_transcript_spans()

    step_print(3, "Calculate peak and transcript summary statistics.")
    workflow.calculate_summary()
    workflow.build_normalized_matrix()

    step_print(4, "Write SeRP peak summary tables.")
    workflow.write_tables()

    if args.plot:
        step_print(5, "Draw SeRP peak summary figures.")
        workflow.draw_figures()


def _run_overlap_summary(args: Namespace) -> None:
    """Run pairwise shared/specific SeRP overlap-result summarization."""
    from utils.serp.OverlapSummary import SeRPOverlapSummary

    workflow = SeRPOverlapSummary(args)

    step_print(2, "Read shared, specific, cluster, and relationship tables.")
    workflow.read_results()

    step_print(3, "Classify significance and summarize group-specific peaks.")
    workflow.calculate_summary()

    step_print(4, "Write overlap, significance, transcript, and gene summary tables.")
    workflow.write_tables()

    if args.plot:
        step_print(5, "Draw shared, specific, significance, and overlap summary figures.")
        workflow.draw_figures()

    result_print(workflow.result_counts())


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for ``serp_summary``."""
    now_time()
    title_print("Summarize selective ribosome profiling peaks.")

    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)

    if args.mode == "overlap":
        _run_overlap_summary(args)
    else:
        _run_peak_summary(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()