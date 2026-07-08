#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Rensc
# Date: 2026-06-24
# Version: 0.2.7.8
# Function: This script is used to detect the P-site offset with RSBM.
# Input: Command-line arguments and input files specified by the user.
# Output: RSBM offset tables and heatmap figures.

"""Command-line entry point for rpf_Offset_RSBM.

This command detects P-site offset using RSBM directly from empirical
``exp_peak`` and ``exp_offset`` values.
"""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="This script is used to detect the P-site offset with RSBM."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-t",
        dest="transcript",
        required=True,
        type=str,
        help="the name of input transcript file in TXT format.",
    )
    required_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="the name of mapping file in BAM format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="the prefix of output file. (prefix + _RSBM_offset.txt)",
    )

    parser.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help=(
            "only retain the transcript with longest CDS of each gene "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=27,
        help="the minimum reads length to keep (default: %(default)s nt).",
    )
    parser.add_argument(
        "-M",
        dest="max",
        required=False,
        type=int,
        default=33,
        help="the maximum reads length to keep (default: %(default)s nt).",
    )
    parser.add_argument(
        "-p",
        dest="exp_peak",
        required=False,
        type=int,
        default=29,
        help="RPFs peak length [~30 nt] (default: %(default)s nt).",
    )
    parser.add_argument(
        "-e",
        dest="exp_offset",
        required=False,
        type=int,
        default=11,
        help="expected offset length (default: %(default)s).",
    )
    parser.add_argument(
        "-s",
        dest="shift",
        required=False,
        type=int,
        default=2,
        help=(
            "P-site shift for different RPF lengths. "
            "Empirical value: 2 nt for eukaryotes, 1 nt for prokaryotes "
            "(default: %(default)s nt)."
        ),
    )
    parser.add_argument(
        "--silence",
        dest="silence",
        required=False,
        action="store_true",
        default=False,
        help="discard warning information (default: %(default)s).",
    )
    parser.add_argument(
        "-d",
        dest="detail",
        action="store_true",
        required=False,
        default=False,
        help="output the details of offset (default: %(default)s).",
    )

    return parser


def _validate_args(parser: argparse.ArgumentParser, args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.transcript, args.bam)

    if args.min <= 0:
        parser.error("-m/--min must be greater than 0.")

    if args.max <= 0:
        parser.error("-M/--max must be greater than 0.")

    if args.max < args.min:
        parser.error("-M/--max must be greater than or equal to -m/--min.")

    if args.exp_peak <= 0:
        parser.error("-p/--exp_peak must be greater than 0.")

    if args.exp_offset <= 0:
        parser.error("-e/--exp_offset must be greater than 0.")

    if args.shift <= 0:
        parser.error("-s/--shift must be greater than 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse and validate command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(parser, args)
    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_offset_rsbm_pipeline(args: Namespace) -> None:
    """Run RSBM offset detection pipeline.

    Notes
    -----
    ``Offset`` is imported lazily here so that ``rpf_Offset_RSBM -h`` does not
    import heavy dependencies before argparse prints help information.
    """
    from utils.ribo.Offset_RSBM import Offset

    offset_attr = Offset(args)

    _print_step(1, "Import the transcripts annotation.")
    offset_attr.read_transcript()

    _print_step(2, "Import the BAM/SAM file.")
    offset_attr.get_mrna_reads()

    _print_step(3, "Detect the RSBM offset of sequence profile.")
    offset_attr.get_frame_offset()
    offset_attr.format_frame_offset()
    offset_attr.adjust_frame_offset()

    _print_step(4, "Output the frame offset.")
    offset_attr.write_frame_offset()

    _print_step(5, "Draw the frame offset heatmap.")
    offset_attr.draw_frame_heatmap()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Offset_RSBM."""
    args = _parse_args(argv)

    now_time()
    print("\nDetect the P-site offset.", flush=True)

    _run_offset_rsbm_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
