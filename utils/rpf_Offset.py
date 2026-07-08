#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Rensc; modified for improved P-site offset detection
# Date: 2026-06-24
# Version: 0.2.7.8
# Function: Detect P-site offset using SSCBM followed by RSBM.
# Input: Command-line arguments and input files specified by the user.
# Output: SSCBM and RSBM offset tables and heatmap figures.

"""Command-line entry point for rpf_Offset.

This command detects P-site offset using a two-step workflow:

1. SSCBM: infer offset from read-end enrichment around TIS/TTS.
2. RSBM: refine offset with reading-frame periodicity.
"""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="This script is used to detect the P-site offset."
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
        help="the prefix of output file. (prefix + _offset.txt)",
    )

    parser.add_argument(
        "-a",
        dest="align",
        required=False,
        type=str,
        default="both",
        choices=["both", "tis", "tts"],
        help=(
            "specify the alignment of reads for offset detection "
            "[both, tis, tts] (default: %(default)s)."
        ),
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
        default=30,
        help=(
            "Expected RPF length fitted to ribosome structure [~30 nt] "
            "(default: %(default)s nt)."
        ),
    )
    parser.add_argument(
        "-s",
        dest="shift",
        required=False,
        type=int,
        default=2,
        help="P-site shift for different RPF lengths (default: %(default)s nt).",
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
    """Print a standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_offset_pipeline(args: Namespace) -> None:
    """Run SSCBM followed by RSBM offset detection.

    Notes
    -----
    ``Offset`` is imported lazily here so that ``rpf_Offset -h`` does not
    import heavy dependencies before argparse prints help information.
    """
    from utils.ribo.Offset import Offset

    offset_attr = Offset(args)

    _print_step(1, "Import the transcripts annotation.")
    offset_attr.read_transcript()

    _print_step(2, "Import the BAM/SAM file.")
    offset_attr.get_mrna_reads()

    _print_step(3, "Detect the SSCBM offset of sequence profile.")
    offset_attr.get_tis_offset()
    offset_attr.adjust_tis_offset()
    offset_attr.write_tis_offset()
    offset_attr.draw_tis_heatmap()

    _print_step(4, "Detect the RSBM offset of sequence profile.")
    offset_attr.get_frame_offset()
    offset_attr.format_frame_offset()
    offset_attr.adjust_frame_offset()
    offset_attr.write_frame_offset()
    offset_attr.draw_frame_heatmap()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Offset."""
    args = _parse_args(argv)

    now_time()
    print("\nDetect the P-site offset.", flush=True)

    _run_offset_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
