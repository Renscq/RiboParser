#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.7.5
# Function: This script is used to draw the periodicity plot.
# Input: RPF density file in JSONL or TXT format.
# Output: Periodicity summary table and periodicity figures.

"""Command-line entry point for 3-nt periodicity analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Periodicity
from utils.ribo.ArgsParser import args_print, file_check, now_time


def periodicity_args_parser(argv: Sequence[str] | None = None) -> Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="This script is used to draw the periodicity plot.")

    input_group = parser.add_argument_group("Required arguments")
    input_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="the name of input RPF density file in JSONL or TXT format.",
    )
    input_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="the prefix of output file.",
    )

    input_group.add_argument(
        "-t",
        dest="transcript",
        required=False,
        type=str,
        default=None,
        help="the name of input transcript filter file in TXT format.",
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=50,
        help="retain transcript with more than minimum RPFs. (default: %(default)s).",
    )
    parser.add_argument(
        "--tis",
        dest="tis",
        required=False,
        type=int,
        default=0,
        help="the number of codons after TIS will be discarded. (default: %(default)s AA).",
    )
    parser.add_argument(
        "--tts",
        dest="tts",
        required=False,
        type=int,
        default=0,
        help="the number of codons before TTS will be discarded. (default: %(default)s AA).",
    )

    args = parser.parse_args(argv)
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)
    args_print(args)

    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Run the 3-nt periodicity analysis workflow."""
    now_time()
    print("\nDraw the periodicity plot.", flush=True)
    print("\nStep1: Checking the input Arguments.", flush=True)
    args = periodicity_args_parser(argv)

    print("\nStep2: Import the RPFs file.", flush=True)
    rpfs = Periodicity.Periodicity(args)
    rpfs.import_rpf()

    print("\nStep3: Calculate the 3nt periodicity.", flush=True)
    rpfs.calc_3nt_period()

    print("\nStep4: Output the 3nt periodicity.", flush=True)
    rpfs.output_meta()

    print("\nStep5: Draw the 3nt periodicity plot.", flush=True)
    rpfs.draw_3nt_period_count()
    rpfs.draw_3nt_period_ratio()
    rpfs.draw_3nt_period_stacked()

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
