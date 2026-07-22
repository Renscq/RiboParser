#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.003
# Function: Quantify RPF abundance at gene level from RPF density files.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Region-level RPF count, RPM, RPKM, TPM tables, QC plots, and summary JSON.

"""Command-line entry point for RPF gene-level quantification."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Quant

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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Quantify gene-level RPF abundance from JSONL or TXT density files.",
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

    input_group = parser.add_argument_group("Input and density arguments")
    input_group.add_argument(
        "-t",
        dest="transcript",
        required=False,
        default=None,
        type=str,
        help="Optional transcript filter table. If provided, transcript_id is used when available.",
    )
    input_group.add_argument(
        "-f",
        dest="frame",
        choices=["0", "1", "2", "all"],
        required=False,
        type=str,
        default="all",
        help="Reading frame used for quantification. Default: %(default)s.",
    )
    input_group.add_argument(
        "--region",
        dest="region",
        choices=Quant.REGION_CHOICES,
        required=False,
        type=str,
        default="cds",
        help="Transcript region to quantify. Use 'all' to quantify 5'UTR, CDS, and 3'UTR. Default: %(default)s.",
    )

    trim_group = parser.add_argument_group("CDS trimming arguments")
    trim_group.add_argument(
        "--tis",
        dest="tis",
        required=False,
        type=int,
        default=0,
        help="Number of codons after TIS to discard for CDS quantification. Default: %(default)s.",
    )
    trim_group.add_argument(
        "--tts",
        dest="tts",
        required=False,
        type=int,
        default=0,
        help="Number of codons before TTS to discard for CDS quantification. Default: %(default)s.",
    )

    filter_group = parser.add_argument_group("Filtering arguments")
    filter_group.add_argument(
        "--remove-outlier",
        "--outlier",
        dest="remove_outlier",
        action="store_true",
        required=False,
        default=False,
        help="Remove isolated codon-level RPF pileups before gene-level quantification. Default: %(default)s.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)

    if args.tis < 0:
        raise ValueError("--tis must be >= 0.")
    if args.tts < 0:
        raise ValueError("--tts must be >= 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_quant_pipeline(args: Namespace) -> None:
    """Run the RPF quantification workflow."""
    quant = Quant.Quant(args)

    step_print(2, "Import the RPF density file.")
    quant.import_rpf()

    step_print(3, "Quantify RPF abundance.")
    quant.quantify_regions()
    quant.output_total_rpf()

    step_print(4, "Draw quantification QC plots.")
    quant.draw_rpf_barplot()
    quant.draw_rpf_cdfplot()
    quant.draw_rpf_pcaplot()
    quant.draw_rpf_heatmap()

    step_print(5, "Write summary JSON.")
    quant.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Quant."""
    now_time()
    title_print('Quantify RPF abundance at gene level.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_quant_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
