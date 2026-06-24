#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Rensc; modified for improved P-site offset detection
# Version: 0.2.7-offset-enhanced
# Function: Detect P-site offset using SSCBM, RSBM, or both.

import argparse

from utils.ribo.Offset import Offset
from utils.ribo.ArgsParser import args_print, file_check, now_time


def offset_args_parser():
    parser = argparse.ArgumentParser(
        description="Detect the P-site offset with SSCBM, RSBM, or both."
    )

    input_group = parser.add_argument_group("Required arguments")
    input_group.add_argument(
        "-t",
        dest="transcript",
        required=True,
        type=str,
        help="Input transcript annotation file in TXT format.",
    )
    input_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="Input mapping file in BAM/SAM format.",
    )
    input_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix. Tables and figures will use this prefix.",
    )

    parser.add_argument(
        "--mode",
        dest="mode",
        required=False,
        type=str,
        default="both",
        choices=["SSCBM", "RSBM", "both"],
        help="Offset detection mode: SSCBM, RSBM, or both. Default: %(default)s.",
    )
    # Deprecated. Kept hidden so old command lines do not fail, but it is no longer used.
    parser.add_argument(
        "-a",
        dest="align",
        required=False,
        type=str,
        default="deprecated",
        choices=["both", "tis", "tts", "deprecated"],
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help="Retain only the transcript with the longest CDS for each gene. Default: %(default)s.",
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=27,
        help="Minimum read length to keep. Default: %(default)s nt.",
    )
    parser.add_argument(
        "-M",
        dest="max",
        required=False,
        type=int,
        default=33,
        help="Maximum read length to keep. Default: %(default)s nt.",
    )
    parser.add_argument(
        "-p",
        dest="exp_peak",
        required=False,
        type=int,
        default=30,
        help="Expected RPF length fitted to ribosome structure. Default: %(default)s nt.",
    )
    parser.add_argument(
        "-s",
        dest="shift",
        required=False,
        type=int,
        default=2,
        help="Expected offset shift per read-length change. Default: %(default)s nt.",
    )
    parser.add_argument(
        "--screen",
        dest="screen",
        action="store_true",
        required=False,
        default=False,
        help="Filter unmapped, secondary, supplementary, duplicate, and low-MAPQ reads before offset detection.",
    )
    parser.add_argument(
        "--min-mapq",
        dest="min_mapq",
        required=False,
        type=int,
        default=10,
        help="Minimum MAPQ used when --screen is enabled. Default: %(default)s.",
    )
    parser.add_argument(
        "--min-offset-rpfs",
        dest="min_offset_rpfs",
        required=False,
        type=int,
        default=50,
        help="Minimum reads supporting one read length for confident offset detection. Default: %(default)s.",
    )
    parser.add_argument(
        "--dp-penalty",
        dest="dp_penalty",
        required=False,
        type=float,
        default=0.18,
        help="Transition penalty used by dynamic-programming smoothing. Default: %(default)s.",
    )
    parser.add_argument(
        "--silence",
        dest="silence",
        required=False,
        action="store_true",
        default=False,
        help="Discard warning information during parsing. Default: %(default)s.",
    )
    parser.add_argument(
        "-d",
        dest="detail",
        action="store_true",
        required=False,
        default=False,
        help="Output detailed offset profiles. Default: %(default)s.",
    )

    args = parser.parse_args()
    file_check(args.transcript, args.bam)
    args_print(args)
    return args


def main():
    now_time()
    print("\nDetect the P-site offset.", flush=True)

    print("\nStep1: Checking the input arguments.", flush=True)
    args = offset_args_parser()
    offset_attr = Offset(args)

    print("\nStep2: Import the transcripts annotation.", flush=True)
    offset_attr.read_transcript()

    print("\nStep3: Import the BAM/SAM file.", flush=True)
    offset_attr.get_mrna_reads()

    if args.mode in {"SSCBM", "both"}:
        print("\nStep4: Detect the SSCBM offset from TIS/TTS profiles.", flush=True)
        offset_attr.get_tis_offset()
        offset_attr.adjust_tis_offset()
        offset_attr.write_tis_offset()
        offset_attr.draw_tis_heatmap()

    if args.mode in {"RSBM", "both"}:
        print("\nStep5: Detect the RSBM offset from CDS frame periodicity.", flush=True)
        offset_attr.get_frame_offset()
        offset_attr.format_frame_offset()
        offset_attr.adjust_frame_offset()
        offset_attr.write_frame_offset()
        offset_attr.draw_frame_heatmap()

    now_time()
    print("\nAll done.", flush=True)


if __name__ == "__main__":
    main()
