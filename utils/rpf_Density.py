#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.008
# Function: Convert Ribo-seq BAM/SAM alignments to multi-sample-compatible self-contained P-site density JSONL.
# Input: RiboParser norm TXT or genePred annotation, transcript FASTA, BAM/SAM alignment, and P-site offset table.
# Output: Gzip-compressed compact RPF density JSONL, summary JSON, and optional legacy TXT table.

"""Command-line entry point for compact RPF density generation."""

from __future__ import annotations

import argparse

from utils.ribo import Ribo
from utils.ribo.ArgsParser import args_print, file_check, now_time


def ribo_args_parser():
    parser = argparse.ArgumentParser(
        description="Convert Ribo-seq BAM/SAM alignments to multi-sample-compatible self-contained P-site density JSONL."
    )

    input_group = parser.add_argument_group("Required arguments")
    input_group.add_argument(
        "-t",
        dest="transcript",
        required=True,
        type=str,
        help="Input transcript annotation file in RiboParser norm TXT or genePred format.",
    )
    input_group.add_argument(
        "-s",
        dest="sequence",
        required=True,
        type=str,
        help="Input transcript sequence file in FASTA format.",
    )
    input_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="Input transcriptome-aligned mapping file in BAM or SAM format.",
    )
    input_group.add_argument(
        "-p",
        dest="psite",
        required=True,
        type=str,
        help="Input P-site offset file in TXT format.",
    )
    input_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix. Default JSON output is prefix + '_rpf.jsonl.gz'. The JSON uses a samples field to support downstream multi-sample merging.",
    )

    parser.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help="Only retain the transcript with longest CDS of each gene. Recommended: True. Default: %(default)s.",
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
        "--period",
        dest="periodicity",
        required=False,
        type=float,
        default=40,
        help="Minimum 3-nt periodicity to keep. Default: %(default)s.",
    )
    parser.add_argument(
        "--min-confidence",
        dest="min_confidence",
        required=False,
        type=float,
        default=50.0,
        help="Minimum offset confidence to keep when the offset file contains a confidence column. Default: %(default)s.",
    )
    parser.add_argument(
        "--drop-warning",
        dest="drop_warning",
        required=False,
        action="store_true",
        default=False,
        help="Drop offset rows whose warning column is not PASS when the offset file contains a warning column. Default: %(default)s.",
    )
    parser.add_argument(
        "--silence",
        dest="silence",
        required=False,
        action="store_true",
        default=False,
        help="Discard warning information. Default: %(default)s.",
    )
    parser.add_argument(
        "--thread",
        dest="thread",
        type=int,
        required=False,
        default=1,
        help="Number of workers for indexed BAM parallel scanning. If BAM is not indexed, stream mode is used. Default: %(default)s.",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["json", "txt", "both"],
        default="json",
        help="Output format. 'json' writes gzip-compressed compact JSONL; 'txt' writes legacy TXT; 'both' writes both. Default: %(default)s.",
    )
    output_group.add_argument(
        "--density-encoding",
        dest="density_encoding",
        choices=["sparse", "dense"],
        default="sparse",
        help="Density encoding in JSONL. Sparse encoding is recommended. Default: %(default)s.",
    )

    args = parser.parse_args()

    file_check(args.transcript, args.sequence, args.bam, args.psite)
    args_print(args)
    return args


def main():
    now_time()
    print("\nConvert reads to multi-sample-compatible compact RPF density.", flush=True)
    print("\nStep1: Checking the input arguments.", flush=True)
    args = ribo_args_parser()

    ribo_attr = Ribo.Ribo(args)

    print("\nStep2: Import the P-site offset.", flush=True)
    ribo_attr.read_offset()

    print("\nStep3: Import the transcripts annotation.", flush=True)
    ribo_attr.read_transcript()
    ribo_attr.check_transcript()

    print("\nStep4: Import the BAM/SAM file and count P-site density.", flush=True)
    ribo_attr.read_bam()

    print("\nStep5: Output the RPF density.", flush=True)
    ribo_attr.output_density()

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
