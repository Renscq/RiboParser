#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Rensc, modified by local fork
# Date: 2026-06-22
# Version: 0.2.7
# Function: This script is used to summary the BAM condition.

import argparse

from utils.ribo.ArgsParser import args_print, file_check, now_time
from utils.ribo.Quality import *


def rpf_bam_check_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="This script is used to summary the BAM condition.",
    )

    # arguments for the Required arguments
    input_group = parser.add_argument_group("Required arguments")
    input_group.add_argument(
        "-t",
        dest="transcript",
        required=True,
        type=str,
        help="the input file name of gene annotation",
    )
    input_group.add_argument(
        "-b",
        dest="bam",
        required=True,
        type=str,
        help="the input file name of bam",
    )
    input_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="the prefix of output file.",
    )

    # arguments for the modification
    parser.add_argument(
        "--thread",
        dest="thread",
        required=False,
        type=int,
        default=1,
        help="""the number of threads (default: %(default)s).
Suitable for large bam files > 1G. It will take a lot of memory.""",
    )
    parser.add_argument(
        "-g",
        dest="tag",
        choices=[0, 1],
        type=int,
        required=False,
        default=0,
        help="""filter the number of reads mapped loci (default: %(default)s).
[0]: all reads will be used;
[1]: reads with unique mapped loci will be used.""",
    )
    parser.add_argument(
        "-a",
        dest="align",
        choices=["star", "hisat2", "bowtie2"],
        type=str,
        required=False,
        default="star",
        help="""aligner used to generate the BAM/SAM file.
This option controls unique-read tag interpretation. (default: %(default)s).""",
    )
    parser.add_argument(
        "-r",
        dest="reverse",
        action="store_true",
        required=False,
        default=False,
        help="reads aligned to negative strand will also be counted. (default: %(default)s).",
    )
    parser.add_argument(
        "-l",
        dest="longest",
        action="store_true",
        required=False,
        default=False,
        help="only keep the representative/longest transcript per gene. (default: %(default)s).",
    )
    parser.add_argument(
        "-s",
        dest="saturation",
        action="store_true",
        required=False,
        default=False,
        help="""whether to calculate RPF saturation. (default: %(default)s).
This step will take extra time and memory.
The saturation random seed is fixed to 5201314 for reproducibility.""",
    )

    args = parser.parse_args()
    file_check(args.transcript, args.bam)
    args_print(args)

    return args


def main():
    now_time()
    print("\nCheck the RPFs mapping condition.", flush=True)

    print("\nStep1: Checking the input Arguments.", flush=True)
    args = rpf_bam_check_parser()

    rpf_quality = Quality(args)

    print("\nStep2: Import the transcripts annotation.", flush=True)
    rpf_quality.read_transcript()

    print("\nStep3: Import the bam file.", flush=True)
    rpf_quality.sort_index_bam()
    rpf_quality.fliter_mrna_reads()

    print("\nStep4: Sort and index the bam file.", flush=True)
    rpf_quality.merge_sort_index_bam()

    print("\nStep5: Detect the type of sequence profile.", flush=True)
    if not rpf_quality.profile:
        rpf_quality.detect_seq_type()
    else:
        print(
            "{bam} is specified as {ribo}-seq.".format(
                bam=rpf_quality.sample_file,
                ribo=rpf_quality.profile,
            ),
            flush=True,
        )

    print("\nStep6: Summary the length distribution of reads aligned to mRNA.", flush=True)
    rpf_quality.write_length_distr()
    rpf_quality.write_filter_stats()
    rpf_quality.write_summary()

    if args.saturation:
        print("\nStep7: Check the RPFs saturation.", flush=True)
        rpf_quality.rpf_saturation()
        rpf_quality.draw_gene_saturation()
        rpf_quality.draw_rpf_saturation()

        # Summary is written again after saturation so saturation-specific fields are complete.
        rpf_quality.write_summary()
    else:
        print("\nStep7: Do not check the RPFs saturation.", flush=True)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
