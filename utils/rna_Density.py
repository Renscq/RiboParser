#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rna_density.py


from utils.ribo import RNA

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def rna_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to convert RNA-seq bam to p-site density.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript file in TXT format.")
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the name of input transcript sequence file in FA format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the name of mapping file in BAM format.")
    input_group.add_argument('-p', dest="psite", required=True, type=str,
                             help="the name of p-site offset file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _rna.txt)")

    # arguments for the RNA-seq parsing
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
                             " Recommended : True")
    parser.add_argument('-m', dest="min", required=False, type=int, default=25,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=150,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-r', dest="rolling", action='store_true', required=False, default=False,
                        help="density is calculated once per ribosome width. (default: %(default)s)."
                             "only suitable for long RNA-seq reads.")
    parser.add_argument('--pe', dest="pair_end", action='store_true', required=False, default=False,
                        help="reads aligned to negative strand will also be counted. (default: %(default)s)."
                             "only suitable for pair-end RNA-seq reads.")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads (default: %(default)s). It will take a lot of memory.")
    args = parser.parse_args()
    file_check(args.transcript, args.sequence, args.bam, args.psite)
    args_print(args)

    return args


def main():
    now_time()
    print('\nConvert reads to reads density.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = rna_args_parser()

    rna_attr = RNA.RNA(args)

    print('\nStep2: Import the P-site offset.', flush=True)
    rna_attr.read_offset()

    print('\nStep3: Import the transcripts annotation.', flush=True)
    rna_attr.read_transcript()

    print('\nStep4: Import the BAM file.', flush=True)
    rna_attr.read_bam()
    rna_attr.calculate_density()

    print('\nStep5: Format the in-frame reads density.', flush=True)
    rna_attr.run_format_reads_with_multi_thread()

    print('\nStep6: Output the reads density.', flush=True)
    rna_attr.output_density()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
