#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_digest.py


from utils.ribo.Digestion import *

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def digestion_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to Detect the digestion sites.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript file in TXT format.")
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the name of input transcript sequence file in FA format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the name of mapping file in BAM format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the name of output file. (prefix + _digestion_sites.txt)")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
                             "Recommended : True")
    parser.add_argument('--scale', dest="scale", action='store_true', required=False, default=False,
                        help="scale the motif matrix (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=20,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=100,
                        help="the maximum reads length to keep (default: %(default)s nt).")

    args = parser.parse_args()
    file_check(args.transcript, args.sequence, args.bam)
    args_print(args)

    return args


def main():
    print('\nDetect the digestion sites.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = digestion_args_parser()
    ribo_attr = Ribo(args)

    print('\nStep2: Import the annotation of transcripts.', flush=True)
    ribo_attr.read_transcript()

    print('\nStep3: Detect the digestion sites.', flush=True)
    ribo_attr.get_digest_sites()

    print('\nStep4: Output the digestion sites.', flush=True)
    ribo_attr.output_digest_sites()

    print('\nStep5: Draw the heatmap of digestion sites.', flush=True)
    ribo_attr.digestion_plot()

    print('\nStep6: Draw the seq logo of digestion sites.', flush=True)
    ribo_attr.output_counts()
    
    ribo_attr.seq_logo_plot2()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
