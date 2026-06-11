#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_density.py


from utils.ribo import Ribo

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def ribo_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to convert Ribo-seq bam to p-site density.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript file in TXT format.")
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the name of input transcript sequence file in FA format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, 
                             help="the name of mapping file in BAM format.")
    input_group.add_argument('-p', dest="psite", required=True, type=str,
                             help="the name of p-site offset file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (output = prefix + _rpf.txt)")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
                             "Recommended : True")
    parser.add_argument('-m', dest="min", required=False, type=int, default=27,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=33,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    input_group.add_argument('--period', dest='periodicity', required=False, type=float, default=40,
                             help="the minimum 3nt periodicity to keep. (default: %(default)s).")
    parser.add_argument('--silence', dest="silence", required=False, action='store_true', default=True,
                        help="discard the warning information. (default: %(default)s).")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads (default: %(default)s). It will take a lot of memory.")
    args = parser.parse_args()
    file_check(args.transcript, args.sequence, args.bam, args.psite)
    args_print(args)

    return args


def main():
    now_time()
    print('\nConvert reads to RPFs density.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ribo_args_parser()

    ribo_attr = Ribo.Ribo(args)

    print('\nStep2: Import the P-site offset.', flush=True)
    ribo_attr.read_offset()

    print('\nStep3: Import the transcripts annotation.', flush=True)
    ribo_attr.read_transcript()
    ribo_attr.check_transcript()

    print('\nStep4: Import the BAM file.', flush=True)
    ribo_attr.read_bam()

    print('\nStep5: Format the in-frame RPFs density.', flush=True)
    ribo_attr.run_format_rpf_with_multi_thread()

    print('\nStep5: Output the RPFs density.', flush=True)
    ribo_attr.output_density()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
