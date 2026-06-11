#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_meta_codon.py


from utils.ribo.MetaCodon import *

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def meta_codon_plot_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the meta codon plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-c', dest="codon", required=False, type=str, help="the codon list in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-a', dest="around", required=False, type=int, default=20,
                        help="retrieve length of codon upstream and downstream. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="The number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="The number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('-u', dest="unique", action='store_true', required=False, default=False,
                        help="delete the cross repetition codon in different window. (default: %(default)s).")
    parser.add_argument('-s', dest="scale", action='store_true', required=False, default=False,
                        help="scale the window density with gene density. (default: %(default)s).")
    parser.add_argument('--smooth', dest="smooth", required=False, default=None, type=str,
                        help="smooth the window density [eg, 3,1]. (default: %(default)s).")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads (default: %(default)s).")
    parser.add_argument('--fig', dest="fig", required=False, action='store_true', default=False,
                        help="output the figure. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():
    now_time()

    print('\nDraw the meta-codon plot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = meta_codon_plot_args_parser()
    meta_codon = MetaCodon(args)

    print('\nStep2: Import the codon list.', flush=True)
    meta_codon.import_codon()

    print('\nStep3: Import the RPFs file.', flush=True)
    meta_codon.import_rpf()
    meta_codon.smooth_rpf_density()

    print('\nStep4: Retrieve specific codon density.', flush=True)
    meta_codon.reterieve_codon_density()

    print('\nStep5: Output the specific codon meta RPFs density.', flush=True)
    meta_codon.output_meta_codon_density()

    print('\nStep6: Output the specific codon sequence.', flush=True)
    meta_codon.output_meta_codon_seq()

    print('\nStep7: Draw the specific codon.', flush=True)
    if args.fig:
        meta_codon.draw_meta_codon()


    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
