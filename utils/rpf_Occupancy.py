#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/11/30 21:46
# @Project : riboParser
# @Script  : rpf_occupancy.py

from utils.ribo import Occupancy

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def rpf_occupancy_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the codon occupancy plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the RPFs calculation
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=30,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=15,
                        help="The number of codons after TIS will be discarded. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=5,
                        help="The number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the occupancy. (default: %(default)s).")
    parser.add_argument('--stop', dest="stop",  action='store_true', required=False, default=False,
                        help="rmove the stop codon. (default: %(default)s).")
    parser.add_argument('--all', dest="all", action='store_true', required=False, default=False,
                        help="output all RPFs density. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():
    now_time()
    print('\nCalculate the codon occupancy.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = rpf_occupancy_args_parser()
    occupancy = Occupancy.Occupancy(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    occupancy.import_rpf()

    print('\nStep3: Calculate the codon occupancy.', flush=True)
    occupancy.codon_occupancy()

    print('\nStep4: Draw the codon occupancy plot.', flush=True)
    occupancy.draw_occupancy_corr()
    occupancy.draw_occupancy_heat()
    occupancy.draw_occupancy_relative_heat()
    occupancy.draw_occupancy_line()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
