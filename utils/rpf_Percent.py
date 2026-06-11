#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_coverage.py


from utils.ribo import Percentage

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def percentage_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the plot of rpf coverage percentage.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript filein TXT format.")
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, help="the prefix of output file.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():
    now_time()
    print('\nDraw the metagene coverage.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = percentage_args_parser()
    percent = Percentage.Percentage(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    percent.read_rpf()
    percent.import_gene()

    print('\nStep3: Calculate the percentage of RPFs coverage.', flush=True)
    percent.calc_density_percent()

    print('\nStep4: Draw the histogram and boxplot of RPFs coverage.', flush=True)
    percent.draw_rpf_histogram()
    percent.draw_rpf_boxplot()

    print('\nStep5: Output the RPFs coverage.', flush=True)
    percent.output_density_percent()

    print('\nAll done.\n', flush=True)
    now_time()


if __name__ == '__main__':
    main()
