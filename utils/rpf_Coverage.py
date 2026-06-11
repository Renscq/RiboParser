#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_coverage.py


from utils.ribo import Coverage

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def coverage_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the coverage meta plot.")

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
    parser.add_argument('-b', dest="bin", required=False, type=str, default='30,100,30',
                        help="adjust the transcript to specified bins. 30 for 5'-UTR"
                             "and 3'-UTR, 100 for CDS. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads. (default: %(default)s).")
    parser.add_argument('--outlier', dest="outlier", action='store_true', required=False, default=False,
                        help="filter the outliers (default: %(default)s).")
    parser.add_argument('--set', dest="set", choices=['intersect', 'union'], required=False, type=str, default='union',
                        help="filter the gene list with 5-UTR / CDS / 3-UTR. (default: %(default)s).")
    parser.add_argument('--heat', dest="heatmap", action='store_true', required=False, default=False,
                        help="draw the coverage heatmap of whole gene. (default: %(default)s).")
    parser.add_argument('--bar', dest="barplot", action='store_true', required=False, default=False,
                        help="draw the coverage barplot of whole gene. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():
    now_time()
    print('\nDraw the metagene coverage.', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = coverage_args_parser()
    meta = Coverage.Coverage(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    meta.read_rpf()
    meta.import_gene()

    print('\nStep3: Adjust mRNAs to the same dimension.', flush=True)
    meta.process_utr5()
    meta.process_cds()
    meta.process_utr3()

    print('\nStep5: Draw the line plot of metagene coverage.', flush=True)
    meta.draw_meta_gene_line()

    print('\nStep6: Draw the heatmap of metagene coverage.', flush=True)
    if args.heatmap:
        meta.draw_meta_gene_heat()

    print('\nStep7: Draw the barplot of metagene coverage.', flush=True)
    if args.barplot:
        meta.draw_meta_gene_bar()

    print('\nStep8: Output the metagene coverage.', flush=True)
    meta.output_meta_gene()

    print('\nAll done.\n', flush=True)
    now_time()


if __name__ == '__main__':
    main()
