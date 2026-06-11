#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_odd_ratio.py


from utils.ribo import Odd_Ratio

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def rpf_odd_ratio_args_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to calculate the odd ratio of stalling codon.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the odd ratio calculation
    parser.add_argument('--thread', dest="thread", required=False, type=int, default=1,
                        help="set the thread number. (default: %(default)s).")
    
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for odd ratio calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for odd ratio calculation. (default: %(default)s).")
    parser.add_argument('-c', dest="control", required=True, type=str,
                        help="specify the name of control samples, separated by commas.")
    parser.add_argument('-t', dest="treat", required=True, type=str,
                        help="specify the name of treat samples, separated by commas.")
    
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('-z', dest="zero", action='store_true', required=False, default=False,
                        help="fill the empty codon with [0.01 * minimal value]. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="specify the minimum count of RPFs per transcript to keep. (default: %(default)s).")
    parser.add_argument('--fdr', dest="fdr", required=False, type=str, choices=['bhfdr', 'pvalue'], default='bhfdr',
                        help="specify the multipletests. (default: %(default)s).")
    parser.add_argument('-v', dest="value", required=False, type=float, default=0.05,
                        help="specify the p value. (default: %(default)s).")
    
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="specify the number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="specify the number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--stop', dest="stop", action='store_true', required=False, default=False,
                        help="remove the stop codon from table. (default: %(default)s).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the codon odd ratio. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():
    now_time()
    print('\nCalculate the relative codon pausing score.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = rpf_odd_ratio_args_parser()
    odd_ratio = Odd_Ratio.OddRatio(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    odd_ratio.read_rpf()

    print('\nStep3: Make two dimensional table for fisher exactly test.', flush=True)
    odd_ratio.make_two_dimensional_table()

    print('\nStep4: Calculate the odd ratio.', flush=True)
    odd_ratio.calc_odd_ratio()
    odd_ratio.calc_chi2_test2()

    print('\nStep5: Output the ribosome pausing odd ratio.', flush=True)
    odd_ratio.output_odd_ratio()
    odd_ratio.summarize_odd_ratio()
    
    print('\nStep6: Draw the ribosome pausing odd ratio.', flush=True)
    odd_ratio.draw_odd_ratio_line()
    odd_ratio.draw_odd_ratio_scatter()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
