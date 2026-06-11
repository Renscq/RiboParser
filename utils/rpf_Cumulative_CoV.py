#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Cumulative_CoV.py


from utils.ribo.Cumulative_CoV import *

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def cumulative_cov_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to calculate the cumulative CoV.")

    # arguments for the required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs density file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, 
                             help="prefix of output file name (default: filename + '_cumulative_CoV.txt'.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="list", required=False, type=str,
                        help="the list of input genes for transcript id.")
    parser.add_argument('-m', dest="min", required=False, type=int, default=0,
                        help="retain transcript with more than minimum RPFs (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM (default: %(default)s).")
    
    parser.add_argument('-t', dest="trim", required=False, type=int, default=50,
                        help="trim transcript with specific length (default: %(default)s nt).")
    
    parser.add_argument('-s', dest="split", action='store_true', required=False, default=False,
                        help="split gene rpf to each TXT file (default: %(default)s).")
    parser.add_argument('-z', dest="zero", action='store_true', required=False, default=False,
                        help="set the start site to zero (default: %(default)s).")
    
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():
    now_time()
    print('\nRetrieve the RPFs with gene list.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = cumulative_cov_args_parser()
    rpfs = CumulativeCoV(args)

    print('\nStep2: Retrieve the gene RPFs.', flush=True)
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    print('\nStep3: Format the RPFs table.', flush=True)
    rpfs.melt_rpf_table()

    print('\nStep4: Calculate the cumulative CoV.', flush=True)
    rpfs.calc_cov()

    print('\nStep5: Output the cumulative CoV meta table.', flush=True)
    rpfs.merge_cov_table()

    print('\nStep5: Output the cumulative CoV table.', flush=True)
    rpfs.output_rpf_table()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
