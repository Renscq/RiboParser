#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : RiboParser
# @Script  : rna_Offset.py


import pandas as pd
import numpy as np

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def rna_offset_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to create the offset table of RNA-seq.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _offset.txt)")

    # arguments for the offset table creation
    parser.add_argument('-m', dest="min", required=False, type=int, default=25,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=151,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-e', dest="exp_offset", required=False, type=int, default=12,
                        help="Expected offset (default: %(default)s nt).")

    args = parser.parse_args()
    args_print(args)

    return args


def main():
    now_time()
    print('\nCreate the p-site offset for RNA-seq.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = rna_offset_args_parser()
    
    # create the Offset table
    min_length = args.min
    max_length = args.max
    length_range = np.arange(min_length, max_length + 1)
    offset_table = pd.DataFrame(length_range, index=length_range, columns=['length'])

    offset_table['frame0'] = args.exp_offset
    offset_table['rpfs0'] = 0

    offset_table['frame1'] = args.exp_offset + 1
    offset_table['rpfs1'] = 0
    
    offset_table['frame2'] = args.exp_offset + 2
    offset_table['rpfs2'] = 0

    offset_table['rpfs'] = 0

    offset_table['p_site'] = args.exp_offset + 1
    offset_table['periodicity'] = 100
    offset_table['ribo'] = 'first'

    # output the offset table
    print('\nStep2: Output the offset table.', flush=True)
    offset_table.to_csv(args.output + '_offset.txt', sep='\t', index=False)
    
    now_time()
    print('\nAll done.', flush=True)


if __name__ == '__main__':
    main()
