#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : fa_gc_sum.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/11/07 19:46:39
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import pandas as pd
# import polars as pl
# import numpy as np
from collections import OrderedDict
from collections import Counter
# from Bio import SeqIO
import argparse


def fa_gc_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to calculate the GC content of fasta file.')

    require_group = parser.add_argument_group('Required arguments')
    require_group.add_argument(
        '-i', dest='input', required=True, type=str, 
        help='input fasta file name.'
    )
    require_group.add_argument(
        '-o', dest='output', required=True, type=str, 
        help='output txt file name.'
    )

    args = parser.parse_args()

    return args


def input_fa(args):
    '''
    @Message  : function for import the fasta file.
    @Input    : fasta --> input file
    @Return   : output --> fasta dict
    '''

    fa_dict = OrderedDict()
    seq_name = None

    with open(args.input, 'r') as fa_in:
        for line in fa_in:
            
            if line.startswith('>'):
                seq_name = line.strip()
                if seq_name in fa_dict:
                    fa_dict[seq_name + 'rep2'] = []
                else:
                    fa_dict[seq_name] = ['', 0, 0, 0, 0]
            else:
                fa_dict[seq_name][0] += line.strip()
    
    return fa_dict


def calc_gc(args, fa_dict):
    '''
    @Message  : function for calculate the GC content of fasta file.
    @Input    : fasta --> input file
    @Return   : output --> GC content of fasta file
    '''

    gc_dict = OrderedDict()

    for seq_name, seq_info in fa_dict.items():
        gc_dict[seq_name] = [seq_name,
                             seq_info[0].count('A'), 
                             seq_info[0].count('T'), 
                             seq_info[0].count('C'), 
                             seq_info[0].count('G')]

    gc_df = pd.DataFrame.from_dict(gc_dict, orient='index', columns=['name', 'A', 'T', 'C', 'G'])

    gc_df.to_csv(args.output, sep='\t', index=False)

    sum_gc = gc_df.iloc[:, 1:5].sum(axis=0)
    print(sum_gc)


def main():
    args = fa_gc_args_parser()
    fa_dict = input_fa(args)
    calc_gc(args, fa_dict)


if __name__ == '__main__':
    main()
