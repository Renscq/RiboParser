#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : fq_trim.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/11/07 17:47:32
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


# import pandas as pd
# import polars as pl
# import numpy as np
# from collections import OrderedDict
# from Bio import SeqIO
import argparse


def fq_trim_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to trim the fastq file.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        '-i', dest='input', required=True, type=str, help='input fq file name.'
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='output fq file name.'
    )

    # optional arguments
    parser.add_argument(
        '-l', dest='length', required=False, type=int, default=30, 
        help='nucleotide length (default: %(default)s nt).'
    )
    parser.add_argument(
        '-d', dest='discard', required=False, action='store_true', default=False, 
        help='discard the reads less than minimum length (default: %(default)s).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def trim_fq(args):
    '''
    @Message  : function for trim the fastq file.
    @Input    : fastq --> input file
                length --> length of nucleotide
    @Return   : output --> trimmed fastq file
    '''
    
    with open(args.input, 'r') as fq_in:
        with open(args.output, 'w') as fq_out:

            for line in fq_in:
                read_id = line.strip()
                read_seq = line.next().strip()
                read_name = line.next().strip()
                read_score = line.next().strip()

                if len(read_seq) >= args.length:
                    fq_out.writelines(read_id + '\n')
                    fq_out.writelines(read_seq[0:args.length] + '\n')
                    fq_out.writelines(read_name + '\n')
                    fq_out.writelines(read_score[0:args.length] + '\n')
                else:
                    pass


def main():
    args = fq_trim_args_parser()
    trim_fq(args)


if __name__ == '__main__':
    main()
