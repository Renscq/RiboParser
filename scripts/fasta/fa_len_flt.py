#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : fa_len_flt.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/11/29 15:44:09
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''

import argparse
import os
import gzip
from Bio import SeqIO


# get arguments and define the scripts usage
def args_parser():
    parser = argparse.ArgumentParser(description='This script is used to filter the reads length.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-i', dest='input', required=True, type=str, help='input fasta file name')
    input_group.add_argument('-o', dest='output', required=True, type=str, help='output fasta file name')

    # optional arguments
    optional_group = parser.add_argument_group('Optional arguments')
    optional_group.add_argument('-m', dest='min', type=int, default=8, help='min length of reads. (default %(default)s)')
    optional_group.add_argument('-M', dest='max', type=int, default=1000, help='max length of reads. (default %(default)s)')

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print("{:<12}:  {:<}".format(k, str(v)), flush=True)

    return args


def import_reads(fa_file):
    '''
    @Message  :   import the fasta file
    @Input    :   fasta filename
    @Return   :   seqio record
    '''
    
    fa_suffix = os.path.splitext(fa_file)[1]
    if fa_suffix == '.gz':
        record = gzip.open(fa_file, "rt")
    else:
        record = open(fa_file, "r")

    return record


def filter_reads(record, output, min_length, max_length):
    '''
    @Message  :   filter the reads length
    @Input    :   SeqIO record
    @Return   :   output the fasta file contain reads length between min and max
    
    '''

    with open(output, "w") as out_flag:
        line = 0
        for reads in SeqIO.parse(record, "fasta"):
            line += 1
            if line % 1000000 == 0:
                print("rows : " + str(line), flush=True)

            length = len(reads)
            if length >= min_length and length <= max_length:
                SeqIO.write(reads, out_flag, "fasta")


def main():
    print("Filter the reads length of fasta file.", flush=True)
    args = args_parser()

    print("\nStep1: import the fasta file.", flush=True)
    record = import_reads(args.input)

    print("Step2: output the fasta.", flush=True)
    filter_reads(record, args.output, args.min, args.max)


if __name__ == "__main__":
    main()
