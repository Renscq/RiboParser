#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      :   riboparser
@Script       :   fa_length.py
@Environment  :   python 3.8.5
@Version      :   1.0
@Author       :   Rensc 
@Time         :   2023/05/29 11:26:28
@E-mail       :   rensc0718@163.com
@License      :   (C)Copyright 2023-2025, Ren Shuchao
'''


import argparse
from Bio import SeqIO
import gzip
import os


def args_parser():

    parser = argparse.ArgumentParser(description='This script is used to summarize the reads length.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-i', dest='input', required=True, type=str, help='input file name')
    input_group.add_argument('-o', dest='output', required=True, type=str, help='output file name')

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


def summarize_reads(record):
    '''
    @Message  :   summarize the reads length
    @Input    :   SeqIO record
    @Return   :   a dict contain length and number
    
    {35: 103, 31: 904, 33: 38391, 32: 391495}
    '''
    
    length_dict = {}
    line = 0
    
    for reads in SeqIO.parse(record, "fasta"):
        line += 1
        if line % 1000000 == 0:
            print("rows : " + str(line), flush=True)

        length = len(reads)
        try:
            length_dict[length] += 1
        except KeyError:
            length_dict[length] = 1
    
    print("rows : " + str(line) + '\n', flush=True)

    return length_dict


def output_length_distr(length_dict, out_file):
    '''
    @Message  :   sort the length dict and output data
    @Input    :   dict contain length and number
    @Return   :   output the txt file
    '''
    length_dict = dict(sorted(length_dict.items(), key=lambda x: x[0]))

    with open(out_file, 'w') as out_flag:
        out_flag.writelines(''.join(['length', '\t', 'count', '\n']))

        for length, count in length_dict.items():
            out_flag.writelines(''.join([str(length), '\t', str(count), '\n']))


def main():
    print("Summarize the reads length of fasta file.", flush=True)
    args = args_parser()

    print("\nStep1: import the fasta file.", flush=True)
    record = import_reads(args.input)

    print("Step2: summarize the reads length.", flush=True)
    length_dict = summarize_reads(record)

    print("Step3: output the length distribution", flush=True)
    output_length_distr(length_dict, args.output)



if __name__ == '__main__':

    main()
