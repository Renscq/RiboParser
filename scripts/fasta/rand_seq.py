#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : rand_seq.py


# import pandas as pd
# import polars as pl
import numpy as np
# from collections import OrderedDict
# from Bio import SeqIO
import argparse


def rand_seq_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to generate random sequence.')

    # needed arguments
    seq_class = parser.add_mutually_exclusive_group(required=True)
    seq_class.add_argument(
        '-d', dest='dna', required=False, action='store_true', default=False, 
        help='generate dna nucleotide sequence.'
    )
    seq_class.add_argument(
        '-r', dest='rna', required=False, action='store_true', default=False,
        help='generate rna nucleotide sequence.'
    )
    seq_class.add_argument(
        '-p', dest='protein', required=False, action='store_true', default=False,
        help='generate protein sequence.'
    )

    # optional arguments
    parser.add_argument(
        '-l', dest='length', required=False, type=int, default=10, 
        help='length of sequeuce. (default: %(default)s nt).'
    )
    parser.add_argument(
        '-n', dest='number', required=False, type=int, default=5, 
        help='number of sequence. (default: %(default)s).'
    )
    parser.add_argument(
        '-f', dest='format', required=False, choices=['fa', 'txt'], type=str, default='fa', 
        help='format of sequence. (default: %(default)s).'
    )

    args = parser.parse_args()
    # args_dict = vars(args)
    # for k, v in args_dict.items():
    #     print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def rand_seq(args):
    '''
    @Message  : function for generate the random sequence.
    @Input    : seq class --> dna/rna/protein
                seq length --> length of sequence
                seq number --> number of sequence
                seq format --> format of sequence

    @Return   : output --> random sequences
    @Flow     : step1 --> check the input arguments
                step2 --> generate the random sequence
                step3 --> output the random sequence
    '''
    
    # set the random seed
    # np.random.seed(1314)

    # check the sequence class
    if args.dna:
        seq_class = 'd'
        seq_list = ['A', 'T', 'C', 'G']

    elif args.rna:
        seq_class = 'r'
        seq_list = ['A', 'U', 'C', 'G']

    else:
        seq_class = 'p'
        seq_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    for i in range(args.number):
        seq = ''.join(np.random.choice(seq_list, args.length))
        if args.format == 'fa':
            print('>{}{}'.format(seq_class, i+1))
            print(seq)
        else:
            print(seq)


def main():
    args = rand_seq_args_parser()
    rand_seq(args)


if __name__ == '__main__':
    main()
