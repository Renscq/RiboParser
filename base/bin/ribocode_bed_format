#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : ribocode_bed_format.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/10/27 15:02:11
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import pandas as pd
# import polars as pl
# import numpy as np
from collections import OrderedDict
# from Bio import SeqIO
import argparse
import copy


def get_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to convert the ribocode bed file to genepred format.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        '-i', dest='input', required=True, type=str, help='input file name.'
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output genepred file name.'
    )

    # optional arguments
    parser.add_argument(
        '-t', dest='type', required=False, type=str, choices=['all', 'smorf', 'annotated'], default='all', 
        help='select the ORF class. (default: %(default)s)'
    )
    parser.add_argument(
        '-l', dest='length', required=False, type=int, default=8, help='length of the amino acids. (default: %(default)s)'
    )
    parser.add_argument(
        '-c', dest='clean', required=False, action='store_true', default=False, help='drop all the fuzzy ORFs. (default: %(default)s)'
    )
    args = parser.parse_args()

    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def filter_the_orf_class(ribocode, args):
    '''
    @Message  : filter the ORFs with different orf class.
    @Input    : ribocode --> filtered pvalue and length of ribocode results
    @Return   : ribocode --> ribocode results
    @Flow     : step1 --> filter the smORF and annotated genes

    '''
    
    tmp_ribocode = copy.deepcopy(ribocode)

    # filter the ORFType
    if args.type == 'smorf':
        # filter the smORF genes
        for index, orf in tmp_ribocode.items():
            if orf[-1] not in ['uORF', 'dORF', 'novel', 'dORF,novel', 'uORF,novel']:
                del ribocode[index]
        
    elif args.type == 'annotated':
        # filter the ORFType contain the string Annotated
        for index, orf in tmp_ribocode.items():
            if orf[-1] not in ['internal', 'annotated']:
                del ribocode[index]

    else:
        pass
    
    return ribocode


def flt_the_orf_length(ribocode, args):
    '''
    @Message  : unique the ORFs with the same gene-blocks.
    @Input    : ribocode --> filtered ribocode results
    @Return   : ribocode --> unique ribocode results
    @Flow     : step1 --> sort the ORFs with genome position and length of amino acids
                step2 --> drop the ORFs with the same gene-blocks
    '''
    
    tmp_ribocode = copy.deepcopy(ribocode)

    for index, orf in tmp_ribocode.items():
        orf_length = int(index.split('_')[-1])

        if orf_length < args.length:
            del ribocode[index]

    return ribocode


def fix_stop_codon(ribocode):
    '''
    @Message  : fix the stop codon of the ORFs.
    @Input    : ribocode --> filtered ribocode results
    @Return   : ribocode --> fixed stop codon of ORFs
    @Flow     : step1 --> fix the stop codon of the ORFs
    '''

    for index, orf in ribocode.items():
        if orf[2] == '+':
            orf[9][-1] = str(int(orf[9][-1]) + 3)

        else:
            orf[8][0] = str(int(orf[8][0]) - 3)

        orf[3] = orf[8][0]
        orf[4] = orf[9][-1]
        orf[5] = orf[8][0]
        orf[6] = orf[9][-1]

    return ribocode


def input_the_ribocode(args):
    '''
    @Message  : filter the pvalue and length of ORFs.
    @Input    : args --> input arguments
    @Return   : ribocode --> the ribocode result
    @Flow     : step1 --> input the ribocode result
                step2 --> drop the fuzzy ORFs
                step3 --> filter the length of ORFs
                step4 --> filter the ORFs with different orf class
                step5 --> output the ribocode result
    '''

    # input the ribocode result
    ribocode = OrderedDict()

    with open(args.input, 'r') as ribocode_input:
        for line in ribocode_input:
            if line.startswith('#'):
                continue
            
            record = line.strip().split('\t')
            
            # orf position
            chrom = record[0]

            # orf name
            orf_mess = record[3].split(';')

            orf_name = orf_mess[0]
            gene = orf_mess[1]
            rna = orf_mess[2]
            orf_type = orf_mess[3]

            # orf position
            txStart = int(orf_name.split('_')[-3])
            txEnd = int(orf_name.split('_')[-2])

            # orf strand
            score = record[4]
            strand = record[5]

            # fix the stop codon to the CDS
            if strand == '-':
                txStart, txEnd = txEnd - 1, txStart
            else:
                txStart, txEnd = txStart - 1, txEnd

            cdsStart = txStart
            cdsEnd = txEnd

            # fix the edge of the orf
            exonStarts = [record[1]]
            exonEnds = [record[2]]

            # store the orf information
            if orf_name not in ribocode.keys():
                exonCount = 1
                name = orf_type + '_' + orf_name

                ribocode[orf_name] = [name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, orf_type]
                
            else:
                exonCount = ribocode[orf_name][7] + 1
                if strand == '-':
                    exonStarts = exonStarts + ribocode[orf_name][8]
                    exonEnds = exonEnds + ribocode[orf_name][9]
                else:
                    exonStarts = ribocode[orf_name][8] + exonStarts
                    exonEnds = ribocode[orf_name][9] + exonEnds

                ribocode[orf_name] = [name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, orf_type]

    ribocode = filter_the_orf_class(ribocode, args)

    ribocode = flt_the_orf_length(ribocode, args)

    ribocode = fix_stop_codon(ribocode)

    return ribocode


def output_genepred(ribocode, args):
    '''
    @Message  : format the gene blocks to genepred.
    @Input    : ribocode --> filtered ribocode results
    @Return   : gene-blocks --> formatted gene-blocks
    @Flow     : step1 --> convert the gene block list to string
                step2 --> convert dict to dataframe, and sort the chrom and start position
                step3 --> output the gene-pred file
    '''
    
    # convet list to string
    for index, orf in ribocode.items():
        orf[8] = ','.join(orf[8]) + ','
        orf[9] = ','.join(orf[9]) + ','

    # convert dict to dataframe
    ribocode = pd.DataFrame.from_dict(ribocode, orient='index')
    ribocode[3] = ribocode[3].astype(int)

    # sort the chrom and start position
    ribocode = ribocode.sort_values(by=[1, 3], ascending=[True, True])

    # delete the last column of the ribocode 
    ribocode = ribocode.iloc[:, :-1]

    # output the genepred file
    ribocode.to_csv(args.output + '.genepred', sep='\t', header=False, index=False)


def main():
    print('Convert ribocode bed file to genepred format.', flush=True)
    args = get_args_parser()

    print('Step1: input the ribocode.', flush=True)
    record = input_the_ribocode(args)

    print('Step2: output the genepred file.', flush=True)
    output_genepred(record, args)

    print('Done.', flush=True)


if __name__ == '__main__':
    main()
