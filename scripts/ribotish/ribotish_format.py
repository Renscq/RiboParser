#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : ribotish_to_gtf.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/10/26 17:40:16
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import pandas as pd
# import polars as pl
# import numpy as np
from collections import OrderedDict
# from Bio import SeqIO
import argparse

def get_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to convert the ribotish result to genepred format.')

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
        '-p', dest='pvalue', required=False, type=float, default=0.05, help='threshold of the ribo pvalue. (default: %(default)s)'
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


def filter_the_pvalue(ribotish, args):
    '''
    @Message  : filter the pvalue and length of ORFs.
    @Input    : ribotish --> the ribotish result
    @Return   : ribotish --> filtered pvalue and length of ribotish results
    @Flow     : step1 --> filter the pvalue
                step2 --> filter the length of ORFs
    '''
    
    # filter the pvalue
    ribotish = ribotish.loc[ribotish['RiboPvalue'] <= args.pvalue, :]
    ribotish = ribotish.loc[ribotish['AALen'] >= args.length, :]

    return  ribotish


def filter_the_orf_class(ribotish, args):
    '''
    @Message  : filter the ORFs with different orf class.
    @Input    : ribotish --> filtered pvalue and length of ribotish results
    @Return   : ribotish --> ribotish results
    @Flow     : step1 --> filter the smORF and annotated genes
    '''
    
    # replace the ['5UTR', '3UTR', 'Novel'] string in the part of TisType to ['uORF', 'dORF', 'Novel']
    ribotish['TisType'] = ribotish['TisType'].str.replace("5'UTR", 'uORF')
    ribotish['TisType'] = ribotish['TisType'].str.replace("3'UTR", 'dORF')

    # clean the TisType
    if args.clean:
        # remove the TisType contain the string 'Known'
        ribotish = ribotish.loc[~ribotish['TisType'].str.contains('Known'), :]

        # remove the TisType contain the string 'CDSFrameOverlap'
        ribotish = ribotish.loc[~ribotish['TisType'].str.contains('CDSFrameOverlap'), :]
        
    else:
        pass

    # filter the smORF and annotated genes
    if args.type == 'smorf':
        # filter the TisType contain the string 5'UTR/3'UTR/Novel
        ribotish = ribotish.loc[ribotish['TisType'].str.contains("uORF|dORF|Novel"), :]

    elif args.type == 'annotated':
        # filter the TisType contain the string Annotated
        ribotish = ribotish.loc[ribotish['TisType'].str.contains("Annotated|Extended|Internal|Truncated"), :]

    else:
        pass

    return ribotish


def unique_the_orf(ribotish):
    '''
    @Message  : unique the ORFs with the same gene-blocks.
    @Input    : ribotish --> filtered ribotish results
    @Return   : ribotish --> unique ribotish results
    @Flow     : step1 --> sort the ORFs with genome position and length of amino acids
                step2 --> drop the ORFs with the same gene-blocks
    '''
    
    ribotish = ribotish.sort_values(by=['GenomePos', 'AALen'], ascending=[True, False])

    ribotish = ribotish.drop_duplicates(subset=['TisType', 'Blocks'], keep='first')

    return ribotish


def input_the_ribotish(args):
    '''
    @Message  : filter the pvalue and length of ORFs.
    @Input    : args --> input arguments
    @Return   : ribotish --> the ribotish result
    @Flow     : step1 --> input the ribotish result
                step2 --> drop the fuzzy ORFs
                step3 --> filter the pvalue and length of ORFs
                step4 --> filter the ORFs with different orf class
                step5 --> unique the ORFs with the same blocks
                step6 --> output the ribotish result
    '''

    ribotish = pd.read_csv(args.input, sep='\t', header=0, names=None)

    ribotish = filter_the_pvalue(ribotish, args)

    ribotish = filter_the_orf_class(ribotish, args)

    ribotish = unique_the_orf(ribotish)

    return ribotish


def retrieve_aa(ribotish, args):
    '''
    @Message  : retrieve the amino acid of the ORFs.
    @Input    : ribotish --> filtered ribotish results
    @Return   : aa_sequence --> the amino acid sequence of the ORFs
    @Flow     : step1 --> retrieve the amino acid sequence
                step2 --> format the amino acid sequence to fasta format
                step3 --> output the amino acid sequence
    '''

    with open(args.output + '.aa', 'w') as aa_output:
        for index, row in ribotish.iterrows():

            aa_name = row['TisType'] + '_' + row['Gid'] + '_' + row['Tid'] + '_' + row['GenomePos'] + ' ' + row['Blocks']
            aa_seq = row['AASeq'][:-1]

            aa_output.writelines('>' + aa_name + '\n')
            aa_output.writelines(aa_seq + '\n')
    

def retrieve_cds(ribotish, args):
    '''
    @Message  : retrieve the amino acid of the ORFs.
    @Input    : ribotish --> filtered ribotish results
    @Return   : aa_sequence --> the amino acid sequence of the ORFs
    @Flow     : step1 --> retrieve the cds sequence
                step2 --> format the cds sequence to fasta format
                step3 --> output the cds sequence
    '''

    with open(args.output + '.cds', 'w') as cds_output:
        for index, row in ribotish.iterrows():
            aa_name = row['Gid'] + '_' + row['Tid'] + ' ' + row['GenomePos'] + ' ' + row['TisType'] + ' ' + row['Blocks']
            aa_seq = row['Seq']

            cds_output.writelines('>' + aa_name + '\n')
            cds_output.writelines(aa_seq + '\n')


def get_cds_frame(strand, gene_blocks):
    '''
    @Message  : get the cds frame.
    @Input    : strand --> gene strand
                gene_blocks --> gene blocks
    @Return   : cds_frame --> frame of cds
    @Flow     : step1 --> check the different strand
                step2 --> get the list of gene blocks
                step3 --> calculate the length of cds
                step4 --> calculate the frame of cds
                
    @notes    : genepred --> cds length = txEnd - txStart
                gtf --> cds length = txEnd - txStart + 1
                bed --> cds length = txEnd - txStart

    '''

    if strand == '+':
        orf_length = 0

        for index, block in enumerate(gene_blocks):
            t_start, t_end = block.split('-')

            if index == 0:
                cds_frame = [0]
                exonStarts = [int(t_start)]
                exonEnds = [int(t_end)]

                continue

            exonStarts.append(int(t_start))
            exonEnds.append(int(t_end))
            cds_len = int(t_end) - int(t_start)
            
            orf_length += cds_len
            frame = orf_length % 3

            cds_frame = cds_frame + [frame]

    elif strand == '-':
        orf_length = 0
        
        exonEnds = []

        for index, block in enumerate(gene_blocks):
            t_start, t_end = block.split('-')

            if index == 0:
                cds_frame = [0]
                exonStarts = [int(t_start)]
                exonEnds = [int(t_end)]

                continue

            exonStarts = [int(t_start)] + exonStarts
            exonEnds = [int(t_end)] + exonEnds
            cds_len = int(t_end) - int(t_start)
            
            orf_length += cds_len
            frame = orf_length % 3

            cds_frame = [frame] + cds_frame

    else:
        pass
    
    cds_frame = ','.join([str(item) for item in cds_frame]) + ','
    exonStarts = ','.join([str(item) for item in exonStarts]) + ','
    exonEnds = ','.join([str(item) for item in exonEnds]) + ','

    return cds_frame, exonStarts, exonEnds


def retrieve_blocks(ribotish, args):
    '''
    @Message  : format the gene blocks to genepred.
    @Input    : ribotish --> filtered ribotish results
    @Return   : gene-blocks --> formatted gene-blocks
    @Flow     : step1 --> retrieve the gene-blocks
                step2 --> format the gene-blocks
                step3 --> output the gene-blocks file
    '''
    
    with open(args.output + '.genepred', 'w') as genepred_output:

        for index, row in ribotish.iterrows():
            # gene names
            name = row['TisType'] + '_' + row['Gid'] + '_' + row['Tid']
            chrom = row['GenomePos'].split(':')[0]
            # source = 'ribotish'
            # type = 'CDS'

            # gene position
            strand = row['GenomePos'].split(':')[-1]

            gene_blocks = row['Blocks'].split(',')
            exonCount = len(gene_blocks)

            exonStarts_list = [item.split('-')[0] for item in gene_blocks]
            exonEnds_list = [item.split('-')[-1] for item in gene_blocks]

            txStart = exonStarts_list[0]
            txEnd = exonEnds_list[-1]
            cdsStart = txStart
            cdsEnd = txEnd

            exonStarts = ','.join(exonStarts_list) + ','
            exonEnds = ','.join(exonEnds_list) + ','

            # output the gene-blocks
            # gene_blocks = [name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, IstringexonFrames]
            gene_blocks = [name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds]

            genepred_output.writelines('\t'.join([str(item) for item in gene_blocks]) + '\n')


def main():
    print('Convert ribotish results to genepred format.', flush=True)
    args = get_args_parser()

    print('Step1: input the ribotish.', flush=True)
    record = input_the_ribotish(args)

    print('Step2: output the cds file.', flush=True)
    retrieve_aa(record, args)

    print('Step3: output the amino acid file.', flush=True)
    retrieve_cds(record, args)

    print('Step4: output the genepred file.', flush=True)
    retrieve_blocks(record, args)

    print('Done.', flush=True)


if __name__ == '__main__':
    main()
