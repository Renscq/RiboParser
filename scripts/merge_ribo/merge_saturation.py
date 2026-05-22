#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_saturation.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_saturation_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the saturation files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for saturation files (e.g., '*_gene_saturation.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name.'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def process_saturation_files(gene_count_list):
    '''
    @Message  : retrieve the length table.
    @Input    : list --> pattern for gene saturation files (e.g., '*_gene_saturation.txt')
    @Return   : 
                output --> output nested dict contain the gene saturation
    @Flow     : step1 --> retrieve the gene saturation
    '''

    # make the dataframe with title
    gene_count_df = pd.DataFrame(columns=['Sample', 'Tercile', 'Covered', 'Uncovered'])

    # for each gne file
    for gene_file in gene_count_list:
        file_prefix = os.path.basename(gene_file).split('.')[0].replace('_gene_saturation', '')

        # import the gene saturation file
        gene_df = pd.read_csv(gene_file, sep='\t', index_col=False)
        gene_df.columns = ['Tercile', 'Covered']
        gene_df['Uncovered'] = gene_df['Covered'].max() - gene_df['Covered']
        gene_df['Sample'] = file_prefix
        gene_df = gene_df[['Sample', 'Tercile', 'Covered', 'Uncovered']]

        # merge the gene saturation file
        gene_count_df = pd.concat([gene_count_df, gene_df], axis=0)

    return gene_count_df


def output_table(gene_count_df, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the gene count
    @Return   : output --> description
    @Flow     : step1 --> convert the number to int
    '''

    # keep the 2 decimal places
    gene_count_df[['Tercile', 'Covered', 'Uncovered']] = gene_count_df[['Tercile', 'Covered', 'Uncovered']].astype('int')

    gene_count_df.to_csv(output_file + '_gene_saturation.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_saturation_args_parser()

    print('Step2: import the gene saturation file.', flush=True)
    gene_count_df = process_saturation_files(args.list)

    print('Step3: output the gene saturation table.', flush=True)
    output_table(gene_count_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
