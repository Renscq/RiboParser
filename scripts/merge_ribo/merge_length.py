#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_length.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_length_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the length distribution files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for length distribution files (e.g., '*length_distribution.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name.'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def process_log_files(length_file_list):
    '''
    @Message  : retrieve the length table.
    @Input    : list --> pattern for reads length files (e.g., '*length_distribution.txt')
    @Return   : 
                output --> output nested dict contain the reads length
    @Flow     : step1 --> retrieve the reads length
    '''

    # make the dataframe
    reads_length_df = pd.DataFrame()
    reads_length_df['Sample'] = []
    reads_length_df['Length'] = []
    reads_length_df['Plus_Count'] = []
    reads_length_df['Minus_Count'] = []
    reads_length_df['Plus_Ratio'] = []
    reads_length_df['Minus_Ratio'] = []

    # for each reads length file
    for length_file in length_file_list:
        file_prefix = os.path.basename(length_file).split('.')[0].replace('_length_distribution', '')

        # import the reads length file
        length_df = pd.read_csv(length_file, sep='\t', index_col=False)

        # rename the columns
        length_df.columns = ['Length', 'Plus_Count', 'Minus_Count']

        # calculate the ratio of reads length
        length_df['Plus_Ratio'] = length_df['Plus_Count'] / length_df['Plus_Count'].sum() * 100
        length_df['Minus_Ratio'] = length_df['Minus_Count'] / length_df['Minus_Count'].sum() * 100

        length_df['Sample'] = file_prefix

        # merge the reads length file
        reads_length_df = pd.concat([reads_length_df, length_df], axis=0)

    return reads_length_df


def output_table(reads_length_df, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the mapped reads
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    # keep the 2 decimal places
    reads_length_df[['Length', 'Plus_Count', 'Minus_Count']] = reads_length_df[['Length', 'Plus_Count', 'Minus_Count']].astype('int')
    reads_length_df[['Plus_Ratio', 'Minus_Ratio']] = reads_length_df[['Plus_Ratio', 'Minus_Ratio']].astype('float').round(4)

    reads_length_df.to_csv(output_file + '_length_distribution.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_length_args_parser()

    print('Step2: import the reads length file.', flush=True)
    reads_length_df = process_log_files(args.list)

    print('Step3: output the length distribution table.', flush=True)
    output_table(reads_length_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
