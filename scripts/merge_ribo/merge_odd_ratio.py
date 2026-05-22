#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_odd_ratio.py


import pandas as pd
import argparse


def merge_odd_ratio_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the codon odd ratio files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for codon odd ratio files (e.g., '*_sum_codon_odd_ratio.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _codon_odd_ratio.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_odd_ratio_files(odd_ratio_list):
    '''
    @Message  : retrieve the codon odd ratio table.
    @Input    : list --> pattern for codon odd ratio files (e.g., '*_sum_codon_odd_ratio.txt')
    @Return   : 
                output --> output dataframe contain the odd ratio
    @Flow     : step1 --> retrieve the codon odd ratio
    '''

    # make the dataframe
    odd_ratio_df = pd.DataFrame()

    # for each odd_ratio file
    for odd_ratio_file in odd_ratio_list:
        # import the odd_ratio file
        odd_ratio = pd.read_csv(odd_ratio_file, sep='\t', index_col=False)

        # merge the gene odd_ratio file
        if odd_ratio_df.empty:
            odd_ratio_df = odd_ratio
        else:
            # delete the columns ['AA', 'Abbr.']
            odd_ratio_df = pd.concat([odd_ratio_df, odd_ratio], axis=0)

    return odd_ratio_df


def output_table(odd_ratio_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the codon odd ratio
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    odd_ratio_merge.to_csv(output_file + '_sum_codon_odd_ratio.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_odd_ratio_args_parser()

    print('Step2: import the codon odd ratio file.', flush=True)
    odd_ratio_df = import_odd_ratio_files(args.list)

    print('Step3: output the codon odd ratio table.', flush=True)
    output_table(odd_ratio_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
