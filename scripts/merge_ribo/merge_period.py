#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_period.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_period_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the 3-nt periodicity files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for periodicity files (e.g., '*_periodicity.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _periodicity.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_period_files(period_list):
    '''
    @Message  : retrieve the periodicity table.
    @Input    : list --> pattern for 3-nt periodicity files (e.g., '*_periodicity.txt')
    @Return   : 
                output --> output dataframe contain the 3-nt periodicity
    @Flow     : step1 --> retrieve the 3-nt periodicity
    '''

    # make the dataframe
    period_colunms = ['Sample', 'Frame', 'Count', 'Ratio']
    period_df = pd.DataFrame(columns=period_colunms)

    # for each gne file
    for period_file in period_list:
        # import the motif file
        period = pd.read_csv(period_file, sep='\t', index_col=False)

        # merge the gene periodion file
        period_df = pd.concat([period_df, period], axis=0)

    return period_df


def output_table(period_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the 3-nt periodicity
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    period_merge.to_csv(output_file + '_periodicity.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_period_args_parser()

    print('Step2: import the 3-nt periodicity file.', flush=True)
    period_df = import_period_files(args.list)

    print('Step3: output the 3-nt periodicity table.', flush=True)
    output_table(period_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
