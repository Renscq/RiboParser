#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_coverage.py


import pandas as pd
import argparse


def merge_coverage_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the coverage files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for coverage files (e.g., '*_mean_coverage.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _mean_coverage.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_coverage_files(coverage_list):
    '''
    @Message  : retrieve the coverage table.
    @Input    : list --> pattern for coverage files (e.g., '*_mean_coverage.txt')
    @Return   : 
                output --> output dataframe contain the coverage
    @Flow     : step1 --> retrieve the coverage
    '''

    # make the dataframe
    coverage_colunms = ['Sample', 'Region', 'Bins', 'Density']
    coverage_df = pd.DataFrame(columns=coverage_colunms)

    # for each gne file
    for coverage_file in coverage_list:

        # import the coverage file
        coverage = pd.read_csv(coverage_file, sep='\t', index_col=False)

        # merge the gene coverage file
        coverage_df = pd.concat([coverage_df, coverage], axis=0)

    return coverage_df


def output_table(coverage_df, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the coverage
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    coverage_df.to_csv(output_file + '_mean_coverage.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_coverage_args_parser()

    print('Step2: import the coverage file.', flush=True)
    coverage_df = import_coverage_files(args.list)

    print('Step3: output the coverage table.', flush=True)
    output_table(coverage_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
