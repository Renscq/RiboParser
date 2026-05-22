#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_rsem.py


import pandas as pd
import argparse

def rsem_merge_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge specified columns from result files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for result files (e.g., '*results')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='output file name.'
    )

    # optional arguments
    parser.add_argument(
        "-c", "--column", required=False, choices=['expected_count', 'TPM', 'FPKM'], default='expected_count',
        help="Column name to merge (e.g., 'expected_count')."
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def merge_files(file_list, column_name):
    '''
    @Message  : merge specified columns from result files.
    @Input    : list --> pattern for result files (e.g., '*results')
                column --> column name to merge (e.g., 'expected_count')
    @Return   : 
                output --> output file name for merged results
    @Flow     : step1 --> run
    '''
    
    merged_data = []
    for file in file_list:
        data = pd.read_csv(file, sep='\t', index_col=0)
        selected_column = data[[column_name]]
        # remove the path and suffix
        filename = file.split('/')[-1]
        selected_column.columns = [filename]
        merged_data.append(selected_column)

    return merged_data


def output_file(merged_data, output_file):
    
    merged_data = pd.concat(merged_data, axis=1)
    merged_data.to_csv(output_file, sep='\t')


def main():
    args = rsem_merge_args_parser()
    merged_data = merge_files(args.list, args.column)
    output_file(merged_data, args.output)

if __name__ == "__main__":
    main()
