#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_quant.py


import pandas as pd
import argparse


def merge_quant_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the quant files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for quant files (e.g., '*_rpf_quant.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _rpf_quant.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_quant_files(quant_list):
    '''
    @Message  : retrieve the quant table.
    @Input    : list --> pattern for quant files (e.g., '*_rpf_quant.txt')
    @Return   : 
                output --> output dataframe contain the quant
    @Flow     : step1 --> retrieve the quant
    '''

    quant_df = None
    # for each gne file
    for quant_file in quant_list:

        # import the quant file
        quant = pd.read_csv(quant_file, sep='\t', index_col=False)

        # merge the gene quant file
        if quant_df is None:
            quant_df = quant
        else:
            quant_df = pd.merge(quant_df, quant, on='name', how='outer')

    return quant_df


def output_table(quant_df, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the quant
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    quant_df.to_csv(output_file + '_rpf_quant.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_quant_args_parser()

    print('Step2: import the quant file.', flush=True)
    quant_df = import_quant_files(args.list)

    print('Step3: output the quant table.', flush=True)
    output_table(quant_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
