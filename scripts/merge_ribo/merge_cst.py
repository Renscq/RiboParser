#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_cst.py


import pandas as pd
import argparse


def merge_cst_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the codon decoding time files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for cst files (e.g., '*_codon_selection_time.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _codon_selection_time.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_cst_files(cst_list):
    '''
    @Message  : retrieve the cst table.
    @Input    : list --> pattern for cst files (e.g., '*_codon_selection_time.txt')
    @Return   : 
                output --> output dataframe contain the cst
    @Flow     : step1 --> retrieve the cst
    '''

    # make the dataframe
    cst_df = pd.DataFrame()

    # for each cst file
    for cst_file in cst_list:
        # import the cst file
        cst = pd.read_csv(cst_file, sep='\t', index_col=False)
        cst.sort_values(by=['Abbr', 'Codon'], ascending=True, inplace=True, ignore_index=True)

        # merge the gene cst file
        if cst_df.empty:
            cst_df = cst
        else:
            # delete the columns ['AA', 'Abbr.']
            cst = cst.drop(['AA', 'Abbr'], axis=1)
            cst_df = pd.merge(left=cst_df, right=cst, how='outer', on='Codon')

    # reorder the cst columns with 
    # ['codon', 'AA', 'Abbr.', *_density, *_absolute_cst, *_relative_cst]
    cst_codon = cst_df[['Codon', 'AA', 'Abbr']]
    absolute_cst = cst_df.filter(regex='_absolute_cst')
    relative_cst = cst_df.filter(regex='_relative_cst')

    cst_merge = pd.concat([cst_codon, absolute_cst, relative_cst], axis=1)

    return cst_merge


def output_table(cst_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the cst
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    cst_merge.to_csv(output_file + '_codon_selection_time.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_cst_args_parser()

    print('Step2: import the cst file.', flush=True)
    cst_df = import_cst_files(args.list)

    print('Step3: output the cst table.', flush=True)
    output_table(cst_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
