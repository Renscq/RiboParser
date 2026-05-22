#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_cdt.py


import pandas as pd
import argparse


def merge_cdt_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the codon decoding time files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for cdt files (e.g., '*_cdt.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _cdt.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_cdt_files(cdt_list):
    '''
    @Message  : retrieve the cdt table.
    @Input    : list --> pattern for cdt files (e.g., '*_codon_cdt.txt')
    @Return   : 
                output --> output dataframe contain the cdt
    @Flow     : step1 --> retrieve the cdt
    '''

    # make the dataframe
    cdt_df = pd.DataFrame()

    # for each cdt file
    for cdt_file in cdt_list:
        # import the cdt file
        cdt = pd.read_csv(cdt_file, sep='\t', index_col=False)
        cdt.sort_values(by=['Abbr', 'Codon'], ascending=True, inplace=True, ignore_index=True)

        # merge the gene cdt file
        if cdt_df.empty:
            cdt_df = cdt
        else:
            # delete the columns ['AA', 'Abbr.']
            cdt = cdt.drop(['AA', 'Abbr'], axis=1)
            cdt_df = pd.merge(left=cdt_df, right=cdt, how='outer', on='Codon')

    # reorder the cdt columns with 
    # ['codon', 'AA', 'Abbr.', *_density, *_absolute_cdt, *_relative_cdt]
    cdt_codon = cdt_df[['Codon', 'AA', 'Abbr']]
    absolute_cdt = cdt_df.filter(regex='_absolute_cdt')
    relative_cdt = cdt_df.filter(regex='_relative_cdt').filter(regex='^(?!.*norm).*$')
    norm_cdt = cdt_df.filter(regex='_norm_cdt')
    relative_norm_cdt = cdt_df.filter(regex='_norm_relative_cdt')

    cdt_merge = pd.concat([cdt_codon, absolute_cdt, norm_cdt, relative_cdt, relative_norm_cdt], axis=1)

    return cdt_merge


def output_table(cdt_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the cdt
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    cdt_merge.to_csv(output_file + '_codon_decoding_time.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_cdt_args_parser()

    print('Step2: import the cdt file.', flush=True)
    cdt_df = import_cdt_files(args.list)

    print('Step3: output the cdt table.', flush=True)
    output_table(cdt_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
