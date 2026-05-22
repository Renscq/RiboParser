#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_pausing.py


import pandas as pd
import argparse


def merge_pausing_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the pausing score files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for pausing score files (e.g., '*sum_codon_pausing_score.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + sum_codon_pausing_score.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_pausing_files(pausing_list):
    '''
    @Message  : retrieve the pausing table.
    @Input    : list --> pattern for pausing files (e.g., '*sum_codon_pausing_score.txt')
    @Return   : 
                output --> output dataframe contain the pausing
    @Flow     : step1 --> retrieve the pausing
    '''

    # make the dataframe
    pausing_df = pd.DataFrame()

    # for each pausing file
    for pausing_file in pausing_list:
        # import the pausing file
        pausing = pd.read_csv(pausing_file, sep='\t', index_col=False)
        pausing.sort_values(by=['Abbr', 'Codon'], ascending=True, inplace=True, ignore_index=True)

        # merge the gene pausing file
        if pausing_df.empty:
            pausing_df = pausing
        else:
            # delete the columns ['AA', 'Abbr.', 'count']
            pausing = pausing.drop(['AA', 'Abbr', 'Count'], axis=1)
            pausing_df = pd.merge(left=pausing_df, right=pausing, how='outer', on='Codon')

    # reorder the pausing columns with 
    # ['codon', 'AA', 'Abbr.', 'count', *_valid_codon, *_rpf_count, *_absolute_total_ps, *_absolute_valid_ps, *_relative_total_ps, *_relative_valid_ps]
    pausing_codon = pausing_df[['Codon', 'AA', 'Abbr', 'Count']]
    pausing_valid_codon = pausing_df.filter(regex='_valid_codon')
    pausing_rpf_count = pausing_df.filter(regex='_rpf_count')
    pausing_total_ps = pausing_df.filter(regex='_absolute_total_ps')
    pausing_valid_ps = pausing_df.filter(regex='_absolute_valid_ps')
    pausing_rel_total_ps = pausing_df.filter(regex='_relative_total_ps')
    pausing_rel_valid_ps = pausing_df.filter(regex='_relative_valid_ps')

    pausing_merge = pd.concat([pausing_codon, pausing_valid_codon, pausing_rpf_count,
                               pausing_total_ps, pausing_valid_ps, pausing_rel_total_ps, pausing_rel_valid_ps], 
                               axis=1)

    return pausing_merge


def output_table(pausing_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the pausing
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    pausing_merge.to_csv(output_file + '_sum_codon_pausing_score.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_pausing_args_parser()

    print('Step2: import the pausing file.', flush=True)
    pausing_df = import_pausing_files(args.list)

    print('Step3: output the pausing table.', flush=True)
    output_table(pausing_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
