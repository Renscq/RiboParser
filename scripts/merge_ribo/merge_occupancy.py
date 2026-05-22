#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_occupancy.py


import pandas as pd
import argparse


def merge_occupancy_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the occupancy files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for occupancy files (e.g., '*_codon_occupancy.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _codon_occupancy.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_occupancy_files(occupancy_list):
    '''
    @Message  : retrieve the occupancy table.
    @Input    : list --> pattern for occupancy files (e.g., '*_codon_occupancy.txt')
    @Return   : 
                output --> output dataframe contain the occupancy
    @Flow     : step1 --> retrieve the occupancy
    '''

    # make the dataframe
    occupancy_df = pd.DataFrame()

    # for each occupancy file
    for occupancy_file in occupancy_list:
        # import the occupancy file
        occupancy = pd.read_csv(occupancy_file, sep='\t', index_col=False)
        occupancy.sort_values(by=['Abbr', 'Codon'], ascending=True, inplace=True, ignore_index=True)

        # merge the gene occupancy file
        if occupancy_df.empty:
            occupancy_df = occupancy
        else:
            # delete the columns ['AA', 'Abbr.']
            occupancy = occupancy.drop(['AA', 'Abbr'], axis=1)
            occupancy_df = pd.merge(left=occupancy_df, right=occupancy, how='outer', on='Codon')

    # reorder the occupancy columns with 
    # ['codon', 'AA', 'Abbr.', *_density, *_absolute_occupancy, *_relative_occupancy]
    occupancy_codon = occupancy_df[['Codon', 'AA', 'Abbr']]
    occupancy_density = occupancy_df.filter(regex='_density')
    occupancy_abs = occupancy_df.filter(regex='_absolute_occupancy')
    occupancy_rel = occupancy_df.filter(regex='_relative_occupancy')

    occupancy_merge = pd.concat([occupancy_codon, occupancy_density, occupancy_abs, occupancy_rel], axis=1)

    return occupancy_merge


def output_table(occupancy_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the occupancy
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    occupancy_merge.to_csv(output_file + '_codon_occupancy.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_occupancy_args_parser()

    print('Step2: import the occupancy file.', flush=True)
    occupancy_df = import_occupancy_files(args.list)

    print('Step3: output the occupancy table.', flush=True)
    output_table(occupancy_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
