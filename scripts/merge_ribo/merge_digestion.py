#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_digestion.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_digestion_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the reads digestion files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for digestion files (e.g., '*_5end_pwm.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name.'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def process_digest_files(digestion_list):
    '''
    @Message  : retrieve the pwm table.
    @Input    : list --> pattern for reads digestion files (e.g., '*_pwm.txt')
    @Return   : 
                output --> output dataframe contain the reads digestion
    @Flow     : step1 --> retrieve the reads digestion
    '''

    # make the dataframe
    digest_colunms = ['Sample', 'Digest', 'Site', 'A', 'T', 'C', 'G']
    digestion_df = pd.DataFrame(columns=digest_colunms)

    # for each gne file
    for motif_file in digestion_list:
        file_name = os.path.basename(motif_file).split('.')[0]

        if file_name.endswith('_5end_pwm'):
            Digest = 'End_5p'
        elif file_name.endswith('_3end_pwm'):
            Digest = 'End_3p'
        else:
            print('Error: the file name is not correct, please check it.')
            exit()

        file_prefix = file_name.replace('_5end_pwm', '').replace('_3end_pwm', '')

        # import the motif file
        motif_df = pd.read_csv(motif_file, sep='\t', index_col=0)

        # rename index to Site
        motif_df.index.name = 'Site'
        motif_df.reset_index(inplace=True)

        # set the site range of 5' end and 3' end
        if Digest == 'End_5p':
            motif_df['Site'] = motif_df['Site'] - 5
        elif Digest == 'End_3p':
            motif_df['Site'] = motif_df['Site'] - 10

        motif_df['Digest'] = Digest
        motif_df['Sample'] = file_prefix
        motif_df = motif_df[digest_colunms]

        # merge the gene digestion file
        digestion_df = pd.concat([digestion_df, motif_df], axis=0)

    return digestion_df


def output_table(digestion_df, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the reads digestion motif
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    # keep the 2 decimal places
    digestion_df['Site'] = digestion_df['Site'].astype('int')

    digestion_df.to_csv(output_file + '_reads_digestion.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_digestion_args_parser()

    print('Step2: improt the reads digestion file.', flush=True)
    digestion_df = process_digest_files(args.list)

    print('Step3: output the reads digestion table.', flush=True)
    output_table(digestion_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
