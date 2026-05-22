#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_metagene.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_metagene_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the metagene files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for metagene files (e.g., '*_metaplot.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _metaplot.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_metagene_files(metagene_list):
    '''
    @Message  : retrieve the metagene table.
    @Input    : list --> pattern for metagene files (e.g., '*tis_tts_metaplot.txt')
    @Return   : 
                output --> output dataframe contain the metagene
    @Flow     : step1 --> retrieve the metagene
    '''

    # make the dataframe
    metagene_colunms = ['Sample', 'Meta', 'Nucleotide', 'Codon', 'Frame', 'Density']
    metagene_df = pd.DataFrame(columns=metagene_colunms)

    # for each gne file
    for metagene_file in metagene_list:

        # import the metagene file
        metagene = pd.read_csv(metagene_file, sep='\t', index_col=False)

        # merge the gene metagene file
        metagene_df = pd.concat([metagene_df, metagene], axis=0)

    return metagene_df


def output_table(metagene_df, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the metagene
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    metagene_df.to_csv(output_file + '_metagene.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_metagene_args_parser()

    print('Step2: import the metagene file.', flush=True)
    metagene_df = import_metagene_files(args.list)

    print('Step3: output the metagene table.', flush=True)
    output_table(metagene_df, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
