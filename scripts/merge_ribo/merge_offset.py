#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_offset.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_offset_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the reads offset files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for RSBM/SSCBM offset files (e.g., '*RSBM_offset.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name (default: prefix + _offset.txt).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def process_offset_files(offset_list):
    '''
    @Message  : retrieve the offset table.
    @Input    : list --> pattern for reads offset files (e.g., '*RSBM_offset.txt')
    @Return   : 
                output --> output dataframe contain the reads offset
    @Flow     : step1 --> retrieve the reads offset
    '''

    # make the dataframe
    offset_colunms = ['Sample', 'Model', 'Length',
                      'Frame0', 'Rpfs0', 'Frame1', 'Rpfs1', 'Frame2', 'Rpfs2',
                      'P_site', 'Rpfs', 'Periodicity', 'Ribo']

    offset_merge = pd.DataFrame(columns=offset_colunms)

    # for each gne file
    for offset_file in offset_list:
        file_name = os.path.basename(offset_file).split('.')[0]

        if file_name.endswith('_RSBM_offset'):
            offset_model = 'RSBM'
        elif file_name.endswith('_SSCBM_offset'):
            offset_model = 'SSCBM'
        else:
            print('Error: the file name is not correct, please check it.')
            exit()

        # import the offset file
        file_prefix = file_name.replace('_RSBM_offset', '').replace('_SSCBM_offset', '')
        offset_df = pd.read_csv(offset_file, sep='\t', index_col=False)

        # convert the first chr of columns to upper
        offset_df.columns = [i[0].upper() + i[1:] for i in offset_df.columns]
        
        offset_df['Model'] = offset_model
        offset_df['Sample'] = file_prefix
        offset_df = offset_df[offset_colunms]

        # merge the gene offset file
        offset_merge = pd.concat([offset_merge, offset_df], axis=0)

    return offset_merge


def output_table(offset_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the reads offset
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    # keep the 2 decimal places
    int_columns = ['Length', 'Frame0', 'Rpfs0', 'Frame1', 'Rpfs1', 'Frame2', 'Rpfs2', 'P_site', 'Rpfs']
    offset_merge[int_columns] = offset_merge[int_columns].astype('int')

    offset_merge.to_csv(output_file + '_offset.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_offset_args_parser()

    print('Step2: import the reads offset file.', flush=True)
    offset_merge = process_offset_files(args.list)

    print('Step3: output the reads offset table.', flush=True)
    output_table(offset_merge, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
