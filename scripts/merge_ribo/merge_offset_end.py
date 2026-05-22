#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_offset_end.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_offset_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the details of tis offset files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for frame/tis offset files (e.g., '*5end.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='output file name.'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def import_offset_files(offset_list):
    '''
    @Message  : retrieve the offset table.
    @Input    : list --> pattern for reads offset files (e.g., '*end.txt')
    @Return   : 
                output --> output dataframe contain the reads offset
    @Flow     : step1 --> retrieve the reads offset
    '''

    offset_dict = OrderedDict()

    # for each gene file
    for offset_file in offset_list:
        file_name = os.path.basename(offset_file).split('.')[0]

        if file_name.endswith('tis_3end'):
            offset_model = 'tis_3end'
        elif file_name.endswith('tis_5end'):
            offset_model = 'tis_5end'
        elif file_name.endswith('tts_3end'):
            offset_model = 'tts_3end'
        elif file_name.endswith('tts_5end'):
            offset_model = 'tts_5end'
        else:
            print('Error: the file name is not correct, please check it.')
            exit()

        # import the offset file
        file_prefix = file_name.replace('_3end', '').replace('_5end', '').replace('_tis', '').replace('_tts', '')
        offset_df = pd.read_csv(offset_file, sep='\t', index_col=0)
        
        try:
            offset_dict[file_prefix][offset_model] = offset_df
        except KeyError:
            offset_dict[file_prefix] = {}
            offset_dict[file_prefix][offset_model] = offset_df

    return offset_dict


def process_offset_dict(offset_dict):
    '''
    @Message  : retrieve the offset table.
    @Input    : list --> pattern for reads offset files (e.g., '*end.txt')
    @Return   : 
                output --> output dataframe contain the reads offset
    @Flow     : step1 --> retrieve the reads offset
    '''

    offset_merge = pd.DataFrame()

    # for each gene file
    for file, offset in offset_dict.items():
        
        offset_df = pd.DataFrame()

        # for each offset model
        tis = pd.concat([offset['tis_5end'], offset['tis_3end']], axis=1)
        tts = pd.concat([offset['tts_5end'], offset['tts_3end']], axis=1)
        
        tis.index.name = 'Length'
        tts.index.name = 'Length'

        tis.reset_index(inplace=True)
        tts.reset_index(inplace=True)

        tis.insert(0, 'Site', 'TIS')
        tts.insert(0, 'Site', 'TTS')

        tts.insert(0, 'Sample', file)
        tis.insert(0, 'Sample', file)

        offset_df = pd.concat([tis, tts], axis=0)

        # merge the gene offset file
        if offset_merge.empty:
            offset_merge = offset_df
        else:
            offset_merge = pd.concat([offset_merge, offset_df], axis=0)

    return offset_merge


def output_table(offset_merge, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> dataframe contain the reads offset
    @Return   : output --> description
    @Flow     : step1 --> output the dataframe to txt file
    '''

    offset_merge.to_csv(output_file + '_offset_end.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_offset_args_parser()

    print('Step2: import the end of offset file.', flush=True)
    offset_dict = import_offset_files(args.list)
    offset_merge = process_offset_dict(offset_dict)

    print('Step3: output the end of offset table.', flush=True)
    output_table(offset_merge, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
