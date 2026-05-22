#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#@Project : riboParser
#@Script  : merge_offset_detail.py


import os
import pandas as pd
from collections import OrderedDict
import argparse


def merge_offset_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to merge the details of offset files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for offset end site files (e.g., '*end.txt')."
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

        offset_site = offset_model.replace('_3end', '').replace('_5end', '').upper()
        # if offset_model.replace('tis_', '').replace('tts_', '') == '3end': offset_end = '3end', else: offset_end = '5end'
        offset_end = "3p-end" if offset_model.replace('tis_', '').replace('tts_', '') == '3end' else "5p-end"

        # import the offset file
        file_prefix = file_name.replace('_3end', '').replace('_5end', '').replace('_tis', '').replace('_tts', '')
        offset_df = pd.read_csv(offset_file, sep='\t', index_col=0)
        
        try:
            offset_dict[file_prefix][offset_model] = offset_df

        except KeyError:
            offset_dict[file_prefix] = {}
            offset_dict[file_prefix][offset_model] = offset_df
            
        offset_dict[file_prefix][offset_model].index.name = 'Length'
        offset_dict[file_prefix][offset_model].reset_index(inplace=True)
        offset_dict[file_prefix][offset_model].insert(0, 'End', offset_end)
        offset_dict[file_prefix][offset_model].insert(0, 'Site', offset_site)
        offset_dict[file_prefix][offset_model].insert(0, 'Sample', file_prefix)

        # convert the table too the long table
        offset_dict[file_prefix][offset_model] = offset_dict[file_prefix][offset_model].melt(id_vars=['Sample', 'Length', 'Site', 'End'], var_name='Offset', value_name='Count')


    return offset_dict


def process_offset_dict(offset_dict):
    '''
    @Message  : retrieve the offset table.
    @Input    : list --> pattern for reads offset files (e.g., '*end.txt')
    @Return   : 
                output --> output dataframe contain the reads offset
    @Flow     : step1 --> retrieve the reads offset
    '''

    offset_5end_merge = pd.DataFrame()
    offset_3end_merge = pd.DataFrame()

    # for each gene file
    for file, offset in offset_dict.items():
        
        # for each offset model
        offset_5end = pd.concat([offset['tis_5end'], offset['tts_5end']], axis=0).reset_index(drop=True)
        offset_3end = pd.concat([offset['tis_3end'], offset['tts_3end']], axis=0).reset_index(drop=True)

        # merge the gene offset file
        if offset_5end_merge.empty:
            offset_5end_merge = offset_5end
        else:
            offset_5end_merge = pd.concat([offset_5end_merge, offset_5end], axis=0)

        if offset_3end_merge.empty:
            offset_3end_merge = offset_3end
        else:
            offset_3end_merge = pd.concat([offset_3end_merge, offset_3end], axis=0)

    offset_3end_merge = offset_3end_merge.dropna(axis=1, how='all')
    offset_5end_merge = offset_5end_merge.dropna(axis=1, how='all')
    offset_merged = pd.concat([offset_5end_merge, offset_3end_merge], axis=0)

    return offset_merged



def output_table(offset_merged, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> dataframe contain the reads offset
    @Return   : output --> description
    @Flow     : step1 --> output the dataframe to txt file
    '''
    
    offset_merged.to_csv(output_file + '_offset_end.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_offset_args_parser()

    print('Step2: import the end of offset file.', flush=True)
    offset_dict = import_offset_files(args.list)
    offset_merged = process_offset_dict(offset_dict)

    print('Step3: output the end of offset table.', flush=True)
    output_table(offset_merged, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
