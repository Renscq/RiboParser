#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : merge_dst_list.py


import os
import pandas as pd
import argparse


def merge_density_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to create the list of density files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for density files (e.g., '*_rpf.txt')."
    )
    input_group.add_argument(
        '-o', dest='output', required=False, type=str, help='output file name (default: RPF.file.list).'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def process_density_files(density_list):
    '''
    @Message  : merge the list file path and name.
    @Input    : list --> pattern for density files (e.g., '*_rpf.txt')
    @Return   : 
                output --> output dataframe contain the reads offset
    @Flow     : step1 --> get the file list
                step2 --> create the file list
    '''

    # create the output list
    density_merge = []

    # for each density file
    for density_file in density_list:
        file_name = os.path.basename(density_file)
        file_path_name = os.path.abspath(density_file)

        file_prefix = file_name.split('.')[0]
        file_prefix = file_prefix.replace('_rpf', '')
        file_prefix = file_prefix.replace('_rna', '')

        if file_name.endswith('_rpf.txt'):
            density_model = 'Ribo'
        elif file_name.endswith('_rna.txt'):
            density_model = 'RNA'
        else:
            print('Error: the file name is not correct, please check it.')
            exit()

        # merge file path and name

        density_merge.append([file_prefix, file_path_name, density_model])

    return density_merge, density_model


def output_table(density_merge, density_model, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the reads offset
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''

    # set the result name
    if output_file:
        output_file = output_file
    else:
        if density_model == 'Ribo':
            output_file = 'RPF.file.list'
        elif density_model == 'RNA':
            output_file = 'RNA.file.list'

    # output the result
    density_merge_df = pd.DataFrame(density_merge, columns=['Name', 'File', 'Type'])

    density_merge_df.to_csv(output_file, sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = merge_density_args_parser()

    print('Step2: list the density file.', flush=True)
    density_file_list, density_model = process_density_files(args.list)

    print('Step3: output the density file list.', flush=True)
    output_table(density_file_list, density_model, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
