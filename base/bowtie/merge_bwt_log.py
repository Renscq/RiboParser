#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : merge_bwt_log.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/01/16 15:47:02
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import os
import pandas as pd
from collections import OrderedDict
import argparse

def stat_bwt_args_parser():

    parser = argparse.ArgumentParser(description='This script is used to statstic the mapped reads from log files')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        "-l", "--list", nargs='+', required=True, help="List for bowtie mapping log files (e.g., '*log')."
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='prefix of output file name.'
    )

    # optional arguments
    parser.add_argument(
        "-n", "--name", required=False, default="rRNA,tRNA,ncRNA,mRNA,Genome", 
        help="set the name of each database (default: %(default)s)."
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def process_log_files(log_file_list, database_list=None):
    '''
    @Message  : retrieve the mapped reads from log files.
    @Input    : list --> pattern for log files (e.g., '*log')
                name --> column name to merge (['rRNA', 'tRNA', 'ncRNA', 'mRNA', 'Genome'])
    @Return   : 
                output --> output nested dict contain the mapped reads
    @Flow     : step1 --> retrieve the mapped reads and unmapped reads
    '''
    
    if database_list is None:
        database_list = ['rRNA', 'tRNA', 'ncRNA', 'mRNA', 'Genome']
    else:
        database_list = database_list.split(',')

    result_dict = OrderedDict()

    # for each log file
    for log_file in log_file_list:
        file_prefix = os.path.basename(log_file).split('.')[0]
        result_dict[file_prefix] = OrderedDict()

        with open(log_file, 'r') as file_in:
            mapped_sum = 0
            flag = -1

            # for each line in log file
            for line in file_in:
                # retrieve the total reads
                if flag == -1 and line.startswith("# reads processed:"):
                    result_dict[file_prefix]['Total'] = int(line.split(' ')[-1])
                    flag += 1
                # retrieve the mapped reads
                elif flag >= 0 and line.startswith("Reported"):
                    result_dict[file_prefix][database_list[flag]] = int(line.split(' ')[-2])
                    mapped_sum += int(line.split(' ')[-2])
                    flag += 1
                # skip the useless lines
                elif flag >= 0 and line.startswith("# reads with") or line.startswith("# reads that") or line.startswith("# reads processed:"):
                    pass
                # retrieve the unmapped reads
                else:
                    print('Error: the log file is not correct.', flush=True)

            result_dict[file_prefix]['Others'] = result_dict[file_prefix]['Total'] - mapped_sum

    return result_dict, database_list


def output_table(result_dict, database_list, output_file):
    '''
    @Message  : function for output.
    @Input    : result_dict --> nested dict contain the mapped reads
    @Return   : output --> description
    @Flow     : step1 --> convert the dict to dataframe
    '''
    
    mapped_df = pd.DataFrame.from_dict(result_dict, orient='index').fillna(0)
    ratio_df = mapped_df.div(mapped_df['Total'], axis=0) * 100

    mapped_df = mapped_df.drop('Total', axis=1).reset_index().melt(id_vars='index', var_name='Database', value_name='Count')
    ratio_df = ratio_df.drop('Total', axis=1).reset_index().melt(id_vars='index', var_name='Database', value_name='Ratio')

    mapped_df = mapped_df.rename(columns={'index': 'Sample'}).set_index(['Sample', 'Database'])
    ratio_df = ratio_df.rename(columns={'index': 'Sample'}).set_index(['Sample', 'Database'])

    result_df = mapped_df.join(ratio_df, how='outer').reset_index()
    result_df['Ratio'] = result_df['Ratio'].astype('float').round(2)
    
    database_list = database_list + ['Others']
    result_df['Database'] = pd.Categorical(result_df['Database'], categories=database_list, ordered=True)
    result_df = result_df.sort_values(by=['Sample', 'Database'])

    result_df.to_csv(output_file + '_mapping.txt', sep='\t', index=False)


def main():

    print('Step1: Checking the input Arguments.', flush=True)
    args = stat_bwt_args_parser()

    print('Step2: improt the bowtie log file.', flush=True)
    result_dict, database_list = process_log_files(args.list, args.name)

    print('Step3: output the mapping table.', flush=True)
    output_table(result_dict, database_list, args.output)

    print('All done.', flush=True)


if __name__ == '__main__':
    main()
