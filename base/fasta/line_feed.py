#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : line_feed.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/11/07 18:00:49
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


# from Bio import SeqIO
import argparse


def line_feed_args_parser():
    
    parser = argparse.ArgumentParser(description='This script is used to convert multiple lines to one line.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        '-i', dest='input', required=True, type=str, help='input fasta file name.'
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='output fasta file name.'
    )

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print('{:<12}:  {:<}'.format(k, str(v)), flush=True)

    return args


def line_feed(args):
    '''
    @Message  : function for convert multiple lines to one line.
    @Input    : fasta --> input file
    @Return   : output --> fasta file
    '''

    with open(args.input, 'r') as fa_in:
        with open(args.output, 'w') as fa_out:
            seq_id = ''
            seq_mess = ''

            for line in fa_in:
                if line.startswith('>') and seq_mess == '':
                    seq_id = line.strip()

                elif line.startswith('>') and seq_mess != '':
                    fa_out.writelines(seq_id + '\n')
                    fa_out.writelines(seq_mess + '\n')

                    seq_id = line.strip()
                    seq_mess = ''
                else:
                    seq_mess += line.strip()
                
            fa_out.writelines(seq_id + '\n')
            fa_out.writelines(seq_mess + '\n')


def main():
    args = line_feed_args_parser()
    line_feed(args)


if __name__ == '__main__':
    main()
