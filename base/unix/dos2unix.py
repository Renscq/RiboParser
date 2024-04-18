#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : dos2unix.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/11/29 16:02:13
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import sys
import argparse


def args_parser():
    parser = argparse.ArgumentParser(description='This script is used to convert the dos format to unix format.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-i', dest='input', required=True, type=str, help='input file name')
    input_group.add_argument('-o', dest='output', required=True, type=str, help='output file name')

    args = parser.parse_args()
    args_dict = vars(args)
    for k, v in args_dict.items():
        print("{:<12}:  {:<}".format(k, str(v)), flush=True)
    
    return args


def dos2unix(input_file, output_file):

    with open(output_file, 'w') as outfile:
        with open(input_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                print(line)
                outfile.writelines(''.join([line, '\n']))


def main():
    print("Convert the dos file to unix format.", flush=True)
    args = args_parser()

    print("\nStep1: output the unix file.", flush=True)
    dos2unix(args.input, args.output)


if __name__ == '__main__':
    main()
