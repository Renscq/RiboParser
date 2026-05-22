#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2020/12/1 20:15
# @Project : mytools
# @Script  : fa_split.py
# @Version : python 3.8.5
# @product : PyCharm
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
import gc
import time
import sys
from argparse import RawTextHelpFormatter


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(description="This script is used to split the fasta files.")
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.1"
    )
    parser.add_argument(
        "-i", "--input", required=True, type=str, help="input cds sequences in fasta format"
    )
    parser.add_argument(
        "-n", "--num", required=False, type=int, default=10,
        help="number of files (default %(default)s)"
    )
    parser.add_argument(
        "-o", "--output", required=False, type=str, default='split',
        help="prefix of output file name (default %(default)s)"
    )

    args = parser.parse_args()

    args_dict = vars(args)
    for k, v in args_dict.items():
        print(k + ": " + str(v))

    sys.stdout.flush()
    return args


fa_list = []


# read the cdna sequence file
def read_fa(input_fa):
    with open(input_fa, 'r') as fa_read:
        ID, seq = '', ''
        for line1 in fa_read:
            if line1.startswith('>') and ID != '':
                fa_list.append(ID)
                fa_list.append(seq)
                ID = line1.strip()
                seq = ''
            elif line1.startswith('>'):
                ID = line1.strip()
            else:
                seq = seq + line1.strip()
        fa_list.append(ID)
        fa_list.append(seq)


def split_fa(num, output_fa):
    lines = len(fa_list)
    flag = 0

    for file_num in range(num):
        with open(output_fa + str(file_num + 1) + ".fa", 'w') as fa_write:
            for line_num in range(0, int(lines / num), 2):
                if flag < lines:
                    fa_write.writelines('\n'.join([fa_list[flag], fa_list[flag + 1]]) + '\n')
                    flag += 2
                else:
                    break


# gc collect
def clear():
    gc.collect()
    print("All done.")
    sys.stdout.flush()


# get time
def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sys.stdout.flush()


# main
def main():
    now_time()
    args = parse_args()
    read_fa(args.input)
    if args.output:
        split_fa(args.num, args.output)
    else:
        output_fa = args.input.split('.')[0] + '-part'
        split_fa(args.num, output_fa)
    clear()
    now_time()


if __name__ == '__main__':
    main()
