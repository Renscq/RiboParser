#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2020/8/2 23:17
# @Project : py-scripts
# @Script  : metaBg.py
# @Version : python 3.8.5
# @product : PyCharm
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
import gc
import time
from argparse import RawTextHelpFormatter


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to convert the bedgraph to plot data.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-i", metavar="input", required=True, type=str, help="input file"
    )
    parser.add_argument(
        "-o", metavar="output", required=True, type=str, help="output file"
    )
    # parser.add_argument(
    #     "-t", metavar="type", required=False, type=str, default='t',
    #     help="RPKM or TPM. (default: %(default)s)"
    # )

    args = parser.parse_args()
    return args


# convert the bedgraph to each nucleotide
def convert_bg(in_bg_file, out_bg_file):
    # pos_minus_list = []
    pos_plus_list = []
    read_minus_list = []
    read_plus_list = []
    with open(in_bg_file, 'r') as in_bg:
        for line in in_bg:
            rec = line.strip().split('\t')
            if rec[-1][0] == '-':
                pos_minus = [rec[0]+"\t"+str(i) for i in range(int(rec[1])+1, int(rec[2])+1)]
                read_minus = [rec[3] for i in range(int(rec[1]), int(rec[2]))]
                # pos_minus_list.extend(pos_minus)
                read_minus_list.extend(read_minus)
            else:
                pos_plus = [rec[0]+"\t"+str(i) for i in range(int(rec[1])+1, int(rec[2])+1)]
                read_plus = [rec[3] for i in range(int(rec[1]), int(rec[2]))]
                pos_plus_list.extend(pos_plus)
                read_plus_list.extend(read_plus)

    with open(out_bg_file, 'w') as out_bg:
        for line in range(len(pos_plus_list)):
            out_bg.writelines('\t'.join([pos_plus_list[line], read_plus_list[line]]) + '\n')
            # out_bg.writelines('\t'.join([pos_plus_list[line], read_plus_list[line], read_minus_list[line]]) + '\n')


# gc collect
def clear():
    print("Clear Memory.")
    gc.collect()
    print("All done.")


# get time
def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


# main
def main():
    now_time()
    args = parse_args()
    convert_bg(args.i, args.o)
    clear()
    now_time()


if __name__ == '__main__':
    main()
