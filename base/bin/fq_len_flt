#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2020/4/1 21:06
# @Project : mytools
# @Script  : fq_len_flt.py
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
        description="This script is used to filter the reads length from fastq file.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-m", metavar="min", required=True, type=int, help="the min length",
    )
    parser.add_argument(
        "-M", metavar="max", required=True, type=int, help="the max length",
    )
    parser.add_argument(
        "-i", metavar="input", required=True, type=str, help="input file in fq format"
    )
    parser.add_argument(
        "-o", metavar="output", required=True, type=str, help="output file in fq format"
    )

    args = parser.parse_args()
    return args


# read the fastq file
def read_fq(in_fq_file, min_len, max_len, out_fq_file):
    dropped = out_fq_file + '.drop'
    with open(in_fq_file, 'r') as fq_in:
        with open(out_fq_file, 'w') as fq_out:
            with open(dropped, 'w') as fq_drop:
                for line in fq_in:
                    line1 = line.strip()
                    line2 = fq_in.__next__().strip()
                    line3 = fq_in.__next__().strip()
                    line4 = fq_in.__next__().strip()
                    if min_len <= len(line2) <= max_len:
                        fq_out.writelines('\n'.join([line1, line2, line3, line4]) + '\n')
                    else:
                        fq_drop.writelines('\n'.join([line1, line2, line3, line4]) + '\n')


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
    print(args)
    read_fq(args.i, args.m, args.M, args.o)
    clear()
    now_time()


# run
if __name__ == "__main__":
    main()
