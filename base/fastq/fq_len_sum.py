#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2020/4/1 21:06
# @Project : mytools
# @Script  : fq_len_sum.py
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
        description="This script is used to sum the reads length.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
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
def read_fq(in_fq_file, length_out_file):
    length_dict = {}
    sum_num = 0
    with open(in_fq_file, 'r') as fq_in:
        with open(length_out_file, 'w') as length_out:
            for line in fq_in:
                line2 = fq_in.__next__().strip()
                reads_length = len(line2)
                if length_dict.get(reads_length):
                    length_dict[reads_length] += 1
                else:
                    length_dict[reads_length] = 1
                fq_in.__next__()
                fq_in.__next__()

            length_out.writelines('length\treads_number\n')

            for length in sorted(length_dict.keys()):
                length_out.writelines(''.join([str(length), '\t', str(length_dict[length]), '\n']))
                sum_num += length_dict[length]
            length_out.writelines('Total reads number\t' + str(sum_num) + '\n')


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
    read_fq(args.i, args.o)
    clear()
    now_time()


# run
if __name__ == "__main__":
    main()
