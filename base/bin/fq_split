#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : fq_split.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/05/14 11:07:12
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import argparse
import gc
import time
import sys
import gzip
from argparse import RawTextHelpFormatter


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to split the fastq file with UMI-adapter file.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-i", metavar="input", required=True, type=str, help="input fastq file"
    )
    parser.add_argument(
        "-a", metavar="adapter", required=True, type=str, help="input adapter table file, ID for output fastq file, Adapter for retrieve."
    )
    parser.add_argument(
        "-o", metavar="output", required=True, type=str, help="prefix of output file"
    )
    # adapter like this
    
    args = parser.parse_args()
    return args


# output splitd sequences
def fqsplit(fq_input, fq_output):
    r1_name = fq_output + ".R1.fq"

    with open(r1_name, 'w') as r1_out:
        with open(r2_name, 'w') as r2_out:

                for line in fq_input:
                    r1_out.writelines(line)
                    r1_out.writelines(fq_input.__next__())
                    r1_out.writelines(fq_input.__next__())
                    r1_out.writelines(fq_input.__next__())

                    r2_out.writelines(fq_input.__next__())
                    r2_out.writelines(fq_input.__next__())
                    r2_out.writelines(fq_input.__next__())
                    r2_out.writelines(fq_input.__next__())
                    
# gc collect
def clear():
    print("Clear Memory.")
    gc.collect()
    print("All done.\n")


# get time
def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sys.stdout.flush()


# main
def main():
    now_time()
    args = parse_args()
    print(args)
    if args.i.split('.')[-1] == 'gz':
        with gzip.open(args.i, 'rt') as fq_input:
            fqsplit(fq_input, args.o)
    else:
        with open(args.i, 'r') as fq_input:
            fqsplit(fq_input, args.o)
    clear()
    now_time()


if __name__ == '__main__':
    main()
