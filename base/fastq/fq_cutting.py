#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2020/8/21 18:49
# @Project : py-scripts
# @Script  : fq_cutting.py
# @Version : python 3.8.5
# @product : PyCharm
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
import gc
import time
import sys
import gzip
from argparse import RawTextHelpFormatter


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to cut fastq file to many small files.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-i", metavar="input", required=True, type=str, help="input fastq file"
    )
    parser.add_argument(
        "-o", metavar="output", required=True, type=str, help="prefix of output file"
    )

    args = parser.parse_args()
    return args


# output splitd sequences
def fqsplit(fq_input, fq_output):
    r1_name = fq_output + ".R1.fq"
    r2_name = fq_output + ".R2.fq"

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
