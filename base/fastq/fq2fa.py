#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Script  : fq2fa.py
# @Time    : 2020/7/10 16:00
# @project : riboparser
# @Version : python 3.7
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
import gc
import time
from argparse import RawTextHelpFormatter


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to convert the FQ to FA format.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    required = parser.add_argument_group(title="Required parameters")
    required.add_argument(
        "-i", metavar="input", required=True, type=str, help="input FQ sequence file"
    )
    required.add_argument(
        "-n", action='store_true',
        help="assign new id number for reads"
    )
    required.add_argument(
        "-o", metavar="out",  required=True, type=str,
        help="output file in FA format",
    )

    args = parser.parse_args()
    return args


# read the fastq file
def fq2fa(fq, new_id, fa):

    flag = 1
    with open(fq, "r") as inputs:
        with open(fa, "w") as outputs:
            if new_id:
                for line in inputs:
                    outputs.writelines(''.join(['>t', str(flag), '\n']))
                    outputs.writelines(''.join(inputs.__next__()))
                    inputs.__next__()
                    inputs.__next__()
                    flag += 1
            else:
                for line in inputs:
                    outputs.writelines(''.join(['>', line]))
                    outputs.writelines(''.join(inputs.__next__()))
                    inputs.__next__()
                    inputs.__next__()


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
    fq2fa(args.i, args.n, args.o)
    clear()
    now_time()


# run
if __name__ == "__main__":
    main()
