#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Script  : fq2txt.py
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
        description="This script is used to convert the FQ to TXT format.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    required = parser.add_argument_group(title="Required parameters")
    required.add_argument(
        "-i", metavar="fq", required=True, type=str, help="input FQ sequence file"
    )
    required.add_argument(
        "-o", metavar="txt",  required=True, type=str,
        help="output file in TXT format",
    )

    args = parser.parse_args()
    return args


# read the database file of SNPs in VCF format
def fq2txt(fq, txt):

    with open(fq, "r") as inputs:
        with open(txt, "w") as outputs:
            for line in inputs:
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
    fq2txt(args.i, args.o)
    clear()
    now_time()


# run
if __name__ == "__main__":
    main()
