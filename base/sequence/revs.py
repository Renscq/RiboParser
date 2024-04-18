#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : revs.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/11/07 20:25:23
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


# import pandas as pd
# import polars as pl
# import numpy as np
# from collections import OrderedDict
# from Bio import SeqIO
import argparse

# get arguments and dim the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(description="converted sequence into the reverse complementary strand")
    parser.add_argument('seq',type=str)
    arguments = parser.parse_args()
    return arguments

def revSeq(seq):
    revsd = ""
    for nt in seq.strip():
        if nt == "A":
            nt = "T"
            revsd =  nt + revsd
        elif nt == "a":
            nt = "t"
            revsd =  nt + revsd
        elif nt == "T":
            nt = "A"
            revsd =  nt + revsd
        elif nt == "t":
            nt = "a"
            revsd =  nt + revsd
        elif nt == "G":
            nt = "C"
            revsd =  nt + revsd
        elif nt == "g":
            nt = "c"
            revsd =  nt + revsd
        elif nt == "C":
            nt = "G"
            revsd =  nt + revsd
        elif nt == "c":
            nt = "g"
            revsd =  nt + revsd
        elif nt == "U":
            nt = "A"
            revsd =  nt + revsd
        elif nt == "u":
            nt = "a"
            revsd =  nt + revsd
        elif nt == "N":
            revsd =  nt + revsd
        else:
            print('illegal nulcleotide in your sequence: %s'%(nt))
            revsd = nt + revsd
    print('reversed  :\t%s'%seq[::-1])
    print('complement:\t%s'%revsd[::-1])
    print('rev + comp:\t%s'%(revsd))

def main():
    args = parse_args()
    revSeq(args.seq)

if __name__ == '__main__':
    main()
