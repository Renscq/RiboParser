#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : site2base.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/03/19 09:17:21
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


# import pandas as pd
# import polars as pl
# import numpy as np
# from collections import OrderedDict
# from Bio import SeqIO
import argparse



def rpm_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to retrieve the nucleotides from genome file with bedgraph.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-b', dest="input", required=True, type=str,
                             help="input bedgraph file name.")
    input_group.add_argument('-g', dest="input", required=True, type=str,
                             help="input genome file name.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="output nucleotides file name.")
    
    # arguments for the RPM calculation
    parser.add_argument('-u', dest="upstream", required=False, type=int, default=30,
                        help="upstream offset. (default: %(default)s ).")
    parser.add_argument('-d', dest="downstream", required=False,  type=int, default=0,
                        help="downstream offset. (default: %(default)s ).")
    parser.add_argument('-m', dest="min", required=False,  type=int, default=1,
                        help="min density of bedgraph (default: %(default)s ).")
    parser.add_argument('-r', dest="reshape", required=False,  action="store_true", default=False,
                        help="reshape the sequence to txt (default: %(default)s ).")
    args = parser.parse_args()

    return args
