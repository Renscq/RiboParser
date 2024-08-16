#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboparser
@Script       : rpf_Cumulative_CoV.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/06/02 16:41:55
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


from foo import ArgsParser
from foo.Cumulative_CoV import *


def main():
    ArgsParser.now_time()
    print('Retrieve the RPFs with gene list.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.cumulative_cov_args_parser()
    rpfs = CumulativeCoV(args)

    print('Step2: Retrieve the gene RPFs.\n', flush=True)
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    print('Step3: Format the RPFs table.\n', flush=True)
    rpfs.melt_rpf_table()

    print('Step4: Calculate the cumulative CoV.\n', flush=True)
    rpfs.calc_cov()

    print('Step5: Output the cumulative CoV meta table.\n', flush=True)
    rpfs.merge_cov_table()

    print('Step5: Output the cumulative CoV table.\n', flush=True)
    rpfs.output_rpf_table()

    print('All done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
