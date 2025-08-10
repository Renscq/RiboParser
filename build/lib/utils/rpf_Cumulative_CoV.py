#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Cumulative_CoV.py


from utils.ribo import ArgsParser
from utils.ribo.Cumulative_CoV import *


def main():
    ArgsParser.now_time()
    print('\nRetrieve the RPFs with gene list.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.cumulative_cov_args_parser()
    rpfs = CumulativeCoV(args)

    print('\nStep2: Retrieve the gene RPFs.', flush=True)
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    print('\nStep3: Format the RPFs table.', flush=True)
    rpfs.melt_rpf_table()

    print('\nStep4: Calculate the cumulative CoV.', flush=True)
    rpfs.calc_cov()

    print('\nStep5: Output the cumulative CoV meta table.', flush=True)
    rpfs.merge_cov_table()

    print('\nStep5: Output the cumulative CoV table.', flush=True)
    rpfs.output_rpf_table()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
