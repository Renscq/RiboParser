#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_odd_ratio.py


from utils.ribo import ArgsParser
from utils.ribo import Odd_Ratio


def main():
    ArgsParser.now_time()
    print('\nCalculate the relative codon pausing score.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_odd_ratio_args_parser()
    odd_ratio = Odd_Ratio.OddRatio(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    odd_ratio.read_rpf()

    print('\nStep3: Make two dimensional table for fisher exactly test.', flush=True)
    odd_ratio.make_two_dimensional_table()

    print('\nStep4: Calculate the odd ratio.', flush=True)
    odd_ratio.calc_odd_ratio()
    odd_ratio.calc_chi2_test2()

    print('\nStep5: Output the ribosome pausing odd ratio.', flush=True)
    odd_ratio.output_odd_ratio()
    odd_ratio.summarize_odd_ratio()
    
    print('\nStep6: Draw the ribosome pausing odd ratio.', flush=True)
    odd_ratio.draw_odd_ratio_line()
    odd_ratio.draw_odd_ratio_scatter()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
