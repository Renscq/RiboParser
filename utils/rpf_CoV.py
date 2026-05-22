#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_CoV.py


from utils.ribo import ArgsParser
from utils.ribo import Coefficient_of_Variation


def main():
    ArgsParser.now_time()
    print('\nCalculate CoV in the CDS region.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_cov_args_parser()
    cov = Coefficient_of_Variation.CoV(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    cov.read_rpf()

    print('\nStep3: Calculate CoV.', flush=True)
    cov.cal_CoV()

    print('\nStep4: Output CoV.', flush=True)
    cov.output_CoV()

    if args.group:
        print('\nStep5: Compare CoV of each group.', flush=True)
        cov.read_group()
        cov.test_group_CoV()
    else:
        print('\nStep5: Skip CoV comparation.', flush=True)
        pass
    
    print('\nStep6: Draw CoV.', flush=True)
    cov.draw_fit_plot()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
