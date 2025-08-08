#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_CST.py


from utils.ribo import ArgsParser
from utils.ribo import CST


def main():
    ArgsParser.now_time()
    print('\nCalculate the codon decoding time.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_cst_args_parser()

    cst = CST.CodonSelectiveTime(args)

    print('\nStep2: Import the Ribo and mRNA density file.', flush=True)
    cst.import_density()

    # print('\nStep3: Combine the Ribo and mRNA density file.', flush=True)
    # cst.merge_rpf_rna()

    print('\nStep3: Calculate the codon selection time.', flush=True)
    cst.calc_cst()

    print('\nStep4: Output the codon selection time.', flush=True)
    cst.format_cst_results()
    cst.output_cst()

    print('\nStep5: Draw the codon decoding time plot.', flush=True)
    # cst.draw_cst_line()
    cst.draw_cst_corr()
    cst.draw_cst_heat()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
