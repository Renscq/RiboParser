#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/11/30
# @Project : riboParser
# @Script  : rpf_CDT.py


from utils.ribo import ArgsParser
from utils.ribo import CDT


def main():
    ArgsParser.now_time()
    print('\nCalculate the codon decoding time.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_cdt_args_parser()
    cdt = CDT.CodonDecodingTime(args)

    print('\nStep2: Import the gene annotation.', flush=True)
    cdt.import_gene()

    print('\nStep3: Import the RPFs file.', flush=True)
    cdt.import_density()
    cdt.merge_rpf_rna()
    cdt.normalize_density()

    print('\nStep4: Calculate the codon decoding time.', flush=True)
    cdt.calc_rpf_norm()
    cdt.codon_decoding_time()

    print('\nStep5: Draw the codon decoding time plot.', flush=True)
    cdt.draw_cdt_corr()
    cdt.draw_cdt_heat()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
