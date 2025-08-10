#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_meta_codon.py


from utils.ribo.ArgsParser import *
from utils.ribo.MetaCodon import *


def main():
    now_time()

    print('\nDraw the meta-codon plot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = meta_codon_plot_args_parser()
    meta_codon = MetaCodon(args)

    print('\nStep2: Import the codon list.', flush=True)
    meta_codon.import_codon()

    print('\nStep3: Import the RPFs file.', flush=True)
    meta_codon.import_rpf()
    meta_codon.smooth_rpf_density()

    print('\nStep4: Retrieve specific codon density.', flush=True)
    meta_codon.reterieve_codon_density()

    print('\nStep5: Output the specific codon meta RPFs density.', flush=True)
    meta_codon.output_meta_codon_density()

    print('\nStep6: Output the specific codon sequence.', flush=True)
    meta_codon.output_meta_codon_seq()

    print('\nStep7: Draw the specific codon.', flush=True)
    if args.fig:
        meta_codon.draw_meta_codon()


    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
