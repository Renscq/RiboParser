#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_quant.py


from utils.ribo import ArgsParser
from utils.ribo import Quant


def main():

    ArgsParser.now_time()
    print('\nQuantify the RPFs in the different region.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_quant_args_parser()

    quant = Quant.Quant(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    quant.read_rpf()

    print('\nStep3: Quantify the RPFs in different region.', flush=True)
    quant.quant_region()
    quant.output_total_rpf()

    print('\nStep4: Draw the RPFs bar plot of different region.', flush=True)
    quant.draw_rpf_barplot()

    print('\nStep5: Draw the RPFs cumulative plot of gene rpm.', flush=True)
    quant.draw_rpf_cdfplot()

    print('\nStep6: Draw the RPFs PCA plot of gene rpm.', flush=True)
    quant.draw_rpf_pcaplot()

    print('\nStep7: Draw the RPFs heatmap of gene rpm.', flush=True)
    quant.draw_rpf_heatmap2()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
