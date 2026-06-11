#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_quant.py


from utils.ribo import Quant

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def rpf_quant_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to quantify RPFs in the CDS region.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (default: prefix + _rpf_quant.txt)")

    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="The number of codons after TIS will be discarded. (default: %(default)s).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="The number of codons before TES will be discarded. (default: %(default)s).")
    parser.add_argument('--utr5', dest="utr5", action='store_true', required=False, default=False,
                        help="quantification of 5'-utr. (default: %(default)s).")
    parser.add_argument('--utr3', dest="utr3", action='store_true', required=False, default=False,
                        help="quantification of 3'-utr. (default: %(default)s).")

    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def main():

    now_time()
    print('\nQuantify the RPFs in the different region.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = rpf_quant_args_parser()

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
    now_time()


if __name__ == '__main__':
    main()
