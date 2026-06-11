#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : find_peak.py


import time

from utils.serp.SeRP import *

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def serp_peak_args_parser():
    parser = argparse.ArgumentParser(
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=5, width=90),
        description="This script is used to detect the selective translatome.")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    # arguments for the input files
    input_group = parser.add_argument_group('Input files arguments')
    input_group.add_argument("-r", dest="rpf", required=True, type=str, help="input the RPF coverage file")
    input_group.add_argument("-n", dest="norm", required=False, type=str, default=None,
                             help="specify the file contain RPFs number for normalise. "
                                  "(default: The total RPFs count of input samples is used.)")
    input_group.add_argument("--scale", dest="scale", required=False, type=float, default=1e6,
                             help="specify the scale level for reads normalise (default RPM %(default)s )")
    input_group.add_argument("-a", dest="anno", required=False, type=str,
                             help="specify the gene annotation file in txt format.")
    input_group.add_argument("--ck", dest="control", required=True, type=str,
                             help="specify the name of control samples, separated by commas.")
    input_group.add_argument("--ip", dest="ip", required=True, type=str,
                             help="specify the name of immunoprecipitation samples, separated by commas.")

    # arguments for the data filtering
    data_group = parser.add_argument_group('Data filtering arguments')
    data_group.add_argument("-m", dest="min", required=False, type=int, default=50,
                            help='''specify the minimum number of RPFs coverage for genes to be retained
             (default %(default)s)''')
    data_group.add_argument("--corr", dest="corr", required=False, type=float, default=0.3,
                            help="specify the minimum correlation of replicates (default %(default)s)")

    # arguments for the peaks scanning, gaps and ratio filtering
    peak_group = parser.add_argument_group('Peak scanning arguments')
    peak_group.add_argument("-f", dest="fill", required=False, type=int, default=30,
                            help='''Use the mean value of gene RPF to fill the missing values.
        (Commonly used first 30 AA, default %(default)s codon). Because the
         ribo-seq data is too noisy in single codon level and many locations 
         are not covered by RPF, it is necessary to fill the blank position 
         use the mean value of RPF as the background. Three modes are provided: 
         (1): [30] average of the first 30 AA  is used (or other length); 
         (2): [-1] the average RPFs value of current gene is used; 
         (3): [0] the average RPF of total genes is used.''')
    
    peak_group.add_argument("-s", dest="size", required=False, type=int, default=3,
                            help="Specifies the window size, the numbers must be odd. "
                                 "(default %(default)s AA)."
                                 "The sequencing data is usually noisy, and need to be smoothed."
                                 "savgol_filter is used here. ")
    peak_group.add_argument("-k", dest="k", required=False, type=int, default=1,
                            help="Specifies the polyorder, the numbers must be odd. "
                                 "(default %(default)s AA)."
                                 "The sequencing data is usually noisy, and need to be smoothed."
                                 "savgol_filter is used here. ")
    # peak_group.add_argument("--mode", dest="k", required=False, type=str, default="nearest", 
    #                         choices=["mirror", "nearest", "constant", "wrap"],
    #                         help="Specifies the smooth mode options, (default %(default)s AA).")
    
    peak_group.add_argument("-w", dest="width", required=False, type=int, default=5,
                            help="specify the width of binding peaks (default %(default)s AA)."
                                 "By default, the width of each peak is not less than 5 amino acids.")
    peak_group.add_argument("-e", dest="enrich", required=False, type=float, default=2.0,
                            help="specify the enrich threshold of each peak height (default %(default)s)")
    peak_group.add_argument("-c", dest="collision", required=False, type=float, default=1.5,
                            help="specify the enrich threshold of collision region of binding peak (default %(default)s)")
    peak_group.add_argument("-g", dest="gaps", required=False, type=int, default=1,
                            help='''specify the gap width inside each peak. (default %(default)s AA).
        The ribo-seq profiling is noisy at the single gene level (for example,
        the ribosomal density along a gene may change from a positive number
        to zero and, again, to a positive number). Therefore, gaps in the peak
        are allowed. By default, the gaps there can be no more than 2 consecutive
        amino acids in a peak.''')
    peak_group.add_argument("-p", dest="proportion", required=False, type=float, default=0.2,
                            help='''specify the proportion of gap in each peak (default %(default)s).
        Data with too much noise has low credibility, so it is necessary to
        add a limit to the gap proportion in each peak. By default, the total gaps
        length within a peak cannot exceed 20%%.''')
    peak_group.add_argument("--back", dest="background", required=False, type=int, default=0,
                            help='''use the 5 prime AA to filter the enrichment fold (default %(default)s codon,
        not applicable). Two modes are provided: 
        (1): [30](or other number > 0) Some proteins that bind to the nascent polypeptide 
        chain cannot affect the translation process of the first 30 AA. Theoretically, 
        the RPFs of the unbound area should be smaller than the RPFs of the bound area. 
        (2): [0] However, some proteins only affect the elongation of the ribosome, so 
        the peak may exist in any location.''')
    peak_group.add_argument("--bf", dest="backFold", required=False, action='store_true', default=True,
                            help='''specify the fold of background (default %(default)s).
        Enrichment after background region need greater than before (HSP70 binding model). 
        backFold used to filter the strong binding peak.
        ''')
    
    peak_group.add_argument("--up", dest="upstream", required=False, type=int, default=10,
                            help='''retrieve the upstream sequence of binding peaks, (default %(default)s codon.)''')
    peak_group.add_argument("--down", dest="downstream", required=False, type=int, default=10,
                            help='''retrieve the downstream sequence of binding peaks, (default %(default)s codon.)''')
    
    # arguments for output figures, ratio, peak annotation
    output_group = parser.add_argument_group('Output files arguments')
    output_group.add_argument("-o", dest="output", required=False, type=str, default='results',
                              help='''prefix of output file name (default %(default)s_ratio.txt).
        This prefix will be used as the output of Peak / Ratio / Matlab Script / Bed files.''')
    output_group.add_argument("--all", dest="all", required=False, action='store_true', default=False,
                              help='''output all the peak region. If this option on, all peak 
        regions include overlapped items in same gene are retained. Instead, 
        only the optimal peak for non-overlapping regions is output. (default %(default)s)''')
    output_group.add_argument("--rpm", dest="rpm", required=False, action='store_true', default=False,
                              help="output the rpm of each peak. (default %(default)s)")
    output_group.add_argument("--ratio", dest="ratio", required=False, action='store_true', default=False,
                              help="output the original ratio. (default %(default)s)")
    # output_group.add_argument(
    #     "--script", dest="script", required=False, action='store_true', default=True,
    #     help="output the MATLAB scripts for each figure. (default %(default)s)"
    # )
    output_group.add_argument("--fig", dest="fig", required=False, action='store_true', default=False,
                              help="draw the demo graph of peak scan results. (default %(default)s)."
                                   "This step may takes a lot of time.")

    args = parser.parse_args()
    file_check(args.rpf)
    file_check(args.anno)
    args_print(args)

    return args


# main programme is here
def main():
    now_time()
    
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = serp_peak_args_parser()
    serp = SeRP(args)
    serp.args_check()

    print('\nStep2: Import the RPF data.', flush=True)
    serp.rpf_txt_read()

    print('\nStep3: Import the gene annotation.', flush=True)
    serp.gene_anno()

    print('\nStep4: Scan peaks from the data.', flush=True)
    serp.detect_binding_peaks()

    print('\nStep5: Output peaks results.', flush=True)
    serp.output_peak()

    print("\nAll done!")
    now_time()


if __name__ == '__main__':
    main()
