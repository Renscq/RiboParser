#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : find_peak.py


import time

from utils.ribo import ArgsParser
from utils.serp.SeRP import *


# main programme is here
def main():
    ArgsParser.now_time()
    
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.serp_peak_args_parser()
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
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
