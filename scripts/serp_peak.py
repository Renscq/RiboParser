#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : find_peak.py


import time

from foo import ArgsParser
from foo.SeRP import *


# main programme is here
def main():
    ArgsParser.now_time()
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.serp_peak_args_parser()
    serp = SeRP(args)
    serp.args_check()

    print('\nStep2: Import the RPF data.\n', flush=True)
    serp.rpf_txt_read()

    print('\nStep3: Import the gene annotation.\n', flush=True)
    serp.gene_anno()

    print('\nStep4: Scan peaks from the data.\n', flush=True)
    serp.detect_binding_peaks()

    print('\nStep5: Output peaks results.\n', flush=True)
    serp.output_peak()

    print("\nAll done!")
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
