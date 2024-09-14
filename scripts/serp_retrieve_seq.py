#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : serp_retrieve_seq.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/05/11 18:35:32
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


from foo import ArgsParser
from foo import SeRP_Seq

def main():
    ArgsParser.now_time()
    print('\nRetrieve the sequence of peak region.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.serp_retrieve_seq()
    peaks_region = SeRP_Seq.SeRPSeq(args)

    print('\nStep2: Import the peak region.\n', flush=True)
    peaks_region.import_peak_region()

    print('\nStep3: Import the rpf file.\n', flush=True)
    peaks_region.import_rpf()

    print('\nStep4: Retrieve the sequence of peak region.\n', flush=True)
    peaks_region.retrieve_seq()

    print('\nStep5: Output the sequence of peak region.\n', flush=True)
    peaks_region.output_seq()

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
