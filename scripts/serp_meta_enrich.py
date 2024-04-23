#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : serp_enrich.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/04/20 16:52:56
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


from foo import ArgsParser
from foo import SeRP_Enrich


def main():
    ArgsParser.now_time()
    print('\nRetrieve the sequence of peak region.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.serp_retrieve_seq()
    enrich = SeRP_Enrich.SeRP_Enrich(args)

    print('\nStep2: Import the peak region.\n', flush=True)
    enrich.import_peak_region()

    print('\nStep3: Import the rpf file.\n', flush=True)
    enrich.import_rpf()

    print('\nStep4: Retrieve the sequence of peak region.\n', flush=True)
    enrich.retrieve_seq()

    print('\nStep5: Output the sequence of peak region.\n', flush=True)
    enrich.output_seq()

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
