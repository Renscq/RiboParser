#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_digest.py


from utils.ribo.ArgsParser import *
from utils.ribo.Digestion import *


def main():
    print('\nDetect the digestion sites.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = digestion_args_parser()
    ribo_attr = Ribo(args)

    print('\nStep2: Import the annotation of transcripts.', flush=True)
    ribo_attr.read_transcript()

    print('\nStep3: Detect the digestion sites.', flush=True)
    ribo_attr.get_digest_sites()

    print('\nStep4: Output the digestion sites.', flush=True)
    ribo_attr.output_digest_sites()

    print('\nStep5: Draw the heatmap of digestion sites.', flush=True)
    ribo_attr.digestion_plot()

    print('\nStep6: Draw the seq logo of digestion sites.', flush=True)
    ribo_attr.output_counts()
    
    ribo_attr.seq_logo_plot2()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == '__main__':
    main()
