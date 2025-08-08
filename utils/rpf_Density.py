#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_density.py


from utils.ribo import ArgsParser
from utils.ribo import Ribo


def main():
    ArgsParser.now_time()
    print('\nConvert reads to RPFs density.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.ribo_args_parser()

    ribo_attr = Ribo.Ribo(args)

    print('\nStep2: Import the P-site offset.', flush=True)
    ribo_attr.read_offset()

    print('\nStep3: Import the transcripts annotation.', flush=True)
    ribo_attr.read_transcript()
    ribo_attr.check_transcript()

    print('\nStep4: Import the BAM file.', flush=True)
    ribo_attr.read_bam()

    print('\nStep5: Format the in-frame RPFs density.', flush=True)
    ribo_attr.run_format_rpf_with_multi_thread()

    print('\nStep5: Output the RPFs density.', flush=True)
    ribo_attr.output_density()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
