#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_density.py


from foo import ArgsParser


def main():
    ArgsParser.now_time()
    print('\nConvert reads to RPFs density.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.ribo_args_parser()

    import foo.Ribo
    ribo_attr = foo.Ribo.Ribo(args)

    print('\nStep2: Import the P-site offset.\n', flush=True)
    ribo_attr.read_offset()

    print('\nStep3: Import the transcripts annotation.', flush=True)
    ribo_attr.read_transcript()

    print('\nStep4: Import the BAM file.\n', flush=True)
    ribo_attr.read_bam()

    print('\nStep5: Format the in-frame RPFs density.', flush=True)
    ribo_attr.run_format_rpf_with_multi_thread()

    print('\nStep5: Output the RPFs density.', flush=True)
    ribo_attr.output_density()

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
