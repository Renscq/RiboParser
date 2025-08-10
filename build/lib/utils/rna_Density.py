#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rna_density.py


from utils.ribo import ArgsParser
from utils.ribo import RNA


def main():
    ArgsParser.now_time()
    print('\nConvert reads to reads density.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rna_args_parser()

    rna_attr = RNA.RNA(args)

    print('\nStep2: Import the P-site offset.', flush=True)
    rna_attr.read_offset()

    print('\nStep3: Import the transcripts annotation.', flush=True)
    rna_attr.read_transcript()

    print('\nStep4: Import the BAM file.', flush=True)
    rna_attr.read_bam()
    rna_attr.calculate_density()

    print('\nStep5: Format the in-frame reads density.', flush=True)
    rna_attr.run_format_reads_with_multi_thread()

    print('\nStep6: Output the reads density.', flush=True)
    rna_attr.output_density()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
