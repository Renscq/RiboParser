#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboparser
# @Script  : rpf_Bam_Filter.py


from utils.ribo import ArgsParser


def main():
    ArgsParser.now_time()
    print('\nFilter the specific length reads from bam file.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.bam_filter_args_parser()

    from ribo import BamFilter
    bam_attr = BamFilter.BamFilter(args)

    print('\nStep2: Filter the alignment bam file.', flush=True)
    bam_attr.import_bam()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
