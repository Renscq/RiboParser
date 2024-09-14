#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboparser
# @Script  : rpf_Bam_Filter.py


from foo import ArgsParser


def main():
    ArgsParser.now_time()
    print('Filter the specific length reads from bam file.\n', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = ArgsParser.bam_filter_args_parser()

    from foo import BamFilter
    bam_attr = BamFilter.BamFilter(args)

    print('\nStep2: Filter the alignment bam file.', flush=True)
    bam_attr.import_bam()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

