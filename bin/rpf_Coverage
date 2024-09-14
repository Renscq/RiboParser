#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_coverage.py


from foo import ArgsParser
from foo import Coverage


def main():
    ArgsParser.now_time()
    print('\nDraw the metagene coverage.\n', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = ArgsParser.coverage_args_parser()
    meta = Coverage.Coverage(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    meta.read_rpf()
    meta.import_gene()

    print('\nStep3: Adjust mRNAs to the same dimension.', flush=True)
    meta.process_utr5()
    meta.process_cds()
    meta.process_utr3()

    print('\nStep5: Draw the line plot of metagene coverage.', flush=True)
    meta.draw_meta_gene_line()

    print('\nStep6: Draw the heatmap of metagene coverage.', flush=True)
    if args.heatmap:
        meta.draw_meta_gene_heat()

    print('\nStep7: Output the metagene coverage.', flush=True)
    meta.output_meta_gene()

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
