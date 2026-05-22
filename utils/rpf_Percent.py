#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_coverage.py


from utils.ribo import ArgsParser
from utils.ribo import Percentage


def main():
    ArgsParser.now_time()
    print('\nDraw the metagene coverage.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.percentage_args_parser()
    percent = Percentage.Percentage(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    percent.read_rpf()
    percent.import_gene()

    print('\nStep3: Calculate the percentage of RPFs coverage.', flush=True)
    percent.calc_density_percent()

    print('\nStep4: Draw the histogram and boxplot of RPFs coverage.', flush=True)
    percent.draw_rpf_histogram()
    percent.draw_rpf_boxplot()

    print('\nStep5: Output the RPFs coverage.', flush=True)
    percent.output_density_percent()

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
