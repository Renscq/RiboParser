#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/11/30 21:46
# @Project : riboParser
# @Script  : rpf_occupancy.py

from utils.ribo import ArgsParser
from utils.ribo import Occupancy


def main():
    ArgsParser.now_time()
    print('\nCalculate the codon occupancy.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_occupancy_args_parser()
    occupancy = Occupancy.Occupancy(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    occupancy.import_rpf()

    print('\nStep3: Calculate the codon occupancy.', flush=True)
    occupancy.codon_occupancy()

    print('\nStep4: Draw the codon occupancy plot.', flush=True)
    occupancy.draw_occupancy_corr()
    occupancy.draw_occupancy_heat()
    occupancy.draw_occupancy_relative_heat()
    occupancy.draw_occupancy_line()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
