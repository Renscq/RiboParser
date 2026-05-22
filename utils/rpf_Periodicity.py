#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_periodicity.py


from utils.ribo import ArgsParser
from utils.ribo import Periodicity

def main():
    ArgsParser.now_time()
    print('\nDraw the periodicity plot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.periodicity_args_parser()

    print('\nStep2: Import the RPFs file.', flush=True)
    rpfs = Periodicity.Periodicity(args)
    rpfs.import_rpf()

    print('\nStep3: Calculate the 3nt periodicity.', flush=True)
    rpfs.calc_3nt_period()

    print('\nStep4: Ouput the 3nt periodicity.', flush=True)
    rpfs.output_meta()

    print('\nStep5: Draw the 3nt periodicity plot.', flush=True)
    rpfs.draw_3nt_period_count()
    rpfs.draw_3nt_period_ratio()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
