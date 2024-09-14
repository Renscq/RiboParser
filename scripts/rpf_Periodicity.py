#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_periodicity.py

from collections import OrderedDict

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from foo import ArgsParser
from foo import Periodicity

def main():
    ArgsParser.now_time()
    print('Draw the periodicity plot.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.periodicity_args_parser()

    print('Step2: Import the RPFs file.\n', flush=True)
    rpfs = Periodicity.Periodicity(args)
    rpfs.import_rpf()

    print('Step3: Calculate the 3nt periodicity.\n', flush=True)
    rpfs.calc_3nt_period()

    print('Step4: Ouput the 3nt periodicity.\n', flush=True)
    rpfs.output_meta()

    print('Step5: Draw the 3nt periodicity plot.\n', flush=True)
    rpfs.draw_3nt_period_count()
    rpfs.draw_3nt_period_ratio()

    print('All done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
