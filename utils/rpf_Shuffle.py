#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Script  : rpf_Shuffle.py


from utils.ribo import ArgsParser
from utils.ribo.Shuffle import *


def main():
    ArgsParser.now_time()
    print('\nShuffle the RPFs data.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.shuffle_args_parser()
    rpfs = Shuffle(args)

    print('\nStep2: Import the RPFs.', flush=True)
    rpfs.import_rpf()

    print('\nStep3: Shuffle the RPFs table.', flush=True)
    rpfs.shuffle_rpfs()

    print('\nStep4: Output the RPFs table.', flush=True)
    rpfs.output_rpfs()

    ArgsParser.now_time()
    print('\nAll done.', flush=True)


if __name__ == '__main__':
    main()
