#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : rpf_Shuffle.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/04/09 17:44:44
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


from foo import ArgsParser
from foo.Shuffle import *


def main():
    ArgsParser.now_time()
    print('Shuffle the RPFs data.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.shuffle_args_parser()
    rpfs = Shuffle(args)

    print('Step2: Import the RPFs.\n', flush=True)
    rpfs.import_rpf()

    print('Step3: Shuffle the RPFs table.\n', flush=True)
    rpfs.shuffle_rpfs()

    print('Step4: Output the RPFs table.\n', flush=True)
    rpfs.output_rpfs()

    ArgsParser.now_time()
    print('All done.', flush=True)


if __name__ == '__main__':
    main()
