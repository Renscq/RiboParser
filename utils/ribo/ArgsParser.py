#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : ArgsParser.py


import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')


def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)


def args_print(args):
    args_dict = vars(args)

    for k, v in args_dict.items():
        print("{:<12}:  {:<}".format(k, str(v)), flush=True)

    sys.stdout.flush()


def file_check(*files):
    for my_file in files:
        if os.path.exists(my_file):
            continue
        else:
            print("\nFile {file_name} is not exists!\n".format(file_name=my_file), flush=True)
            sys.exit()
