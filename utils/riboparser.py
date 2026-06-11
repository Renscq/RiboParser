#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : ribo_parser.py


from utils.data import RiboParser

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def ribo_parser():
    parser = argparse.ArgumentParser(description="Check the information of RiboParser.")

    # arguments for the modification
    parser.add_argument('-v', dest="version", required=False, action="store_true", default=False,
                        help="show the current version of RiboParser. (default: %(default)s).")
    parser.add_argument('-c', dest="citation", required=False, action="store_true", default=False,
                        help="show the citation of RiboParser. (default: %(default)s).")
    parser.add_argument('-d', dest="dependency", required=False, action="store_true", default=False,
                        help="check the dependency of RiboParser. (default: %(default)s).")
    parser.add_argument('-m', dest="module", required=False, action="store_true", default=False,
                        help="check the modules of RiboParser. (default: %(default)s).")

    args = parser.parse_args()
    # args_print(args)

    return args


def main():
    now_time()
    args = ribo_parser()

    if args.version:
        print('\nShow the version of RiboParser.', flush=True)
        RiboParser.RiboParserInfo.show_version()

    if args.citation:
        print('\nShow the citation of RiboParser.', flush=True)
        RiboParser.RiboParserInfo.show_citation()
    
    if args.dependency:
        print('\nShow the dependency of RiboParser.', flush=True)
        RiboParser.RiboParserInfo.check_dependencies()
    
    if args.module:
        print('\nShow the modules of RiboParser.', flush=True)
        RiboParser.RiboParserInfo.check_package_modules()

    print('')
    now_time()


if __name__ == "__main__":
    main()
    
    
