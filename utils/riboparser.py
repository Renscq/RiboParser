#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : ribo_parser.py


from utils.ribo import ArgsParser
from utils.data import RiboParser


def main():
    ArgsParser.now_time()
    args = ArgsParser.ribo_parser()

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
    ArgsParser.now_time()


if __name__ == "__main__":
    main()
    
    
