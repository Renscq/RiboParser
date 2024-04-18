#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_RSCU.py


from foo import RSCU
from foo import ArgsParser


def main():
    ArgsParser.now_time()
    print('\nEvaluate the different properties of sequence.', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = ArgsParser.codon_args_parser()

    print('\nStep2: Import the sequence.', flush=True)
    codon = RSCU.Codon(args)
    codon.read_rpf()

    print('\nStep3: Calculate the codon usage.', flush=True)
    codon.get_all_codon_usage()
    codon.get_codon_usage()

    print('\nStep4: Draw the codon usage plot.', flush=True)
    codon.draw_codon_usage()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

