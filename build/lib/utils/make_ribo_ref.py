#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : make_ribo_ref.py


from utils.ribo import ArgsParser
from utils.ribo.GenePred import *


def main():
    ArgsParser.now_time()
    print('\nMake the reference for riboParser.', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = ArgsParser.make_ribo_ref()
    ribo_ref = GenePred(args)

    print('\nStep2: Import genome sequence.', flush=True)
    ribo_ref.read_genome()

    print('\nStep3: Format gtf annotation.', flush=True)
    ribo_ref.gtf2gp()
    ribo_ref.read_genepred()
    ribo_ref.get_rep_transcript()
    ribo_ref.add_utr()

    print('\nStep4: Output gtf annotation.', flush=True)
    ribo_ref.write_txt()
    ribo_ref.gp2gtf()

    print('\nStep5: Retrieve mRNA sequence.', flush=True)
    ribo_ref.get_seq()
    ribo_ref.write_seq()

    print('\nAll done.', flush=True)


if __name__ == '__main__':
    main()
