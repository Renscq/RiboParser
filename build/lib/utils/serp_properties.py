#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : serp_properties.py


from utils.ribo import ArgsParser
from utils.serp import Properties


def main():
    ArgsParser.now_time()
    
    print('\nEvaluate the different properties of sequence.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.serp_properties()

    print('\nStep2: Import the sequence.', flush=True)
    seq = Properties.Sequence(args)

    print('\nStep3: Calculate the gene codon usage.', flush=True)
    seq.create_codon_table()
    seq.calc_gene_codon_usage()

    print('\nStep4: Calculate the whole codon usage.', flush=True)
    seq.calc_whole_codon_usage()

    print('\nStep4: Calculate the properties of sequence.', flush=True)
    seq.protein_analysis()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
