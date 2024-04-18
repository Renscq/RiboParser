#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_TE.py


from foo.ArgsParser import *
from foo.TranslationEfficiency import *


def main():
    print('\nCalculate translation efficiency.\n', flush=True)
    print('\nStep1: Checking the input Arguments.\n', flush=True)
    args = calc_te_args_parser()
    gene_te = TranslationEfficiency(args)

    print('\nStep2: Import the Ribo-seq results generate from DESeq2.\n', flush=True)
    gene_te.import_riboseq()

    print('\nStep3: Import the RNA-seq results generate from DESeq2.\n', flush=True)
    gene_te.import_rnaseq()

    print('\nStep4: Draw the correlation of Ribo-seq and RNA-seq fold change.\n', flush=True)
    gene_te.merge_logfc()
    gene_te.draw_co_diff_genes('total')
    gene_te.draw_co_diff_genes('sig')

    print('\nStep5: Calculate the translation efficiency.\n', flush=True)
    gene_te.import_design()
    gene_te.calc_te()
    gene_te.output_te()

    print('\nStep6: Draw the cumulative of translation efficiency.\n', flush=True)
    gene_te.draw_cdf()

    print('\nStep6: Draw the volcano of translation efficiency.\n', flush=True)
    gene_te.draw_volcano()

    print('\nAll done.\n', flush=True)
    now_time()


if __name__ == '__main__':
    main()
