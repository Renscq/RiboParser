#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Retrieve.py


from utils.ribo import ArgsParser
from utils.ribo.Retrieve import *


def main():
    ArgsParser.now_time()
    print('\nRetrieve the RPFs with gene list.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.retrieve_args_parser()
    rpfs = Retrieve(args)

    # print('Step2: Import gene list.\n', flush=True)
    # rpfs.import_gene_list()

    print('\nStep2: Retrieve the gene RPFs.', flush=True)
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    print('\nStep3: Format the RPFs table.', flush=True)
    rpfs.melt_rpf_table()

    print('\nStep4: Output the RPFs table.', flush=True)
    rpfs.output_rpf_table()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
