#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_merge.py


from collections import OrderedDict
import pandas as pd
import polars as pl

from utils.ribo import ArgsParser


def read_sample_list(sample_list):
    sp_list = pd.read_csv(sample_list, header=0, names=None, sep = '\t')

    sp_dict = OrderedDict()
    
    for _, rec in sp_list.iterrows():
        sp_dict[rec.iloc[0]] = rec.iloc[1]

    return sp_dict


def merge_rpf(sp_dict, output_prefix):
    total_rpf = None
    now_sample_num = 0
    for sp_name, sp_file in sp_dict.items():
        print('Import the RPFs file: {file_name}.'.format(file_name=sp_name), flush=True)
        if now_sample_num == 0:
            rpf = pl.read_csv(sp_file, separator='\t', has_header=True)
            rpf.columns = ["name", "now_nt", "from_tis", "from_tts", "region", "codon",
                           sp_name + "_f0", sp_name + "_f1", sp_name + "_f2"]
            total_rpf = rpf
        else:
            rpf = pl.read_csv(sp_file, separator='\t', has_header=True)
            rpf = rpf[:, [0, 1, 6, 7, 8]]
            rpf.columns = ["name", "now_nt", sp_name + "_f0", sp_name + "_f1", sp_name + "_f2"]
            total_rpf = total_rpf.join(rpf, on=['name', 'now_nt'])
        now_sample_num += 1

    output_name = output_prefix + "_merged.txt"
    total_rpf = total_rpf.with_columns(total_rpf['codon'].str.to_uppercase())
    total_rpf = total_rpf.filter(~total_rpf['codon'].str.contains('N'))

    total_rpf.write_csv(output_name, separator = '\t', include_header = True)


def main():
    ArgsParser.now_time()
    print('\nMerge RPFs files from different samples.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_merge_args_parser()

    print('\nStep2: Import the sample list.', flush=True)
    sp_list = read_sample_list(args.list)

    print('\nStep3: Merge the RPFs file.', flush=True)
    merge_rpf(sp_list, args.output)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
