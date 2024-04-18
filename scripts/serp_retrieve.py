#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : serp_retrieve_seq.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/05/11 18:35:32
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


from collections import OrderedDict
import pandas as pd
import polars as pl
from foo import ArgsParser


def import_bed(peak_file):
    peak_df = pl.read_csv(peak_file)
    peak_df = peak_df.to_pandas()
    
    return gene


def import_rpf(rpf_file):
    rpf = pl.read_csv(rpf_file)
    rpf = rpf.to_pandas()



    return gene_list


def retrieve_seq(gene_list, bed_in, up, down, bed_out):

    bed_dict = OrderedDict()
    with open(bed_in, 'r') as bed:
        for line in bed:
            rec = line.strip().split('\t')
            if rec[0] in bed_dict.keys():
                bed_dict[rec[0]].append([rec[1], rec[2]])
            else:
                bed_dict[rec[0]] = [[rec[1], rec[2]]]

    bed_fa = OrderedDict()
    for gene, ranges in bed_dict.items():
        now_gene = gene_list.get_group(gene)
        for now_range in ranges:
            start = int(now_range[0]) - up
            end = int(now_range[1]) + down
            now_seq = now_gene.loc[(now_gene['from_tis'] >= start) & (now_gene['from_tis'] < end), 'aa_list'].sum()
            seq_name = ' '.join([gene, str(start), str(end)])
            bed_fa[seq_name] = now_seq

    if bed_out == "input":
        bed_out = bed_out[:-3] + 'fa'
    else:
        bed_out = bed_out

    with open(bed_out, 'w') as fa_out:
        for seq_name, seqs in bed_fa.items():
            fa_out.writelines(''.join(['>', seq_name, '\n', str(seqs)]) + '\n')


def main():
    ArgsParser.now_time()
    print('\nRetrieve the sequence of peak region.', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = ArgsParser.serp_retrieve()

    print('\nStep2: Import the rpf file.', flush=True)
    gene_list = import_rpf(args.rpf)

    print('\nStep2: Import the rpf file.', flush=True)
    gene_list = import_rpf(args.rpf)

    print('\nStep4: output the peak sequence.', flush=True)
    retrieve_seq(gene_list, args.bed, args.up, args.down, args.output)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
