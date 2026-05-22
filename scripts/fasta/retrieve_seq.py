#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/8/28 23:00
# @Project : riboParser
# @Script  : retrieve_seq.py
# @Version : python 3.8.5
# @product : PyCharm
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
import time
import sys
from Bio import SeqIO


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(description="This script is used to retrieve the fasta sequence by name.")
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-i", dest="input", required=True, type=str,
        help="input the fasta file"
    )
    parser.add_argument(
        "-n", dest="name", required=True, type=str,
        help="gene ids in txt format"
    )
    parser.add_argument(
        "-u", dest="unmapped", required=False, type=str,default='unmapped.ids',
        help="output the unmapped gene ids"
    )
    parser.add_argument(
        "-o", dest="output", required=True, type=str, default='results',
        help="prefix of output file name (default %(default)s_peaks.txt)"
    )
    args = parser.parse_args()

    args_dict = vars(args)
    for k, v in args_dict.items():
        print("%12s : %s" % (k, str(v)))

    sys.stdout.flush()
    return args


def import_gene_name(gene_name_file):
    gene_name_list = []
    with open(gene_name_file, 'r') as gene_name_in:
        for line in gene_name_in:
            gene_name_list.append(line.strip('>').strip())

    return gene_name_list


def retrieve_seq(gene_name_list, fasta_file, output_name, unmapped):
    out_fa = open(output_name, 'wb')
    out_ids = open(unmapped, 'w')
    fa_seq = SeqIO.index(fasta_file, "fasta")
    for gene in gene_name_list:
        try:
            out_fa.write(fa_seq.get_raw(gene))
        except KeyError:
            out_ids.write(str(gene) + '\n')


# main
def main():

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sys.stdout.flush()
    args = parse_args()

    print('\nStep1: import the gene names.')
    gene_name_list = import_gene_name(args.name)

    print('\nStep2: retrieve the fasta sequences.')
    sys.stdout.flush()
    retrieve_seq(gene_name_list, args.input, args.output, args.unmapped)

    print("\nAll done!")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))



if __name__ == '__main__':
    main()
