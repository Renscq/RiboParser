#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Script  : nt2aa.py
# @Version : python 3.8.5
# @E-mail  : rensc0718@163.com


import argparse
import gc
import re
import sys
import time
from collections import OrderedDict


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(description="This script is used to convert the nucleotide to amino acid.")
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-i", "--input", required=True, type=str, help="input cds sequences in fasta format"
    )
    parser.add_argument(
        "-t", "--tis", required=False, action='store_true', default=False,
        help="start with ATG (default %(default)s)"
    )
    parser.add_argument(
        "-o", "--output", required=False, type=str, default='results-aa.fa',
        help="output file name (default %(default)s)"
    )

    args = parser.parse_args()

    args_dict = vars(args)
    for k, v in args_dict.items():
        print(k + ": " + str(v))
    print()
    sys.stdout.flush()
    return args


codon_dict = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T',
              'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
              'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H',
              'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
              'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L',
              'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
              'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G',
              'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
              'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S',
              'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
              'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}


def nt2aa(nt_file_in, aa_file_out):
    gene = OrderedDict()
    with open(nt_file_in, 'r') as nt_file:
        for line in nt_file:

            if line.startswith('>'):
                gene_name = line.strip()
                gene[gene_name] = ''
            else:
                rec = line.strip().upper().replace('U', 'T')
                gene[gene_name] += rec

    with open(aa_file_out, 'w') as aa_file:
        for gene_name, seq in gene.items():
            aa = ''
            unexp_nt = re.sub(r'[ATGCU]', '', seq)

            if len(seq) % 3 != 0:
                print("length error:\t" + gene_name + '\t' + str(len(seq)))
                continue
            elif len(unexp_nt) > 0:
                print('nucleotide error:\t' + gene_name + '\t' + ''.join(set(list(unexp_nt))))
                continue
            else:
                for num in range(0, len(seq), 3):
                    nt = seq[num:num + 3]
                    aa += codon_dict[nt]
                aa_file.writelines('\n'.join([gene_name, aa]) + '\n')


def nt2aa_atg(nt_file_in, aa_file_out):
    gene = OrderedDict()
    with open(nt_file_in, 'r') as nt_file:
        for line in nt_file:

            if line.startswith('>'):
                gene_name = line.strip()
                gene[gene_name] = ''
            else:
                rec = line.strip().upper().replace('U', 'T')
                gene[gene_name] += rec

    with open(aa_file_out, 'w') as aa_file:
        for gene_name, seq in gene.items():
            aa = ''
            unexp_nt = re.sub(r'[ATGCU]', '', seq)
            stop_codon = seq[-3::]

            if len(seq) % 3 != 0:
                print("length error:\t" + gene_name + '\t' + str(len(seq)))
                continue
            elif not seq.startswith('ATG'):
                print('start codon error:\t' + gene_name + '\t' + seq[0:3])
                continue
            elif len(unexp_nt) > 0:
                print('nucleotide error:\t' + gene_name + '\t' + ''.join(set(list(unexp_nt))))
                continue
            elif stop_codon != 'TGA' and stop_codon != 'TAG' and stop_codon != 'TAA':
                print('stop codon warning:\t' + gene_name + '\t' + stop_codon)
                for num in range(0, len(seq), 3):
                    nt = seq[num:num + 3]
                    aa += codon_dict[nt]
                aa_file.writelines('\n'.join([gene_name, aa]) + '\n')
            else:
                for num in range(0, len(seq), 3):
                    nt = seq[num:num + 3]
                    aa += codon_dict[nt]
                aa_file.writelines('\n'.join([gene_name, aa]) + '\n')


# gc collect
def clear():
    gc.collect()
    print("\nAll done.")
    sys.stdout.flush()


# get time
def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sys.stdout.flush()


# main
def main():
    args = parse_args()
    if args.tis:
        nt2aa_atg(args.input, args.output)
    else:
        nt2aa(args.input, args.output)
    clear()
    now_time()


if __name__ == '__main__':
    main()
