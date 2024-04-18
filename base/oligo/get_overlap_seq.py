#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : get_overlap_seq.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2023/03/10 10:45:16
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import argparse
from Bio import SeqIO
from Bio import Seq

def args_parser():

    parser = argparse.ArgumentParser(description='This script is used to get overlap seq.')

    # needed arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        '-i', dest='input', required=True, type=str, help='input fasta file name.'
    )
    input_group.add_argument(
        '-o', dest='output', required=True, type=str, help='output file name.'
    )

    # optional arguments
    parser.add_argument(
        '-l', dest='length', required=False, type=int, default=6, help='nucleotide length (default: %(default)s nt).'
    )

    args = parser.parse_args()

    return args


def parse_fasta(fasta_file):
    # import the fasta file with biopython and return a dict
    sequences = {}

    record = SeqIO.parse(fasta_file, "fasta")
    for line in record:
        sequences[str(line.id)] = str(line.seq)

    return sequences


def search_overlap(sequences, nt_length):
    matched_sequences = []
    for id1, seq1 in sequences.items():
        for id2, seq2 in sequences.items():
            if id1 == id2:
                continue
            for i in range(nt_length, len(seq1) + 1):
                sub_seq = seq2[:i]
                if seq1.endswith(sub_seq) and len(sub_seq) > nt_length:
                    matched_sequences.append((id1, seq1, len(seq1), id2, seq2, len(seq2), sub_seq))
                    break

    return matched_sequences


def output_overlap_seq(matched_sequences, out_file):
    with open(out_file, 'w', newline='') as output:

        output.writelines('\t'.join(['id1', 'seq1', 'seq1_len', 'id2', 'seq2', 'seq2_len', 'overlap']) + '\n')
        for row in matched_sequences:
            now_row = [str(i) for i in row]
            output.writelines('\t'.join(list(now_row)) + '\n')


def main():
    args = args_parser()
    sequences = parse_fasta(args.input)
    matched_sequences = search_overlap(sequences, args.length)
    output_overlap_seq(matched_sequences, args.output)


if __name__ == '__main__':
    main()
