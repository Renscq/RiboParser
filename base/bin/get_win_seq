#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Script  : get_win_seq.py
# @Time    : 2022/11/1
# @project : py-scripts
# @Version : python 3.7
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse


def get_argparse():
    parser = argparse.ArgumentParser(description='rolling the windows on each sequence')
    parser.add_argument("-i", "--input", dest="fasta", required=False, help="input fasta")
    parser.add_argument("-w", "--win", dest="width", type=int, default=20, required=False, help="width of window")
    parser.add_argument("-o", "--output", dest="txt", required=False, help="output txt")
    args = parser.parse_args()
    return args


def get_win(seq, flag, win_width, outfile):
    length = len(seq)
    if length > win_width:
        b = [seq[i: i + win_width] for i in range(length + 1 - win_width)]
        for items in b:
            outfile.write('\t'.join([str(flag), seq, items]) + '\n')
    else:
        outfile.write('\t'.join([str(flag), seq, seq]) + '\n')


def read_file(inputfile, win_width, out_name):
    if out_name:
        pass
    else:
        out_name = inputfile + '.win'

    outfile = open(out_name, 'w')

    with open(inputfile, 'r') as file_in:
        flag = 0
        for line in file_in:
            if line.startswith('>'):
                continue
            flag += 1
            get_win(line.strip(), flag, win_width, outfile)

    outfile.close()


if __name__ == '__main__':
    args = get_argparse()
    read_file(args.fasta, args.width, args.txt)
