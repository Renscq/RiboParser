#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Script  : getTissueFreqSeq.py
# @Time    : 2020/3/24 11:18
# @project : py-scripts
# @Version : python 3.7
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
from argparse import RawTextHelpFormatter
import gc
import time
import copy


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(description='This script is used to combine the tissue freq.',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    required = parser.add_argument_group(title='Required parameters')
    required.add_argument('-i', dest='input', required=True, type=str,
                          help='input file contain path and name of each sample')
    required.add_argument('-b', dest='barcode', required=False, type=str, default='all_barcodes.txt',
                          help='output file contain reads count of each sample (default: %(default)s)')
    required.add_argument('-o', dest='output', required=False, type=str, default='sRNA_tissue_freq_seq.fa',
                          help='output file contain non-redundant reads of each sample (default: %(default)s)')
    required.add_argument('-t', dest='txt', required=False, default=False, action='store_true',
                            help='output in the table format. (default: %(default)s)')
    args = parser.parse_args()
    return args


# global parameter
sample_list = {}
barcode = []
srna = {}


# read the list file
def list_read(file_list):
    print('Reading file list.')
    with open(file_list, 'r') as ids_input:
        num = 0
        for line in ids_input:
            pre_ids = line.strip().split('/')[-1].split('.')[0]
            flag = [(pre_ids == v.split('\t')[0]) for k, v in sample_list.items()]
            if True in flag:
                print('There is a repeats sample in your list:\n%s ' % pre_ids, flush=True)
                exit(1)
            else:
                sample_list[num] = pre_ids + '\t' + line.strip()
            barcode.append([pre_ids, 0])
            num += 1


# read srna file in fa format
def srna_read():
    total = len(sample_list.keys()) + 1
    temp_list = [0 for i in range(total)]
    
    for num, file in sample_list.items():
        names = file.split('\t')[1]
        print('Reading: %s' % names, flush=True)
        
        with open(names, 'r') as srna_input:
            for seq in srna_input:
                if seq.startswith('>'):
                    pass
                else:
                    seq = seq.strip()
                    if srna.get(seq):
                        srna[seq][num] += 1
                        barcode[num][1] += 1
                    else:
                        srna[seq] = copy.copy(temp_list)
                        srna[seq][num] += 1
                        barcode[num][1] += 1


# sum barcode
def sum_barcode():
    print('Sum total barcode.')
    for seq, num in srna.items():
        srna[seq][-1] = sum(num)


# output nn-redundant srna reads
def srna_write(output_file):
    print('Writing: %s' % output_file, flush=True)
    with open(output_file, 'w') as srna_output:
        num = 0
        for seq, counts in srna.items():
            num += 1
            counts = map(str, counts)
            srna_output.writelines('>t' + str(num) + ' ' + ' '.join(counts) + '\n' + seq.lower() + '\n')


# output nn-redundant srna reads
def barcode_write(barcode_file):
    print('Writing: %s' % barcode_file, flush=True)
    with open(barcode_file, 'w') as barcode_output:
        sum_total = 0
        for line in barcode:
            sum_total += line[1]
            barcode_output.writelines(line[0] + '\t' + str(line[1]) + '\n')
        barcode_output.writelines('total\t' + str(sum_total) + '\n')


# get time
def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)


# main programme
def main():
    now_time()
    args = parse_args()
    list_read(args.i)
    srna_read()
    sum_barcode()
    srna_write(args.o)
    barcode_write(args.b)
    gc.collect()
    now_time()


# run programme here
if __name__ == '__main__':
    main()
