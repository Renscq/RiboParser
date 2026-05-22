#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : phred_quality.py


import sys
import argparse
import gc
import gzip

scType = [
    "Unknown data score format!",
    "Sanger | Phred+33 | Qual[33,73] | Val(0,40)",
    "Solexa | Solexa+64 | Qual[59,104] | Val(-5,40)",
    "Illumina1.3 | Phred+64 | Qual[64,104] | Val(0,40)",
    "Illumina1.5 | Phred+64 | Qual[66,104] | Val(3,40)",
    "Illumina1.8 | Phred+33 | Qual[33,74] | Val(0,41)"
]

def readFile(fq_read, result_list, num):
    flag = 1
    for line in fq_read:
        if flag > num * 4:
            break
        elif flag % 4 == 0:
            result_list.append(line.strip())
        flag += 1

def judgeFile(input_fq, num, result_list):
    basename = str(input_fq).split('.')[-1]
    if basename in ["fq", "fastq"]:
        with open(input_fq, 'r') as fq_read:
            readFile(fq_read, result_list, num)
    elif basename == "gz":
        with gzip.open(input_fq, 'rb') as fq_read:
            readFile(fq_read, result_list, num)
    else:
        print("Please input the FASTQ file!")
        sys.exit()

def checkFq(input_fq, result_list):
    score = []
    for line in result_list:
        for ABC in line:
            score.extend([str(ord(ABC))])
    MIN = int(min(score))
    MAX = int(max(score))

    if MIN < 33 or MAX > 104:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[0]))
    elif MIN >= 33 and MAX <= 73:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[1]))
    elif MIN >= 59 and MAX <= 104:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[2]))
    elif MIN >= 64 and MAX <= 104:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[3]))
    elif MIN >= 66 and MAX <= 104:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[4]))
    elif MIN >= 33 and MAX <= 74:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[5]))
    else:
        print ("%s [%s,%s] ==> %s" % (input_fq, MIN, MAX, scType[0]))


def args_parser():
    parser = argparse.ArgumentParser(description='This script is used to check the quality scores of FASTQ files.')
    parser.add_argument('-i', '--input', required=True, help='Specify the input FASTQ file (gzip format is OK).')
    parser.add_argument('-n', '--number', type=int, default=40000, help='Specify the number of sequences to extract (default: 40000).')
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    input_fq = args.input
    num = args.number
    result_list = []

    judgeFile(input_fq, num, result_list)
    checkFq(input_fq, result_list)
    gc.collect()

if __name__ == '__main__':
    main()
