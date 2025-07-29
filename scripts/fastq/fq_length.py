#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : RiboParser
@Script       : fq_length.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2025/07/29 20:57:55
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import argparse
import gzip
import os
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Count read lengths from FASTQ/FASTQ.GZ files")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="One or more input FASTQ or FASTQ.GZ files")
    parser.add_argument("-o", "--output", default="fq_length_distr.txt",
                        help="Output file name (default: fq_length_distr.txt)")
    return parser.parse_args()

def open_file(filename):
    """Open a file in text mode, using gzip if needed"""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def get_sample_name(filepath):
    """Extract sample name by removing path and file extension"""
    basename = os.path.basename(filepath)
    for ext in [".fq.gz", ".fastq.gz", ".fq", ".fastq"]:
        if basename.endswith(ext):
            return basename.replace(ext, "")
    return basename

def count_read_lengths(fq_path):
    """Count number of reads by length in a FASTQ file"""
    length_counter = defaultdict(int)
    with open_file(fq_path) as f:
        line_count = 0
        for line in f:
            line_count += 1
            if line_count % 4 == 2:  # sequence line
                read_len = len(line.strip())
                length_counter[read_len] += 1
    return length_counter

def main():
    args = parse_args()
    output_lines = []

    for fq_file in args.input:
        sample_name = get_sample_name(fq_file)
        print(f"Processing sample: {sample_name} ({fq_file})")
        len_counts = count_read_lengths(fq_file)
        for length, count in sorted(len_counts.items()):
            output_lines.append(f"{sample_name}\t{length}\t{count}")

    with open(args.output, "w") as out:
        out.write("Sample\tLength\tRead_Count\n")
        out.write("\n".join(output_lines))
    print(f"\nFinished. Results saved to: {args.output}")

if __name__ == "__main__":
    main()
