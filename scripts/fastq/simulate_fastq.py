#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : RiboParser
@Script       : simulate_fastq.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2025/07/29 21:02:35
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import argparse
import random

def parse_args():
    parser = argparse.ArgumentParser(description="Generate simulated FASTQ reads with random sequence and quality")
    parser.add_argument("-q", "--qmin", type=int, required=True,
                        help="Minimum quality score (e.g. 30)")
    parser.add_argument("-Q", "--qmax", type=int, required=True,
                        help="Maximum quality score (e.g. 45)")
    parser.add_argument("-m", "--minlen", type=int, required=True,
                        help="Minimum read length (e.g. 75)")
    parser.add_argument("-M", "--maxlen", type=int, required=True,
                        help="Maximum read length (e.g. 150)")
    parser.add_argument("-n", "--number", type=int, required=True,
                        help="Number of reads to generate")
    parser.add_argument("-o", "--output", default="simulated.fq",
                        help="Output FASTQ file name (default: simulated.fq)")
    return parser.parse_args()

def random_seq(length):
    """Generate random DNA sequence of given length"""
    return ''.join(random.choices("ACGT", k=length))

def random_qual(length, qmin, qmax):
    """Generate random Phred+33 quality string"""
    return ''.join(chr(random.randint(qmin, qmax) + 33) for _ in range(length))

def main():
    args = parse_args()

    if args.qmin > args.qmax:
        raise ValueError("qmin must be less than or equal to qmax")
    if args.minlen > args.maxlen:
        raise ValueError("minlen must be less than or equal to maxlen")

    with open(args.output, "w") as out:
        for i in range(1, args.number + 1):
            length = random.randint(args.minlen, args.maxlen)
            seq = random_seq(length)
            qual = random_qual(length, args.qmin, args.qmax)
            out.write(f"@read_{i}\n{seq}\n+\n{qual}\n")

    print(f"Generated {args.number} reads in: {args.output}")

if __name__ == "__main__":
    main()
