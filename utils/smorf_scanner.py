#!/usr/bin/env python3
# Author: Rensc
# date: 2026-05-21

"""
Command-line interface for smORF scanning.

This script calls classes from utils/smorf and performs transcript-centric
ORF scanning using genome FASTA and genePred annotation.
"""

import argparse
import os
import sys


# Add the project root directory to Python import path.
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(CURRENT_DIR)

if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

from utils.ribo.ArgsParser import args_print, file_check, now_time
from utils.smorf import SmORFPipeline


def parse_args():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="This script is used to scan transcript-centric smORFs."
    )

    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-g', '--genome', dest="genome", required=True, type=str,
                             help="the input genome sequence file in FASTA format.")
    input_group.add_argument('-a', '--annotation', dest="annotation", required=True, type=str,
                             help="the input transcript annotation file in genePred format.")

    parser.add_argument('-o', '--out-prefix', dest="out_prefix", required=False, type=str, default="ORF",
                        help="the prefix of output files. (default: %(default)s).")
    parser.add_argument('-p', '--orf-prefix', dest="orf_prefix", required=False, type=str, default="ORF",
                        help="the prefix of generated ORF identifiers. (default: %(default)s).")
    parser.add_argument('-s', '--start-codons', dest="start_codons", required=False, type=str, default="ATG",
                        help="comma-separated start codons used for ORF scanning. (default: %(default)s).")
    parser.add_argument('-m', '--min-aa', dest="min_aa", required=False, type=int, default=8,
                        help="the minimum ORF length to keep. (default: %(default)s aa).")
    parser.add_argument('-M', '--max-aa', dest="max_aa", required=False, type=int, default=10000,
                        help="the maximum ORF length to keep. (default: %(default)s aa).")
    parser.add_argument('-x', '--scan-strand', dest="scan_strand", required=False, type=str,
                        choices=["sense", "antisense", "both"], default="sense",
                        help="specify which transcript strand will be scanned. (default: %(default)s).")
    parser.add_argument('-u', '--kozak-up', dest="kozak_up", required=False, type=int, default=6,
                        help="the upstream nucleotides extracted for Kozak sequence. (default: %(default)s nt).")
    parser.add_argument('-d', '--kozak-down', dest="kozak_down", required=False, type=int, default=6,
                        help="the downstream nucleotides after start codon extracted for Kozak sequence. (default: %(default)s nt).")
    parser.add_argument('-t', '--threads', dest="threads", required=False, type=int, default=1,
                        help="the number of worker processes for ORF scanning. (default: %(default)s).")
    parser.add_argument('-O', '--mark-overlap', dest="mark_overlap", action="store_true", required=False, default=False,
                        help="mark nested or overlapping ORFs after scanning. (default: %(default)s).")
    parser.add_argument('-R', '--remove-discarded', dest="remove_discarded", action="store_true", required=False, default=False,
                        help="remove same-frame internal ORFs from the final output. (default: %(default)s).")
    parser.add_argument('-I', '--include-stop', dest="include_stop", action="store_true", required=False, default=False,
                        help="keep the stop codon symbol in peptide sequence. (default: %(default)s).")

    args = parser.parse_args()

    return args


def main():
    """
    Run smORF scanner from command-line arguments.
    """

    args = parse_args()

    now_time()

    print('\nScan transcript-centric smORFs.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    file_check(args.genome, args.annotation)
    args_print(args)

    pipeline = SmORFPipeline(
        genome=args.genome,
        annotation=args.annotation,
        out_prefix=args.out_prefix,
        orf_prefix=args.orf_prefix,
        start_codons=args.start_codons,
        min_aa=args.min_aa,
        max_aa=args.max_aa,
        scan_strand=args.scan_strand,
        kozak_up=args.kozak_up,
        kozak_down=args.kozak_down,
        mark_overlap=args.mark_overlap,
        remove_discarded=args.remove_discarded,
        include_stop=args.include_stop,
        threads=args.threads,
    )

    print('\nStep2: Scan ORFs from transcript sequences.', flush=True)
    pipeline.run()

    print('\nAll done.', flush=True)
    now_time()


if __name__ == "__main__":
    main()
