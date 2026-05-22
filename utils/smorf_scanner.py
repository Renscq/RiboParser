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

from utils.smorf import SmORFPipeline


def parse_args():
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Scan transcript-centric smORFs from genome FASTA and genePred annotation."
    )

    parser.add_argument(
        "-g",
        "--genome",
        required=True,
        help="Input genome FASTA file."
    )

    parser.add_argument(
        "-a",
        "--annotation",
        required=True,
        help="Input genePred annotation file."
    )

    parser.add_argument(
        "-o",
        "--out-prefix",
        default="ORF",
        help="Output prefix."
    )

    parser.add_argument(
        "--orf-prefix",
        default="ORF",
        help="Prefix for ORF IDs."
    )

    parser.add_argument(
        "--start-codons",
        default="ATG",
        help="Comma-separated start codons, such as ATG,CTG,GTG,TTG."
    )

    parser.add_argument(
        "--min-aa",
        type=int,
        default=8,
        help="Minimum ORF length in amino acids."
    )

    parser.add_argument(
        "--max-aa",
        type=int,
        default=10000,
        help="Maximum ORF length in amino acids."
    )

    parser.add_argument(
        "--scan-strand",
        choices=["sense", "antisense", "both"],
        default="sense",
        help="Scan sense, antisense, or both strands."
    )

    parser.add_argument(
        "--kozak-up",
        type=int,
        default=6,
        help="Number of upstream nucleotides for Kozak sequence."
    )

    parser.add_argument(
        "--kozak-down",
        type=int,
        default=6,
        help="Number of downstream nucleotides after start codon for Kozak sequence."
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of worker processes for parallel ORF scanning."
    )

    parser.add_argument(
        "--mark-overlap",
        action="store_true",
        help="Mark nested or overlapping ORFs."
    )

    parser.add_argument(
        "--remove-discarded",
        action="store_true",
        help="Remove same-frame internal ORFs."
    )

    parser.add_argument(
        "--include-stop",
        action="store_true",
        help="Keep stop codon symbol in peptide sequence."
    )

    return parser.parse_args()


def main():
    """
    Run smORF scanner from command-line arguments.
    """

    args = parse_args()

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

    pipeline.run()


if __name__ == "__main__":
    main()
