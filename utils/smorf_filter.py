#!/usr/bin/env python3
# Author: Rensc
# date: 2026-05-22

"""
Command-line interface for ORF filtering.

Default Kozak mode:
    annotated

This means the script builds a species-specific Kozak PWM from annotated_ORF
records in ORF.message.txt and then scores all ORFs using this PWM.
"""

import argparse
import os
import sys


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(CURRENT_DIR)

if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

from utils.smorf.smorf_filter import ORFFilter
from utils.smorf.smorf_kozak import (
    KozakPWM,
    BUILTIN_KOZAK_PWMS,
    BUILTIN_KOZAK_CONSENSUS,
)


def parse_args():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Filter ORF.message.txt using rule-based criteria and Kozak PWM scoring."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input ORF.message.txt file."
    )

    parser.add_argument(
        "-o",
        "--out-prefix",
        default="ORF.filtered",
        help="Output prefix. (default: %(default)s)."
    )

    parser.add_argument(
        "--keep-start-codons",
        default="ATG,CTG,GTG,TTG",
        help="Comma-separated start codons to keep. (default: %(default)s)."
    )

    parser.add_argument(
        "--min-aa",
        type=int,
        default=8,
        help="Minimum ORF peptide length. (default: %(default)s)."
    )

    parser.add_argument(
        "--max-aa",
        type=int,
        default=10000,
        help="Maximum ORF peptide length. (default: %(default)s)."
    )

    parser.add_argument(
        "--keep-categories",
        default="uORF,dORF,lncORF,iORF,emORF,overlap_uORF,overlap_dORF,other_ORF,annotated_ORF",
        help="Comma-separated ORF categories to keep. (default: %(default)s)."
    )

    parser.add_argument(
        "--remove-categories",
        default="same_frame_iORF,antisense_ORF",
        help="Comma-separated ORF categories to remove. (default: %(default)s)."
    )

    parser.add_argument(
        "--keep-antisense",
        action="store_true",
        help="Keep antisense ORFs. (default: False)."
    )

    parser.add_argument(
        "--keep-secondary",
        action="store_true",
        help="Keep secondary ORFs. (default: False)."
    )

    parser.add_argument(
        "--keep-partial",
        action="store_true",
        help="Keep incomplete ORFs. (default: False)."
    )

    parser.add_argument(
        "--kozak-mode",
        choices=["none", "annotated", "builtin", "pwm", "sequence"],
        default="annotated",
        help="Kozak PWM mode. (default: %(default)s)."
    )

    parser.add_argument(
        "--builtin-kozak",
        default="plant",
        choices=sorted(BUILTIN_KOZAK_PWMS.keys()),
        help="Built-in Kozak PWM name. (default: %(default)s)."
    )

    parser.add_argument(
        "--kozak-pwm",
        default=None,
        help="Custom Kozak PWM matrix file. (default: %(default)s)."
    )

    parser.add_argument(
        "--kozak-seq",
        default=None,
        help="Aligned Kozak sequence file used to build PWM. (default: %(default)s)."
    )

    parser.add_argument(
        "--annotated-categories",
        default="annotated_ORF",
        help="Categories used to build annotated ORF Kozak PWM. (default: %(default)s)."
    )

    parser.add_argument(
        "--min-annotated-kozak",
        type=int,
        default=100,
        help="Minimum annotated ORF Kozak sequences required to build PWM. (default: %(default)s)."
    )

    parser.add_argument(
        "--fallback-builtin-kozak",
        default="plant",
        choices=sorted(BUILTIN_KOZAK_PWMS.keys()),
        help="Fallback built-in Kozak PWM if annotated mode fails. (default: %(default)s)."
    )

    parser.add_argument(
        "--no-kozak-fallback",
        action="store_true",
        help="Do not fallback to built-in PWM if annotated PWM construction fails. (default: False)."
    )

    parser.add_argument(
        "--min-kozak-pwm-score",
        type=float,
        default=0.0,
        help="Minimum normalized Kozak PWM score. Range: 0-1. (default: %(default)s)."
    )

    parser.add_argument(
        "--export-kozak-pwm",
        default=None,
        help="Export loaded or constructed Kozak PWM. (default: %(default)s)."
    )
    
    parser.add_argument(
        "--list-builtin-kozak",
        action="store_true",
        help="List built-in Kozak PWM models and exit."
    )

    return parser.parse_args()


def print_builtin_kozak_models():
    """
    Print available built-in Kozak models.
    """

    print("name\tconsensus\tlength")

    for name in sorted(BUILTIN_KOZAK_PWMS.keys()):
        pwm = BUILTIN_KOZAK_PWMS[name]
        consensus = BUILTIN_KOZAK_CONSENSUS.get(name, "NA")
        print("{}\t{}\t{}".format(name, consensus, len(pwm)))


def load_kozak_pwm(args, records):
    """
    Load or build Kozak PWM.
    """

    if args.kozak_mode == "none":
        return None

    if args.kozak_mode == "builtin":
        return KozakPWM.from_builtin(args.builtin_kozak)

    if args.kozak_mode == "pwm":
        if args.kozak_pwm is None:
            raise ValueError("--kozak-pwm is required when --kozak-mode pwm")

        return KozakPWM.from_pwm_file(
            args.kozak_pwm,
            name="custom_pwm",
        )

    if args.kozak_mode == "sequence":
        if args.kozak_seq is None:
            raise ValueError("--kozak-seq is required when --kozak-mode sequence")

        sequences = []

        with open(args.kozak_seq, "r") as handle:
            for line in handle:
                line = line.strip()

                if not line or line.startswith(">"):
                    continue

                sequences.append(line)

        return KozakPWM.from_sequences(
            sequences=sequences,
            name="custom_sequence_pwm",
        )

    if args.kozak_mode == "annotated":
        try:
            return KozakPWM.from_annotated_orf_records(
                records=records,
                annotated_categories=args.annotated_categories,
                min_sequences=args.min_annotated_kozak,
                name="annotated_ORF_pwm",
            )
        except ValueError as error:
            if args.no_kozak_fallback:
                raise error

            print(
                "[smORFFilter] Warning: {}. Fallback to built-in Kozak PWM: {}".format(
                    error,
                    args.fallback_builtin_kozak,
                ),
                flush=True,
            )

            return KozakPWM.from_builtin(args.fallback_builtin_kozak)

    return None


def main():
    """
    Run ORF filtering.
    """

    args = parse_args()

    if args.list_builtin_kozak and not args.input and not args.out_prefix:
        print('\nStep1: Loading Kozak PWM', flush=True)
        print_builtin_kozak_models()
        return

    print('\nStep1: Reading ORF records', flush=True)
    header, records = ORFFilter.read_table(args.input)

    print('\nStep2: Loading Kozak PWM', flush=True)
    kozak_pwm = load_kozak_pwm(args, records)

    if kozak_pwm is not None and args.export_kozak_pwm is not None:
        kozak_pwm.export_pwm(args.export_kozak_pwm)

    print('\nStep3: Filtering ORF records', flush=True)
    orf_filter = ORFFilter(
        keep_start_codons=args.keep_start_codons,
        min_aa=args.min_aa,
        max_aa=args.max_aa,
        min_kozak_pwm_score=args.min_kozak_pwm_score,
        kozak_pwm=kozak_pwm,
        keep_categories=args.keep_categories,
        remove_categories=args.remove_categories,
        require_sense=not args.keep_antisense,
        require_primary=not args.keep_secondary,
        require_complete=not args.keep_partial,
    )

    
    header = orf_filter.add_kozak_pwm_fields(header, records)

    passed, removed = orf_filter.filter_records(records)

    print('\nStep4: Format the output', flush=True)
    passed_path = "{}.passed.message.txt".format(args.out_prefix)
    removed_path = "{}.removed.message.txt".format(args.out_prefix)
    all_path = "{}.all.message.txt".format(args.out_prefix)

    print('\nStep5: Writing output files', flush=True)
    ORFFilter.write_table(passed_path, header, passed)
    ORFFilter.write_table(removed_path, header, removed)
    ORFFilter.write_table(all_path, header, records)

    print("[smORFFilter] Total ORFs: {}".format(len(records)), flush=True)
    print("[smORFFilter] Passed ORFs: {}".format(len(passed)), flush=True)
    print("[smORFFilter] Removed ORFs: {}".format(len(removed)), flush=True)

    if kozak_pwm is not None:
        print("[smORFFilter] Kozak PWM name: {}".format(kozak_pwm.name), flush=True)
        print("[smORFFilter] Kozak PWM source: {}".format(kozak_pwm.source), flush=True)
        print("[smORFFilter] Kozak consensus: {}".format(kozak_pwm.consensus), flush=True)

    print("[smORFFilter] Output passed file: {}".format(passed_path), flush=True)


if __name__ == "__main__":
    main()
