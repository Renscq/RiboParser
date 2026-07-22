#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Integrate sample-specific smORF evidence tables.
# Input: Complete sample evidence tables and optional ORF genePred annotation.
# Output: Reduced integrated evidence, metric matrices, and filtered genePred.

"""Command-line entry point for simplified smORF evidence integration."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import (
    args_print,
    file_check,
    now_time,
    result_print,
    step_print,
    title_print,
)
from utils.smorf.smorf_riboseq_integrate import run_integration


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser.

    Returns:
        Configured parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Integrate complete sample-specific smORF evidence tables. "
            "ORFs marked NoEvidence in every sample are removed before "
            "full metric integration."
        )
    )

    input_group = parser.add_argument_group("Input arguments")
    input_group.add_argument(
        "-i",
        "--input",
        dest="input",
        nargs="+",
        required=True,
        type=str,
        help=(
            "One evidence table per sample. Paths and glob patterns are "
            "supported. Every file must contain the same complete ORF set."
        ),
    )
    input_group.add_argument(
        "--orf-genepred",
        dest="orf_genepred",
        default=None,
        type=str,
        help=(
            "Optional scanner ORF genePred keyed by orf_id. Only final "
            "retained ORFs are written."
        ),
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "-o",
        "--output-integrated",
        dest="output_integrated",
        default=None,
        type=str,
        help="Output compact ORF-level integrated evidence table.",
    )
    output_group.add_argument(
        "-m",
        "--output-matrix",
        dest="output_matrix",
        default=None,
        type=str,
        help=(
            "Matrix output template. The RPF-sum matrix is written by "
            "default when this option is provided."
        ),
    )
    output_group.add_argument(
        "--output-density-matrix",
        dest="output_density_matrix",
        action="store_true",
        default=False,
        help="Also output frame0/frame1/frame2 density matrix.",
    )
    output_group.add_argument(
        "--output-ratio-matrix",
        dest="output_ratio_matrix",
        action="store_true",
        default=False,
        help="Also output frame0/frame1/frame2 ratio matrix.",
    )
    output_group.add_argument(
        "--output-covered-codon-ratio-matrix",
        dest="output_covered_codon_ratio_matrix",
        action="store_true",
        default=False,
        help="Also output covered-codon-ratio matrix.",
    )
    output_group.add_argument(
        "--output-translation-evidence-matrix",
        dest="output_translation_evidence_matrix",
        action="store_true",
        default=False,
        help="Also output translation-evidence matrix.",
    )

    integration_group = parser.add_argument_group("Integration arguments")
    integration_group.add_argument(
        "-c",
        "--capture-labels",
        dest="capture_labels",
        default="LowConfidence,MediumConfidence,HighConfidence",
        type=str,
        help=(
            "Labels counted as captured when rpf_sum is positive "
            "(default: %(default)s)."
        ),
    )
    integration_group.add_argument(
        "-p",
        "--pass-labels",
        dest="pass_labels",
        default="MediumConfidence,HighConfidence",
        type=str,
        help=(
            "Labels counted as reliable when rpf_sum is positive "
            "(default: %(default)s)."
        ),
    )
    integration_group.add_argument(
        "-n",
        "--excellent-min-samples",
        dest="excellent_min_samples",
        default=2,
        type=int,
        help=(
            "Reliable sample count defining Excellent "
            "(default: %(default)s)."
        ),
    )
    integration_group.add_argument(
        "--reliable-only",
        dest="reliable_only",
        action="store_true",
        default=False,
        help=(
            "After all-sample NoEvidence removal, retain only ORFs with at "
            "least one reliable sample."
        ),
    )

    runtime_group = parser.add_argument_group("Runtime arguments")
    runtime_group.add_argument(
        "-t",
        "--threads",
        dest="threads",
        default=4,
        type=int,
        help=(
            "Concurrent sample-file readers used in both input passes "
            "(default: %(default)s)."
        ),
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments.

    Args:
        args: Parsed arguments.

    Raises:
        ValueError: If an argument combination is invalid.
    """
    if not args.output_integrated and not args.output_matrix:
        raise ValueError(
            "At least one of --output-integrated or --output-matrix is required."
        )
    if args.excellent_min_samples < 1:
        raise ValueError("--excellent-min-samples must be >= 1.")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1.")
    optional_matrices = (
        args.output_density_matrix,
        args.output_ratio_matrix,
        args.output_covered_codon_ratio_matrix,
        args.output_translation_evidence_matrix,
    )
    if any(optional_matrices) and not args.output_matrix:
        raise ValueError(
            "Optional matrix switches require --output-matrix."
        )
    if args.orf_genepred is not None:
        file_check(args.orf_genepred)


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse and validate arguments.

    Args:
        argv: Optional arguments for tests.

    Returns:
        Validated arguments.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        _validate_args(args)
    except ValueError as error:
        parser.error(str(error))
    args_print(args)
    return args




def main(argv: Sequence[str] | None = None) -> None:
    """Run the smorf_integrate command.

    Args:
        argv: Optional arguments for tests.
    """
    now_time()
    title_print("Integrate smORF Ribo-seq evidence.")
    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)

    step_print(2, "Remove all-sample NoEvidence ORFs and integrate evidence.")
    result = run_integration(args)

    step_print(3, "Report integration outputs.")
    summary_items = [
        ("Input samples", f"{result.sample_count:,}"),
        ("Input ORFs", f"{result.input_orf_count:,}"),
        ("Evidence ORFs", f"{result.evidence_orf_count:,}"),
        ("Final ORFs", f"{result.retained_orf_count:,}"),
    ]
    if result.integrated_output:
        summary_items.append(("Integrated table", result.integrated_output))
    for index, matrix_output in enumerate(result.matrix_outputs, start=1):
        summary_items.append((f"Matrix table {index}", matrix_output))
    if result.genepred_output:
        summary_items.append(("Filtered genePred", result.genepred_output))
    result_print(summary_items)

    title_print("All done.")
    now_time()


if __name__ == "__main__":
    main()
