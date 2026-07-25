#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.007
# Function: Evaluate family-aware smORF evidence with cached parallel execution.
# Input: smorf_cluster outputs, source ORF table, and one or more density tracks.
# Output: Unit, family, member, calibration, and evidence-summary tables.

"""Command-line entry point for family-aware smORF Ribo-seq evidence."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence
from pathlib import Path

from utils.ribo.ArgsParser import (
    args_print,
    file_check,
    now_time,
    result_print,
    step_print,
    title_print,
)
from utils.smorf.smorf_riboseq_family_pipeline import (
    run_family_riboseq_evidence,
)


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Evaluate Ribo-seq translation evidence for fixed smORF families. "
            "Multiple samples support separate, pooled, and hybrid analysis."
        ),
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-i",
        "--family-table",
        dest="family_table",
        required=True,
        type=str,
        help="smorf_cluster <prefix>.family.message.txt file.",
    )
    required_group.add_argument(
        "-M",
        "--family-members",
        dest="family_members",
        required=True,
        type=str,
        help="smorf_cluster <prefix>.family.members.txt file.",
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output prefix for evidence tables.",
    )

    family_group = parser.add_argument_group("Family arguments")
    family_group.add_argument(
        "--orf-source",
        dest="orf_source",
        default=None,
        type=str,
        help=(
            "Full ORF message table used as smorf_cluster input. Required "
            "when the compact family-members table contains alternative-start "
            "representatives not present in the family-primary table."
        ),
    )
    family_group.add_argument(
        "--positive-category",
        dest="positive_category",
        default="annotated_ORF",
        type=str,
        help=(
            "Category used as positive controls for library calibration "
            "(default: %(default)s)."
        ),
    )

    density_group = parser.add_argument_group("P-site density arguments")
    density_group.add_argument(
        "-l",
        "--density-list",
        dest="density_list",
        default=None,
        type=str,
        help=(
            "Density list with sample/path and optional strand/format/group "
            "columns."
        ),
    )
    density_group.add_argument(
        "-d",
        "--density",
        dest="density",
        action="append",
        default=None,
        type=str,
        help="Direct density file; repeat for a plus/minus pair.",
    )
    density_group.add_argument(
        "-s",
        "--sample",
        dest="sample",
        default="sample1",
        type=str,
        help="Sample name for direct density input (default: %(default)s).",
    )
    density_group.add_argument(
        "--group",
        dest="group",
        default=None,
        type=str,
        help="Group name for direct density input (default: sample name).",
    )
    density_group.add_argument(
        "-f",
        "--density-format",
        dest="density_format",
        choices=["auto", "wig", "bedgraph"],
        default="auto",
        help="Direct density format (default: %(default)s).",
    )
    density_group.add_argument(
        "--analysis-mode",
        "--aggregation-mode",
        dest="analysis_mode",
        choices=[
            "auto",
            "single",
            "separate",
            "individual",
            "pooled",
            "merge",
            "merged",
            "hybrid",
        ],
        default="auto",
        help=(
            "Evidence aggregation mode. auto selects single for one sample "
            "and hybrid for multiple samples (default: %(default)s)."
        ),
    )

    annotation_group = parser.add_argument_group("Annotation arguments")
    annotation_group.add_argument(
        "-g",
        "--genepred",
        dest="genepred",
        default=None,
        type=str,
        help=(
            "Original transcript genePred used by smorf_scanner. It enables "
            "transcript-spliced post-stop release analysis."
        ),
    )
    annotation_group.add_argument(
        "-c",
        "--coord-mode",
        dest="coord_mode",
        choices=["0based-half-open", "1based-closed"],
        default="0based-half-open",
        type=str,
        help="ORF and genePred coordinate convention (default: %(default)s).",
    )
    annotation_group.add_argument(
        "-r",
        "--post-stop-codons",
        dest="post_stop_codons",
        default=10,
        type=int,
        help="Post-stop codons used for release analysis (default: %(default)s).",
    )

    evidence_group = parser.add_argument_group("Evidence arguments")
    evidence_group.add_argument(
        "--evidence-mode",
        dest="evidence_mode",
        choices=["sensitive", "balanced", "strict"],
        default="balanced",
        help="Base evidence preset before calibration (default: %(default)s).",
    )
    evidence_group.add_argument(
        "--min-rpf-sum",
        dest="min_rpf_sum",
        default=None,
        type=float,
        help="Override minimum coding-region P-site sum.",
    )
    evidence_group.add_argument(
        "--min-covered-codon",
        dest="min_covered_codon",
        default=None,
        type=int,
        help="Override minimum covered codon count.",
    )
    evidence_group.add_argument(
        "--min-codon-coverage",
        dest="min_codon_coverage",
        default=None,
        type=float,
        help="Override minimum whole-ORF codon coverage fraction.",
    )
    evidence_group.add_argument(
        "--moderate-periodicity",
        dest="moderate_periodicity",
        default=None,
        type=float,
        help="Override moderate frame-0 ratio.",
    )
    evidence_group.add_argument(
        "--strong-periodicity",
        dest="strong_periodicity",
        default=None,
        type=float,
        help="Override strong frame-0 ratio.",
    )

    adaptive_group = parser.add_argument_group("Length-adaptive arguments")
    adaptive_group.add_argument(
        "--short-max-codons",
        dest="short_max_codons",
        default=30,
        type=int,
        help="Maximum codon length using the short-ORF model (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--long-min-codons",
        dest="long_min_codons",
        default=100,
        type=int,
        help="Minimum codon length using the long-ORF model (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--window-codons",
        dest="window_codons",
        default=20,
        type=int,
        help="Sliding-window size in codons (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--window-step-codons",
        dest="window_step_codons",
        default=5,
        type=int,
        help="Sliding-window step in codons (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--min-window-rpf",
        dest="min_window_rpf",
        default=3.0,
        type=float,
        help="Minimum P-site sum for a supported window (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--min-window-covered-codons",
        dest="min_window_covered_codons",
        default=3,
        type=int,
        help="Minimum covered codons per supported window (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--min-supported-windows",
        dest="min_supported_windows",
        default=2,
        type=int,
        help="Distributed windows required for long ORFs (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--min-window-gap-codons",
        dest="min_window_gap_codons",
        default=15,
        type=int,
        help="Minimum start separation between independent windows (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--min-coverage-span",
        dest="min_coverage_span",
        default=0.35,
        type=float,
        help="Minimum supported span for whole long-ORF evidence (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--localized-span-max",
        dest="localized_span_max",
        default=0.20,
        type=float,
        help="Maximum span considered localized-only support (default: %(default)s).",
    )
    adaptive_group.add_argument(
        "--localized-top-window-fraction",
        dest="localized_top_window_fraction",
        default=0.70,
        type=float,
        help="Single-window RPF fraction defining localized signal (default: %(default)s).",
    )

    calibration_group = parser.add_argument_group("Positive-control calibration")
    calibration_group.add_argument(
        "--positive-quantile",
        dest="positive_quantile",
        default=0.20,
        type=float,
        help="Lower annotated-mORF quantile used for calibration (default: %(default)s).",
    )
    calibration_group.add_argument(
        "--positive-min-controls",
        dest="positive_min_controls",
        default=30,
        type=int,
        help="Minimum eligible annotated mORFs for calibration (default: %(default)s).",
    )
    calibration_group.add_argument(
        "--positive-min-rpf-sum",
        dest="positive_min_rpf_sum",
        default=10.0,
        type=float,
        help="Minimum annotated-mORF RPF sum used for calibration (default: %(default)s).",
    )
    calibration_group.add_argument(
        "--positive-max-controls",
        dest="positive_max_controls",
        default=5000,
        type=int,
        help=(
            "Maximum annotated mORFs sampled for threshold calibration "
            "(default: %(default)s)."
        ),
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--unit-member-mode",
        dest="unit_member_mode",
        choices=["family", "all", "none"],
        default="family",
        help=(
            "Write unit-level member evidence for multi-start families only, "
            "all families, or none (default: %(default)s)."
        ),
    )

    runtime_group = parser.add_argument_group("Runtime arguments")
    runtime_group.add_argument(
        "-t",
        "--thread",
        dest="threads",
        default=1,
        type=int,
        help=(
            "Maximum processes used for density caching, calibration, and "
            "chromosome evidence evaluation (default: %(default)s)."
        ),
    )
    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    files = [args.family_table, args.family_members]
    for path in (args.orf_source, args.density_list, args.genepred):
        if path is not None:
            files.append(path)
    files.extend(args.density or [])
    file_check(*files)

    if args.density_list and args.density:
        raise ValueError("--density-list cannot be combined with --density.")
    if not args.density_list and not args.density:
        raise ValueError("Provide --density-list or at least one --density file.")
    if args.group is None:
        args.group = args.sample
    if not str(args.sample).strip() or not str(args.group).strip():
        raise ValueError("--sample and --group must not be empty.")
    if args.post_stop_codons < 0:
        raise ValueError("--post-stop-codons must be >= 0.")
    if args.threads < 1:
        raise ValueError("--thread must be >= 1.")
    if args.short_max_codons < 3:
        raise ValueError("--short-max-codons must be >= 3.")
    if args.long_min_codons <= args.short_max_codons:
        raise ValueError("--long-min-codons must exceed --short-max-codons.")
    if args.window_codons < 6 or args.window_step_codons < 1:
        raise ValueError("Window size must be >= 6 and step must be >= 1.")
    if args.min_supported_windows < 1:
        raise ValueError("--min-supported-windows must be >= 1.")
    if args.positive_max_controls < args.positive_min_controls:
        raise ValueError(
            "--positive-max-controls must be >= --positive-min-controls."
        )
    for name in (
        "min_codon_coverage",
        "moderate_periodicity",
        "strong_periodicity",
        "min_coverage_span",
        "localized_span_max",
        "localized_top_window_fraction",
        "positive_quantile",
    ):
        value = getattr(args, name)
        if value is not None and not 0 <= value <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be in [0, 1].")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        _validate_args(args)
    except ValueError as error:
        parser.error(str(error))
    args_print(args)
    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Run family-aware smORF Ribo-seq evidence analysis."""
    now_time()
    title_print("Evaluate family-aware smORF translation evidence.")

    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)

    step_print(2, "Evaluate adaptive Ribo-seq evidence for smORF families.")
    result = run_family_riboseq_evidence(args)

    step_print(3, "Report family-aware evidence outputs.")
    result_print(
        [
            ("Analysis mode", result.mode),
            ("Samples", f"{result.sample_count:,}"),
            ("Families", f"{result.family_count:,}"),
            ("Representatives", f"{result.representative_count:,}"),
        ]
    )
    for output in result.outputs:
        result_print([("Output", output)])

    title_print("All done.")
    now_time()


if __name__ == "__main__":
    main()
