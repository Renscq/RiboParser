#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev002
# Function: Calculate and draw sample-first SeRP enrichment metaplots around TIS and TTS.
# Input: TXT/JSON RPF density data and matched control/IP sample groups.
# Output: Sample metaprofiles, paired enrichment profiles, summary tables, and figures.

"""Command-line entry point for sample-first SeRP enrichment metaplots."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence
from pathlib import Path

from utils.ribo.ArgsParser import (
    args_print,
    complete_print,
    file_check,
    now_time,
    result_print,
    step_print,
    title_print,
)
from utils.serp.Metaplot import SeRPMetaplot


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build one normalized TIS/TTS metaprofile per sample, calculate matched "
            "IP/control enrichment between sample metaprofiles, and summarize the "
            "enrichment across biological replicate pairs. TXT and JSON/JSONL RPF "
            "density formats are supported."
        ),
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density TXT, JSON, JSONL, or compressed JSON density file.",
    )
    required_group.add_argument(
        "--ck",
        dest="control",
        required=True,
        type=str,
        help=(
            "Comma-separated control sample names. Sample order defines biological "
            "replicate pairing with --ip."
        ),
    )
    required_group.add_argument(
        "--ip",
        dest="ip",
        required=True,
        type=str,
        help=(
            "Comma-separated IP sample names. Sample order defines biological "
            "replicate pairing with --ck."
        ),
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output file prefix.",
    )

    input_group = parser.add_argument_group("Optional input arguments")
    input_group.add_argument(
        "-n",
        "--norm",
        dest="norm",
        default=None,
        type=str,
        help="Optional two-column file containing sample-level total RPF counts.",
    )
    input_group.add_argument(
        "--label",
        dest="label",
        default=None,
        type=str,
        help="Optional comparison label shown above the metaplot figure.",
    )

    filtering_group = parser.add_argument_group("Filtering arguments")
    filtering_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=50,
        type=int,
        help=(
            "Minimum gene-level CDS RPF count required in every control/IP sample. "
            "The same retained transcript set is used for all sample metaprofiles."
        ),
    )
    filtering_group.add_argument(
        "--scale",
        dest="scale",
        default=1e6,
        type=float,
        help="Normalization scale used to calculate sample-level RPM-like density.",
    )

    enrichment_group = parser.add_argument_group("Enrichment arguments")
    enrichment_group.add_argument(
        "--window",
        dest="window",
        default=5,
        type=int,
        help=(
            "Centered rolling-sum window in codons applied to sample metaprofiles "
            "before matched IP/control ratio calculation."
        ),
    )
    enrichment_group.add_argument(
        "--pseudocount",
        dest="pseudocount",
        default=0.1,
        type=float,
        help="Normalized-density pseudocount used for sample-metaprofile ratios.",
    )
    enrichment_group.add_argument(
        "--aggregate",
        dest="aggregate",
        choices=["median", "mean"],
        default="median",
        help="Statistic used to combine enrichment curves across replicate pairs.",
    )

    window_group = parser.add_argument_group("Metaplot window arguments")
    window_group.add_argument(
        "--tis",
        dest="tis",
        default=150,
        type=int,
        help="Number of CDS codons displayed from the start codon (0 to TIS-1).",
    )
    window_group.add_argument(
        "--tts",
        dest="tts",
        default=150,
        type=int,
        help="Number of CDS codons displayed before and including the stop codon.",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--detail",
        dest="detail",
        action="store_true",
        default=False,
        help=(
            "Write the per-transcript normalized sample densities used to build "
            "sample metaprofiles as gzip-compressed text."
        ),
    )
    output_group.add_argument(
        "--show-replicates",
        dest="show_replicates",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show matched biological-replicate enrichment curves behind the aggregate curve.",
    )
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["pdf", "png", "both"],
        default="pdf",
        help="Metaplot figure output format.",
    )
    output_group.add_argument(
        "--font-size",
        dest="font_size",
        default=9.0,
        type=float,
        help="Base figure font size.",
    )
    output_group.add_argument(
        "--dpi",
        dest="dpi",
        default=300,
        type=int,
        help="PNG output resolution.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments before starting the workflow."""
    file_check(args.rpf)
    if args.norm:
        file_check(args.norm)

    output_parent = Path(args.output).expanduser().resolve().parent
    if not output_parent.exists():
        raise FileNotFoundError(
            "Output directory does not exist: {0}".format(output_parent)
        )
    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.scale <= 0:
        raise ValueError("--scale must be > 0.")
    if args.window < 1:
        raise ValueError("--window must be >= 1.")
    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be > 0.")
    if args.tis < 1:
        raise ValueError("--tis must be >= 1.")
    if args.tts < 1:
        raise ValueError("--tts must be >= 1.")
    if args.font_size <= 0:
        raise ValueError("--font-size must be > 0.")
    if args.dpi < 1:
        raise ValueError("--dpi must be >= 1.")

    control_samples = [item.strip() for item in args.control.split(",") if item.strip()]
    ip_samples = [item.strip() for item in args.ip.split(",") if item.strip()]
    if len(control_samples) != len(ip_samples):
        raise ValueError(
            "--ck and --ip must contain equal numbers of samples for matched "
            "biological replicate enrichment."
        )


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        _validate_args(args)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    args_print(args)
    return args


def _run_metaplot_pipeline(args: Namespace) -> None:
    """Run the complete sample-first SeRP enrichment metaplot workflow."""
    workflow = SeRPMetaplot(args)

    step_print(2, "Import the RPF density data.")
    workflow.import_rpf()

    step_print(3, "Select one common transcript set for all sample metaprofiles.")
    workflow.select_genes()

    step_print(4, "Build sample metaprofiles and calculate matched enrichment.")
    workflow.calculate_metaplot()

    step_print(5, "Write sample and enrichment metaplot tables.")
    workflow.write_tables()

    step_print(6, "Draw TIS and TTS enrichment metaplots.")
    workflow.draw_metaplot()

    result_print(workflow.result_counts())


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for ``serp_metaplot``."""
    now_time()
    title_print("Draw selective ribosome profiling enrichment metaplots.")

    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)
    _run_metaplot_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
