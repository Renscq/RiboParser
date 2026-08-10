#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-08
# Version: dev003
# Function: Call SeRP peaks with legacy or replicate-consensus algorithms without drawing figures.
# Input: TXT/JSON RPF density data, control/IP samples, and optional annotation/normalization.
# Output: Peak tables, BED files, peak sequences, and reusable enrichment profiles.

"""Command-line entry point for SeRP peak calling."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo.ArgsParser import (
    args_print,
    complete_print,
    file_check,
    now_time,
    step_print,
    title_print,
)
from utils.serp.SeRP import SeRP


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Call selective ribosome profiling peaks from control and IP RPF "
            "density profiles using legacy or replicate-consensus algorithms. "
            "TXT and JSON/JSONL density formats are supported."
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
            "Comma-separated control sample names. In consensus mode, sample order "
            "defines biological replicate pairing with --ip."
        ),
    )
    required_group.add_argument(
        "--ip",
        dest="ip",
        required=True,
        type=str,
        help=(
            "Comma-separated immunoprecipitation sample names. In consensus mode, "
            "sample order defines biological replicate pairing with --ck."
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
        "-a",
        "--annotation",
        dest="anno",
        default=None,
        type=str,
        help="Optional transcript annotation table used to retrieve gene names.",
    )

    method_group = parser.add_argument_group("Peak caller selection")
    method_group.add_argument(
        "--method",
        dest="method",
        choices=["legacy", "consensus"],
        default="legacy",
        help=(
            "Peak-calling method. legacy uses the corrected historical RiboParser "
            "algorithm; consensus uses matched biological replicates and peak-overlap support."
        ),
    )

    filtering_group = parser.add_argument_group("Data filtering arguments")
    filtering_group.add_argument(
        "-m",
        "--min",
        dest="min",
        default=50,
        type=int,
        help="Minimum gene-level RPF count required in every sample.",
    )
    filtering_group.add_argument(
        "--corr",
        dest="corr",
        default=None,
        type=float,
        help=(
            "Minimum pairwise replicate correlation. Method defaults are 0.3 for "
            "legacy and 0.5 for consensus."
        ),
    )
    filtering_group.add_argument(
        "--scale",
        dest="scale",
        default=1e6,
        type=float,
        help="Normalization scale used to calculate RPM-like abundance.",
    )

    background_group = parser.add_argument_group("Background arguments")
    background_group.add_argument(
        "-f",
        "--fill",
        dest="fill",
        default=30,
        type=int,
        help=(
            "Legacy-only control zero-fill mode: 0 uses global CDS background, "
            "-1 uses current-gene CDS mean, and positive values use a leading-codon window."
        ),
    )
    background_group.add_argument(
        "--back",
        dest="background",
        default=0,
        type=int,
        help=(
            "Leading CDS codons treated as background/non-callable region. Legacy "
            "mode supports 0 or 30; consensus mode accepts any non-negative value."
        ),
    )
    background_group.add_argument(
        "--bf",
        dest="back_fold",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Legacy-only switch for maximum-background fold filtering.",
    )

    peak_group = parser.add_argument_group("Peak calling arguments")
    peak_group.add_argument(
        "-s",
        "--smooth-window",
        dest="size",
        default=3,
        type=int,
        help="Legacy-only Savitzky-Golay smoothing window; use 0 to disable smoothing.",
    )
    peak_group.add_argument(
        "-k",
        "--polyorder",
        dest="k",
        default=1,
        type=int,
        help="Legacy-only Savitzky-Golay polynomial order.",
    )
    peak_group.add_argument(
        "-w",
        "--width",
        dest="width",
        default=5,
        type=int,
        help="Minimum binding-peak width in amino acids.",
    )
    peak_group.add_argument(
        "-e",
        "--enrich",
        dest="enrich",
        default=2.0,
        type=float,
        help="Strong-binding enrichment threshold.",
    )
    peak_group.add_argument(
        "-c",
        "--edge",
        "--collision",
        dest="collision",
        default=1.5,
        type=float,
        help=(
            "Lower enrichment threshold used to bridge/extend peak edges. "
            "--collision is retained as a compatibility alias."
        ),
    )
    peak_group.add_argument(
        "-g",
        "--gaps",
        dest="gaps",
        default=1,
        type=int,
        help="Maximum consecutive gap length retained within a candidate peak.",
    )
    peak_group.add_argument(
        "-p",
        "--proportion",
        dest="proportion",
        default=0.2,
        type=float,
        help="Maximum gap proportion retained within a candidate peak.",
    )
    peak_group.add_argument(
        "--all",
        dest="keep_all",
        action="store_true",
        default=False,
        help="Legacy-only: retain all qualified peak-region permutations instead of the best set.",
    )

    consensus_group = parser.add_argument_group("Consensus peak calling arguments")
    consensus_group.add_argument(
        "--consensus-window",
        dest="consensus_window",
        default=5,
        type=int,
        help="Centered rolling-sum window in codons for each matched replicate pair.",
    )
    consensus_group.add_argument(
        "--pseudocount",
        dest="pseudocount",
        default=0.1,
        type=float,
        help="RPM pseudocount added before local IP/control ratio calculation.",
    )
    consensus_group.add_argument(
        "--min-support",
        dest="min_support",
        default=1.0,
        type=float,
        help="Minimum fraction of matched replicate pairs supporting each consensus peak.",
    )
    consensus_group.add_argument(
        "--min-overlap",
        dest="min_overlap",
        default=0,
        type=int,
        help="Minimum replicate/consensus peak overlap in codons; 0 uses ceil(width/2).",
    )
    consensus_group.add_argument(
        "--stop-trim",
        dest="stop_trim",
        default=5,
        type=int,
        help="Number of terminal CDS codons excluded from consensus QC and peak calling.",
    )
    consensus_group.add_argument(
        "--min-codon-rpf",
        dest="min_codon_rpf",
        default=0.0,
        type=float,
        help=(
            "Consensus-only minimum mean raw RPF count per analyzed CDS codon "
            "required in every CK/IP sample; 0 disables this filter."
        ),
    )
    consensus_group.add_argument(
        "--max-edge-extension",
        dest="max_edge_extension",
        default=10,
        type=int,
        help="Consensus-only maximum lower-threshold edge extension per peak side in codons.",
    )

    sequence_group = parser.add_argument_group("Peak sequence arguments")
    sequence_group.add_argument(
        "--up",
        dest="upstream",
        default=10,
        type=int,
        help="Upstream codons retained around each peak sequence.",
    )
    sequence_group.add_argument(
        "--down",
        dest="downstream",
        default=10,
        type=int,
        help="Downstream codons retained around each peak sequence.",
    )

    output_group = parser.add_argument_group("Additional output arguments")
    output_group.add_argument(
        "--ratio",
        dest="ratio",
        action="store_true",
        default=False,
        help=(
            "Also output method-specific replicate enrichment ratios: all-pairwise "
            "ratios for legacy mode or matched local ratios for consensus mode."
        ),
    )

    return parser


def _parse_sample_names(value: str) -> list[str]:
    """Parse comma-separated sample names."""
    return [item.strip() for item in value.split(",") if item.strip()]


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.norm:
        file_check(args.norm)
    if args.anno:
        file_check(args.anno)

    control_samples = _parse_sample_names(args.control)
    ip_samples = _parse_sample_names(args.ip)
    if not control_samples:
        raise ValueError("--ck must contain at least one sample name.")
    if not ip_samples:
        raise ValueError("--ip must contain at least one sample name.")
    if len(set(control_samples)) != len(control_samples):
        raise ValueError("--ck contains duplicated sample names.")
    if len(set(ip_samples)) != len(ip_samples):
        raise ValueError("--ip contains duplicated sample names.")
    overlap = sorted(set(control_samples).intersection(ip_samples))
    if overlap:
        raise ValueError(
            "Control and IP groups share sample(s): {samples}".format(
                samples=", ".join(overlap)
            )
        )

    if args.method == "consensus":
        if len(control_samples) != len(ip_samples):
            raise ValueError(
                "--method consensus requires equal numbers of --ck and --ip samples."
            )
        if len(control_samples) < 2:
            raise ValueError(
                "--method consensus requires at least two matched biological replicate pairs."
            )

    if args.corr is None:
        args.corr = 0.5 if args.method == "consensus" else 0.3

    if args.scale <= 0:
        raise ValueError("--scale must be > 0.")
    if args.min < 0:
        raise ValueError("--min must be >= 0.")
    if not -1 <= args.corr <= 1:
        raise ValueError("--corr must be in [-1, 1].")
    if args.fill < -1:
        raise ValueError("--fill must be -1, 0, or a positive integer.")
    if args.background < 0:
        raise ValueError("--back must be >= 0.")
    if args.method == "legacy" and args.background not in {0, 30}:
        raise ValueError("--method legacy supports --back 0 or --back 30 only.")
    if args.size < 0 or (args.size > 0 and args.size % 2 == 0):
        raise ValueError("--smooth-window must be 0 or a positive odd integer.")
    if args.size == 0 and args.k != 0:
        args.k = 0
    if args.size > 0 and (args.k < 0 or args.k >= args.size):
        raise ValueError("--polyorder must satisfy 0 <= k < smooth-window.")
    if args.width < 1:
        raise ValueError("--width must be >= 1.")
    if args.enrich <= 0:
        raise ValueError("--enrich must be > 0.")
    if args.collision <= 0:
        raise ValueError("--edge/--collision must be > 0.")
    if args.method == "consensus" and args.collision > args.enrich:
        raise ValueError("--edge must be <= --enrich for consensus peak calling.")
    if args.gaps < 0:
        raise ValueError("--gaps must be >= 0.")
    if not 0 <= args.proportion <= 1:
        raise ValueError("--proportion must be in [0, 1].")
    if args.consensus_window < 1 or args.consensus_window % 2 == 0:
        raise ValueError("--consensus-window must be a positive odd integer.")
    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be > 0.")
    if not 0 < args.min_support <= 1:
        raise ValueError("--min-support must be in (0, 1].")
    if args.min_overlap < 0:
        raise ValueError("--min-overlap must be >= 0.")
    if args.stop_trim < 0:
        raise ValueError("--stop-trim must be >= 0.")
    if args.min_codon_rpf < 0:
        raise ValueError("--min-codon-rpf must be >= 0.")
    if args.max_edge_extension < 0:
        raise ValueError("--max-edge-extension must be >= 0.")
    if args.upstream < 0 or args.downstream < 0:
        raise ValueError("--up and --down must be >= 0.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for ``serp_peak``."""
    now_time()
    title_print("Call selective ribosome profiling peaks.")
    step_print(1, "Checking the input arguments.")
    args = _parse_args(argv)

    workflow = SeRP(args)
    workflow.prepare_output()

    step_print(2, "Import the RPF density data.")   
    workflow.import_rpf()

    step_print(3, "Import gene annotation.")
    workflow.import_annotation()

    step_print(4, "Call selective ribosome profiling peaks.")
    workflow.call_peaks()

    step_print(5, "Write peak-calling results.")
    workflow.write_results()

    complete_print()
    now_time()


if __name__ == "__main__":
    main()