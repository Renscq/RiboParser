#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-11
# Version: 0.2.8-dev.005
# Function: Draw IGV-like gene-level RPF/RNA density profiles.
# Input: RPF/RNA density file in JSONL or TXT format, target ID(s), and optional genome annotation.
# Output: Transcript-coordinate and optional genome-coordinate gene-level density profiles and figures.

"""Command-line entry point for gene-level RPF density plotting."""

from __future__ import annotations

import argparse
import copy
import re
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Geneplot

from utils.ribo.ArgsParser import (
    args_print,
    complete_print,
    file_check,
    now_time,
    step_print,
    title_print,
)


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Draw IGV-like gene-level RPF/RNA density profiles from JSONL or TXT density files."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        "--rpf",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF/RNA density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-g",
        "--target",
        dest="target",
        required=False,
        type=str,
        default=None,
        help="Target gene ID or transcript ID to plot. Required unless --target-list is provided.",
    )
    required_group.add_argument(
        "--target-list",
        dest="target_list",
        required=False,
        type=str,
        default=None,
        help=(
            "Target list file for batch plotting. The first non-empty column is used. "
            "Lines beginning with # are ignored. Required unless --target/-g is provided."
        ),
    )
    required_group.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    selection_group = parser.add_argument_group("Target and sample selection")
    selection_group.add_argument(
        "--id-type",
        dest="id_type",
        choices=["auto", "gene", "transcript"],
        default="auto",
        help="How to match target IDs in JSON input. TXT input can only match the name column. Default: %(default)s.",
    )
    selection_group.add_argument(
        "--sample",
        dest="sample",
        required=False,
        default=None,
        type=str,
        help="Optional comma-separated sample names to plot. Default: all samples.",
    )
    selection_group.add_argument(
        "--select-transcript",
        dest="select_transcript",
        choices=["first", "longest", "highest"],
        default="longest",
        help="Transcript selection method when one target matches multiple JSON records. Default: %(default)s.",
    )

    coordinate_group = parser.add_argument_group("Coordinate and density options")
    coordinate_group.add_argument(
        "--coordinate",
        dest="coordinate",
        choices=["auto", "both", "genome", "transcript"],
        default="auto",
        help=(
            "Plot coordinate system. 'auto' draws transcript coordinates and additionally draws genome "
            "coordinates when embedded JSONL genome_mapping or --annotation is available. "
            "'both' has the same plotting behavior but is explicit. Default: %(default)s."
        ),
    )
    coordinate_group.add_argument(
        "--annotation",
        dest="annotation",
        required=False,
        type=str,
        default=None,
        help=(
            "Optional transcript annotation in genePred/genePredExt or RiboParser norm TXT format. "
            "This enables genome-coordinate plotting for TXT density files and can rescue JSONL files "
            "without embedded genome_mapping. Default: not used."
        ),
    )

    coordinate_group.add_argument(
        "--data-type",
        dest="data_type",
        choices=["ribo", "rna"],
        type=str.lower,
        default="ribo",
        help=(
            "Density data type. 'ribo' keeps frame-colored bars in --mode bar --frame all; "
            "'rna' draws all density bars with one color and automatically uses --frame all. Default: %(default)s."
        ),
    )
    coordinate_group.add_argument(
        "-f",
        "--frame",
        dest="frame",
        choices=["all", "0", "1", "2"],
        default="all",
        help="Frame density to plot. Bar mode with all shows three frame colors. Default: %(default)s.",
    )
    norm_group = coordinate_group.add_mutually_exclusive_group()
    norm_group.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        default=True,
        help="Normalize RPF/RNA counts to RPM before plotting. Default: True.",
    )
    norm_group.add_argument(
        "--raw-count",
        dest="normal",
        action="store_false",
        help="Plot raw RPF/RNA counts instead of RPM.",
    )

    plot_group = parser.add_argument_group("Density plot style")
    plot_group.add_argument(
        "--mode",
        dest="mode",
        choices=["line", "bar"],
        default="bar",
        help="Density plot style. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--plot-transform",
        dest="plot_transform",
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        default="none",
        help=(
            "Transform plotted density values without changing output tables. "
            "'log' is an alias of log1p. Default: %(default)s."
        ),
    )
    plot_group.add_argument(
        "--utr-gray",
        dest="utr_gray",
        action="store_true",
        default=False,
        help="Draw UTR density and UTR structure in gray. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--line-width",
        dest="line_width",
        type=float,
        default=1.2,
        help="Line width for line mode. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--bar-width",
        dest="bar_width",
        type=float,
        default=0.85,
        help="Bar width for bar mode. Default: %(default)s.",
    )
    plot_group.add_argument(
        "--y-max",
        dest="y_max",
        type=float,
        default=None,
        help="Optional shared y-axis maximum after plot transformation. Default: auto.",
    )

    clipping_group = parser.add_argument_group("Spike clipping options")
    clipping_group.add_argument(
        "--spike-clip",
        dest="spike_clip",
        action="store_true",
        default=False,
        help=(
            "Clip extreme plotted density values after plot transformation. "
            "Only the figure is clipped; profile.txt keeps original values. Default: %(default)s."
        ),
    )
    clipping_group.add_argument(
        "--clip-quantile",
        dest="clip_quantile",
        type=float,
        default=0.995,
        help="Quantile used as spike clipping cutoff when --spike-clip is enabled. Default: %(default)s.",
    )
    clipping_group.add_argument(
        "--clip-value",
        dest="clip_value",
        type=float,
        default=None,
        help="Fixed upper cutoff after plot transformation. Overrides --clip-quantile. Default: auto.",
    )

    structure_group = parser.add_argument_group("Gene structure display")
    structure_group.add_argument(
        "--intron-scale",
        dest="intron_scale",
        type=float,
        default=1.0,
        help="Scale factor for intron display length in genome mode. Use 1 for true genomic scale. Default: %(default)s.",
    )
    structure_group.add_argument(
        "--x-margin",
        dest="x_margin",
        type=float,
        default=0.02,
        help="Fractional x-axis margin on both sides. Default: %(default)s.",
    )

    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--export-track",
        dest="export_track",
        action="store_true",
        default=False,
        help="Export a multi-sample bedGraph-like table in addition to raw profile.txt. Default: %(default)s.",
    )
    output_group.add_argument(
        "--output-format",
        dest="output_format",
        choices=["pdf", "png", "both"],
        default="both",
        help="Figure output format. Default: %(default)s.",
    )
    output_group.add_argument(
        "--dpi",
        dest="dpi",
        type=int,
        default=300,
        help="PNG output DPI. Default: %(default)s.",
    )

    layout_group = parser.add_argument_group("Figure size and layout")
    layout_group.add_argument(
        "--width",
        dest="width",
        type=float,
        default=11.0,
        help="Figure width in inches. Default: %(default)s.",
    )
    layout_group.add_argument(
        "--per-sample-height",
        dest="per_sample_height",
        type=float,
        default=1.0,
        help="Figure height per sample panel in inches. Default: %(default)s.",
    )
    layout_group.add_argument(
        "--structure-height",
        dest="structure_height",
        type=float,
        default=0.45,
        help="Relative height of the gene structure panel. Default: %(default)s.",
    )
    layout_group.add_argument(
        "--max-height",
        dest="max_height",
        type=float,
        default=18.0,
        help="Maximum figure height in inches. Default: %(default)s.",
    )

    text_group = parser.add_argument_group("Text and font options")
    text_group.add_argument(
        "--title",
        dest="title",
        type=str,
        default=None,
        help="Optional figure title. Default: gene/transcript information.",
    )
    text_group.add_argument(
        "--font-size",
        dest="font_size",
        type=float,
        default=9.0,
        help="Base font size for geneplot figures. Default: %(default)s.",
    )
    text_group.add_argument(
        "--title-size",
        dest="title_size",
        type=float,
        default=None,
        help="Title font size. Default: --font-size + 2.",
    )
    text_group.add_argument(
        "--label-size",
        dest="label_size",
        type=float,
        default=None,
        help="Axis-label font size. Default: --font-size.",
    )
    text_group.add_argument(
        "--tick-size",
        dest="tick_size",
        type=float,
        default=None,
        help="Tick-label font size. Default: --font-size - 1.",
    )
    text_group.add_argument(
        "--sample-label-size",
        dest="sample_label_size",
        type=float,
        default=None,
        help="Right-side sample-label font size. Default: --font-size.",
    )
    text_group.add_argument(
        "--legend-size",
        dest="legend_size",
        type=float,
        default=None,
        help="Legend font size. Default: --font-size - 1.",
    )

    # help_group = parser.add_argument_group("Help")
    # help_group.add_argument(
    #     "-h",
    #     "--help",
    #     action="help",
    #     default=argparse.SUPPRESS,
    #     help="Show this help message and exit.",
    # )

    return parser

def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.target_list:
        file_check(args.target_list)
    if args.annotation:
        file_check(args.annotation)
    if not args.target and not args.target_list:
        raise ValueError("Either --target/-g or --target-list must be provided.")
    if args.target and args.target_list:
        raise ValueError("Please provide only one of --target/-g and --target-list.")

    if args.data_type == "rna" and args.frame != "all":
        print("Warning: --data-type rna ignores frame-specific coloring; --frame was reset to all.", flush=True)
        args.frame = "all"

    if args.line_width <= 0:
        raise ValueError("--line-width must be > 0.")
    if args.bar_width <= 0:
        raise ValueError("--bar-width must be > 0.")
    if args.y_max is not None and args.y_max <= 0:
        raise ValueError("--y-max must be > 0 when provided.")
    if args.clip_quantile <= 0 or args.clip_quantile > 1:
        raise ValueError("--clip-quantile must be in the interval (0, 1].")
    if args.clip_value is not None and args.clip_value <= 0:
        raise ValueError("--clip-value must be > 0 when provided.")
    if args.intron_scale <= 0:
        raise ValueError("--intron-scale must be > 0.")
    if args.x_margin < 0:
        raise ValueError("--x-margin must be >= 0.")
    if args.width <= 0:
        raise ValueError("--width must be > 0.")
    if args.per_sample_height <= 0:
        raise ValueError("--per-sample-height must be > 0.")
    if args.structure_height <= 0:
        raise ValueError("--structure-height must be > 0.")
    if args.max_height <= 0:
        raise ValueError("--max-height must be > 0.")
    if args.dpi <= 0:
        raise ValueError("--dpi must be > 0.")
    if args.font_size <= 0:
        raise ValueError("--font-size must be > 0.")
    for name in ["title_size", "label_size", "tick_size", "sample_label_size", "legend_size"]:
        value = getattr(args, name)
        if value is not None and value <= 0:
            raise ValueError("--" + name.replace("_", "-") + " must be > 0 when provided.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _safe_name(name: str) -> str:
    """Return a filesystem-safe target identifier."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_") or "target"


def _read_target_list(target_list: str) -> list[str]:
    """Read target IDs from a plain text target-list file."""
    targets: list[str] = []
    header_tokens = {"target", "target_id", "gene", "gene_id", "transcript", "transcript_id", "name", "id"}

    with open(target_list, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first_token = line.split()[0].strip()
            if line_number == 1 and first_token.lower() in header_tokens:
                continue
            targets.append(first_token)

    targets = list(dict.fromkeys(targets))
    if not targets:
        raise ValueError("No target IDs were found in --target-list: {file}".format(file=target_list))
    return targets


def _coordinate_jobs(args: Namespace) -> list[tuple[str, str | None, bool]]:
    """Return coordinate jobs as (coordinate, output_suffix, skip_if_missing_genome)."""
    if args.coordinate in {"auto", "both"}:
        return [
            ("transcript", "transcript", False),
            ("genome", "genome", True),
        ]
    return [(args.coordinate, None, False)]


def _run_one_geneplot(args: Namespace, output_suffix: str | None = None, skip_if_missing_genome: bool = False) -> bool:
    """Run one geneplot task and return whether output was written."""
    now_args = copy.copy(args)
    if output_suffix:
        now_args.output = args.output + "_" + output_suffix

    geneplot = Geneplot.Geneplot(now_args)
    geneplot.import_rpf()

    if skip_if_missing_genome and now_args.coordinate == "genome" and geneplot.resolved_coordinate != "genome":
        print(
            "Skip genome-coordinate output for {target}: no usable genome annotation was detected.".format(
                target=now_args.target
            ),
            flush=True,
        )
        return False

    geneplot.output_profile()
    geneplot.draw_geneplot()
    return True


def _run_target(args: Namespace) -> None:
    """Run all requested coordinate jobs for one target."""
    jobs = _coordinate_jobs(args)
    for coordinate, output_suffix, skip_if_missing_genome in jobs:
        now_args = copy.copy(args)
        now_args.coordinate = coordinate
        _run_one_geneplot(
            now_args,
            output_suffix=output_suffix,
            skip_if_missing_genome=skip_if_missing_genome,
        )


def _run_geneplot_pipeline(args: Namespace) -> None:
    """Run the gene-level plotting workflow."""
    if args.target_list:
        targets = _read_target_list(args.target_list)
        step_print(2, "Batch geneplot for {count:,} target(s).".format(count=len(targets)))
        for index, target in enumerate(targets, start=1):
            print("\nTarget {index}/{total}: {target}".format(index=index, total=len(targets), target=target), flush=True)
            now_args = copy.copy(args)
            now_args.target = target
            now_args.output = args.output + "_" + _safe_name(target)
            _run_target(now_args)
        return

    step_print(2, "Import target density profile and draw requested coordinate plot(s).")
    _run_target(args)


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Geneplot."""
    now_time()
    title_print('Draw gene-level RPF/RNA density plot.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_geneplot_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
