#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-09
# Version: 0.2.8-dev.004
# Function: Draw RPF metaplot profiles around start and stop codons.
# Input: RPF density file in JSONL or TXT format and optional transcript filter.
# Output: Metaplot summary table, outlier table, sample-level plots, and heatmap figures.

"""Command-line entry point for RPF metaplot analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.ribo import Metaplot

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
        description="Draw RPF metaplot profiles around start and stop codons from JSONL or TXT density files."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-r",
        dest="rpf",
        required=True,
        type=str,
        help="Input RPF density file in JSONL or TXT format.",
    )
    required_group.add_argument(
        "-o",
        dest="output",
        required=True,
        type=str,
        help="Output prefix.",
    )

    parser.add_argument(
        "-t",
        dest="transcript",
        required=False,
        type=str,
        default=None,
        help="Optional transcript filter table in TXT format. If provided, transcript_id is used when available.",
    )
    parser.add_argument(
        "-m",
        dest="min",
        required=False,
        type=int,
        default=50,
        help="Retain transcripts with at least this sample-specific CDS RPF count. Default: %(default)s.",
    )
    parser.add_argument(
        "--utr5",
        dest="utr5",
        type=int,
        required=False,
        default=20,
        help="Number of codons upstream of the start codon. Default: %(default)s AA.",
    )
    parser.add_argument(
        "--cds",
        dest="cds",
        type=int,
        required=False,
        default=50,
        help="Number of CDS codons used around start and stop codons. Default: %(default)s AA.",
    )
    parser.add_argument(
        "--utr3",
        dest="utr3",
        type=int,
        required=False,
        default=20,
        help="Number of codons downstream of the stop codon. Default: %(default)s AA.",
    )
    parser.add_argument(
        "-n",
        "--normal",
        dest="normal",
        action="store_true",
        required=False,
        default=False,
        help="Normalize RPF counts to RPM before metaplot aggregation. Default: %(default)s.",
    )
    parser.add_argument(
        "--mode",
        dest="mode",
        required=False,
        choices=["line", "bar", "both"],
        type=str,
        default="bar",
        help="Metaplot style for per-sample figures. Default: %(default)s.",
    )
    parser.add_argument(
        "--plot-transform",
        dest="plot_transform",
        required=False,
        choices=["none", "sqrt", "log", "log1p", "log2", "log10"],
        type=str,
        default="none",
        help=(
            "Transform plotted density values without changing output tables. "
            "'log' is an alias of log1p. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--scale",
        dest="scale",
        required=False,
        choices=["none", "row"],
        type=str,
        default="none",
        help=(
            "Heatmap scaling method. 'row' applies row-wise min-max scaling "
            "to emphasize each sample profile shape. Default: %(default)s."
        ),
    )

    outlier_group = parser.add_argument_group("Outlier arguments")
    outlier_group.add_argument(
        "--remove-outlier",
        dest="remove_outlier",
        action="store_true",
        required=False,
        default=False,
        help=(
            "Remove extreme gene-position RPF pileups before aggregation. "
            "This is useful for rRNA-like or mis-mapped fragments. Default: %(default)s."
        ),
    )
    outlier_group.add_argument(
        "--outlier-iqr",
        dest="outlier_iqr",
        type=float,
        required=False,
        default=8.0,
        help="IQR multiplier for extreme pileup detection. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-window",
        dest="outlier_window",
        type=int,
        required=False,
        default=5,
        help="Number of neighboring codons on each side used to estimate local background. Default: %(default)s.",
    )
    outlier_group.add_argument(
        "--outlier-local-fold",
        dest="outlier_local_fold",
        type=float,
        required=False,
        default=10.0,
        help="Minimum fold over local background for dynamic isolated-pileup filtering. Default: %(default)s.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.rpf)
    if args.transcript:
        file_check(args.transcript)

    if args.utr5 < 0:
        raise ValueError("--utr5 must be >= 0.")
    if args.cds <= 0:
        raise ValueError("--cds must be > 0.")
    if args.utr3 < 0:
        raise ValueError("--utr3 must be >= 0.")
    if args.min < 0:
        raise ValueError("-m/--min must be >= 0.")
    if args.outlier_iqr <= 0:
        raise ValueError("--outlier-iqr must be > 0.")
    if args.outlier_window < 1:
        raise ValueError("--outlier-window must be >= 1.")
    if args.outlier_local_fold <= 1:
        raise ValueError("--outlier-local-fold must be > 1.")


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _run_metaplot_pipeline(args: Namespace) -> None:
    """Run the metaplot analysis workflow."""
    metaplot = Metaplot.Metaplot(args)

    step_print(2, "Import the RPFs file.")
    metaplot.import_rpf()

    step_print(3, "Calculate metaplot profiles.")
    metaplot.calc_metaplot()

    step_print(4, "Output metaplot tables.")
    metaplot.output_meta()

    step_print(5, "Draw metaplot figures.")
    if args.mode in {"bar", "both"}:
        metaplot.draw_metaplot("bar")
    if args.mode in {"line", "both"}:
        metaplot.draw_metaplot("line")
    metaplot.draw_heatmap()

    step_print(6, "Write metaplot summary.")
    metaplot.write_summary()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for rpf_Metaplot."""
    now_time()
    title_print('Draw the metaplot.')
    step_print(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_metaplot_pipeline(args)

    complete_print()
    now_time()


if __name__ == "__main__":
    main()
