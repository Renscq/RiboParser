#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Compare significant peak overlap between mock-IP and flag-IP samples.
# Input: Mock-IP and flag-IP SeRP peak tables.
# Output: Annotated overlap tables and Venn diagrams.

"""Command-line entry point for SeRP peak-overlap analysis."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd

from utils.ribo.ArgsParser import args_print, file_check, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Compare significant peak overlap between mock-IP and flag-IP."
    )

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument(
        "-m",
        dest="mock",
        required=True,
        type=str,
        help="Input mock-IP peak table.",
    )
    required_group.add_argument(
        "-f",
        dest="flag",
        required=True,
        type=str,
        help="Input flag-IP peak table.",
    )

    output_group = parser.add_argument_group("Output arguments")
    output_group.add_argument(
        "--om",
        dest="out_mock",
        default=None,
        type=str,
        help="Mock-IP output prefix. Default: input filename without suffix.",
    )
    output_group.add_argument(
        "--of",
        dest="out_flag",
        default=None,
        type=str,
        help="Flag-IP output prefix. Default: input filename without suffix.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    file_check(args.mock, args.flag)


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)

    if args.out_mock is None:
        args.out_mock = str(Path(args.mock).with_suffix(""))
    if args.out_flag is None:
        args.out_flag = str(Path(args.flag).with_suffix(""))

    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print a standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _import_peak_region(
    mock_file: str,
    flag_file: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Import peak tables and retain significant peaks."""
    mock_peak = pd.read_csv(mock_file, sep="\t")
    flag_peak = pd.read_csv(flag_file, sep="\t")

    required_columns = {"BHFDR", "transcripts", "peak_start", "peak_end"}
    for label, table in (("mock-IP", mock_peak), ("flag-IP", flag_peak)):
        missing = required_columns.difference(table.columns)
        if missing:
            raise ValueError(
                f"{label} peak table is missing columns: {', '.join(sorted(missing))}"
            )

    mock_fdr = pd.to_numeric(mock_peak["BHFDR"].replace("-", pd.NA), errors="coerce")
    flag_fdr = pd.to_numeric(flag_peak["BHFDR"].replace("-", pd.NA), errors="coerce")

    significant_mock = mock_peak.loc[mock_fdr < 0.05].copy()
    significant_flag = flag_peak.loc[flag_fdr < 0.05].copy()
    return significant_mock, significant_flag


def _annotate_overlap(
    left_peak: pd.DataFrame,
    right_peak: pd.DataFrame,
    marker: str,
) -> pd.DataFrame:
    """Annotate interval overlap against the second peak table."""
    annotated = left_peak.copy()
    annotated["overlap"] = marker

    right_by_transcript = {
        transcript: group[["peak_start", "peak_end"]].to_numpy()
        for transcript, group in right_peak.groupby("transcripts", sort=False)
    }

    for row_index, row in annotated.iterrows():
        candidate_intervals = right_by_transcript.get(row["transcripts"])
        if candidate_intervals is None:
            continue

        left_start = int(row["peak_start"])
        left_end = int(row["peak_end"])
        for right_start, right_end in candidate_intervals:
            if int(right_start) <= left_end and left_start <= int(right_end):
                annotated.at[row_index, "overlap"] = "overlap"
                break

    return annotated


def _merge_overlap(
    significant_mock: pd.DataFrame,
    significant_flag: pd.DataFrame,
    out_mock: str,
    out_flag: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Annotate and write the two peak-overlap tables."""
    mock_annotated = _annotate_overlap(
        significant_mock,
        significant_flag,
        "mock-IP",
    )
    flag_annotated = _annotate_overlap(
        significant_flag,
        significant_mock,
        "flag-IP",
    )

    mock_annotated.to_csv(
        f"{out_mock}.overlap.txt",
        sep="\t",
        header=True,
        index=False,
    )
    flag_annotated.to_csv(
        f"{out_flag}.overlap.txt",
        sep="\t",
        header=True,
        index=False,
    )
    return mock_annotated, flag_annotated


def _plot_venn(
    mock_annotated: pd.DataFrame,
    flag_annotated: pd.DataFrame,
    out_mock: str,
    out_flag: str,
) -> None:
    """Draw gene-level and peak-level overlap diagrams."""
    mock_genes = set(mock_annotated["transcripts"])
    flag_genes = set(flag_annotated["transcripts"])

    mock_only = int((mock_annotated["overlap"] == "mock-IP").sum())
    flag_only = int((flag_annotated["overlap"] == "flag-IP").sum())
    shared_peak_count = int(
        round(
            (
                (mock_annotated["overlap"] == "overlap").sum()
                + (flag_annotated["overlap"] == "overlap").sum()
            )
            / 2
        )
    )

    figure, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=300)

    venn2(
        [mock_genes, flag_genes],
        set_labels=["mock-IP", "flag-IP"],
        ax=axes[0],
    )
    axes[0].set_title("Gene overlap")

    venn2(
        subsets=(mock_only, flag_only, shared_peak_count),
        set_labels=["mock-IP", "flag-IP"],
        ax=axes[1],
    )
    axes[1].set_title("Peak-region overlap")

    figure.tight_layout()
    output_prefix = f"{out_mock}_vs_{Path(out_flag).name}_venn"
    figure.savefig(f"{output_prefix}.pdf")
    figure.savefig(f"{output_prefix}.png", dpi=300)
    plt.close(figure)


def _run_overlap_pipeline(args: Namespace) -> None:
    """Run the SeRP peak-overlap workflow."""
    _print_step(2, "Import significant peak regions.")
    significant_mock, significant_flag = _import_peak_region(
        args.mock,
        args.flag,
    )

    _print_step(3, "Annotate overlapping peak regions.")
    mock_annotated, flag_annotated = _merge_overlap(
        significant_mock,
        significant_flag,
        args.out_mock,
        args.out_flag,
    )

    _print_step(4, "Draw the overlap Venn diagrams.")
    _plot_venn(
        mock_annotated,
        flag_annotated,
        args.out_mock,
        args.out_flag,
    )


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for serp_overlap."""
    now_time()
    print("\nCompare mock-IP and flag-IP peak overlap.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_overlap_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
