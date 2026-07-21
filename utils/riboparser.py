#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Display RiboParser version, citation, dependency, and module information.
# Input: Command-line options selecting the requested package information.
# Output: RiboParser package information printed to standard output.

"""Command-line entry point for the RiboParser package information utility."""

from __future__ import annotations

import argparse
from argparse import Namespace
from collections.abc import Sequence

from utils.data import RiboParser
from utils.ribo.ArgsParser import args_print, now_time


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Display RiboParser package information."
    )

    information_group = parser.add_argument_group("Information arguments")
    information_group.add_argument(
        "-v",
        "--version",
        dest="version",
        action="store_true",
        default=False,
        help="Display the installed RiboParser version.",
    )
    information_group.add_argument(
        "-c",
        "--citation",
        dest="citation",
        action="store_true",
        default=False,
        help="Display the recommended RiboParser citation.",
    )
    information_group.add_argument(
        "-d",
        "--dependency",
        dest="dependency",
        action="store_true",
        default=False,
        help="Check external software dependencies.",
    )
    information_group.add_argument(
        "-m",
        "--module",
        dest="module",
        action="store_true",
        default=False,
        help="Check the installed RiboParser Python modules.",
    )

    return parser


def _validate_args(args: Namespace) -> None:
    """Validate command-line arguments."""
    if not any((args.version, args.citation, args.dependency, args.module)):
        args.version = True


def _parse_args(argv: Sequence[str] | None = None) -> Namespace:
    """Parse, validate, and print command-line arguments."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    args_print(args)
    return args


def _print_step(step: int, message: str) -> None:
    """Print a standardized pipeline step message."""
    print(f"\nStep{step}: {message}", flush=True)


def _run_information_pipeline(args: Namespace) -> None:
    """Display the requested RiboParser package information."""
    if args.version:
        _print_step(2, "Display RiboParser version information.")
        RiboParser.RiboParserInfo.show_version()

    if args.citation:
        _print_step(3, "Display RiboParser citation information.")
        RiboParser.RiboParserInfo.show_citation()

    if args.dependency:
        _print_step(4, "Check external software dependencies.")
        RiboParser.RiboParserInfo.check_dependencies()

    if args.module:
        _print_step(5, "Check installed RiboParser modules.")
        RiboParser.RiboParserInfo.check_package_modules()


def main(argv: Sequence[str] | None = None) -> None:
    """Command-line entry point for riboparser."""
    now_time()
    print("\nDisplay RiboParser package information.", flush=True)
    _print_step(1, "Checking the input arguments.")

    args = _parse_args(argv)
    _run_information_pipeline(args)

    print("\nAll done.", flush=True)
    now_time()


if __name__ == "__main__":
    main()
