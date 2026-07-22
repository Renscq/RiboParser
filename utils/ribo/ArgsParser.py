#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Provide standardized console-output and argument-validation helpers.
# Input: argparse namespaces, workflow messages, result values, and file paths.
# Output: Consistently formatted RiboParser console messages.

"""Shared console-output helpers for RiboParser command-line programs.

The module keeps command output intentionally plain. It does not add module
names, log levels in brackets, or per-class prefixes. All command-line programs
can use the same title, step, argument, progress, warning, and result formats.
"""

from __future__ import annotations

import logging
import os
import sys
import time
from argparse import Namespace
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any, TextIO


logger = logging.getLogger("RiboParser")
MIN_LABEL_WIDTH = 12


def now_time() -> None:
    """Print the current local time."""
    print(
        time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        flush=True,
    )


def title_print(message: str) -> None:
    """Print one workflow title.

    Args:
        message: Title text.
    """
    print(f"\n{message}", flush=True)


def step_print(step: int, message: str) -> None:
    """Print one numbered workflow step.

    Args:
        step: One-based workflow step number.
        message: Step description.
    """
    print(f"\nStep{int(step)}: {message}", flush=True)


def message_print(
    message: str,
    *,
    stream: TextIO | None = None,
) -> None:
    """Print one plain informational message.

    Args:
        message: Message text.
        stream: Optional output stream. Defaults to standard output.
    """
    print(
        message,
        file=sys.stdout if stream is None else stream,
        flush=True,
    )


def progress_print(message: str) -> None:
    """Print one standardized progress message.

    Args:
        message: Progress text without a prefix.
    """
    print(f"Progress: {message}", flush=True)


def warning_print(message: str) -> None:
    """Print one standardized warning message.

    Args:
        message: Warning text without a prefix.
    """
    print(f"Warning: {message}", file=sys.stderr, flush=True)


def _display_value(value: Any) -> str:
    """Convert one scalar value to stable display text."""
    if value is None:
        return "None"
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _is_sequence_value(value: Any) -> bool:
    """Return whether a value should use multiline sequence formatting."""
    return isinstance(value, Sequence) and not isinstance(
        value,
        (str, bytes, bytearray),
    )


def key_value_print(
    items: Mapping[str, Any] | Iterable[tuple[str, Any]],
    *,
    min_width: int = MIN_LABEL_WIDTH,
) -> None:
    """Print aligned key-value rows.

    Sequence values are printed as an item count followed by indented rows.
    This prevents long input-file lists from appearing as a Python list on one
    unreadable line.

    Args:
        items: Mapping or iterable of ``(label, value)`` pairs.
        min_width: Minimum label width.
    """
    pairs = list(items.items() if isinstance(items, Mapping) else items)
    if not pairs:
        return

    labels = [str(label) for label, _ in pairs]
    width = max(int(min_width), max(len(label) for label in labels))
    continuation = " " * (width + 3)

    for label, value in pairs:
        label_text = str(label)
        if _is_sequence_value(value):
            values = list(value)
            print(
                f"{label_text:<{width}} : {len(values)} item(s)",
                flush=True,
            )
            for index, item in enumerate(values, start=1):
                print(
                    f"{continuation}{index:>2}. {_display_value(item)}",
                    flush=True,
                )
            continue

        print(
            f"{label_text:<{width}} : {_display_value(value)}",
            flush=True,
        )

    sys.stdout.flush()


def args_print(args: Namespace) -> None:
    """Print parsed command-line arguments with dynamic alignment.

    Args:
        args: Parsed argparse namespace.
    """
    key_value_print(vars(args))


def result_print(
    items: Mapping[str, Any] | Iterable[tuple[str, Any]],
) -> None:
    """Print an aligned result summary.

    Args:
        items: Mapping or iterable of result-label/value pairs.
    """
    key_value_print(items)


def file_check(*files: str | os.PathLike[str]) -> None:
    """Check that all supplied paths exist.

    Args:
        *files: Input file paths.

    Raises:
        SystemExit: If any path does not exist.
    """
    missing = [str(path) for path in files if not os.path.exists(path)]
    if not missing:
        return

    for path in missing:
        print(f"File not found: {path}", file=sys.stderr, flush=True)
    raise SystemExit(1)
