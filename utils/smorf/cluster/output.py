#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Write deterministic cluster shards.
# Input: Filtered candidates and family records.
# Output: Family, member, removed, and summary shards.

"""Write deterministic cluster shards."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from types import TracebackType
from typing import BinaryIO

from .config import _BUFFER_BYTES, _FLUSH_BYTES


class _ShardWriter:
    """Write one headerless worker shard with large byte buffers."""

    def __init__(self, directory: Path, shard_id: int) -> None:
        self.paths = {
            "family": directory / f"family.{shard_id:04d}.part",
            "members": directory / f"members.{shard_id:04d}.part",
            "removed": directory / f"removed.{shard_id:04d}.part",
        }
        self.handles: dict[str, BinaryIO] = {}
        self.buffers = {
            "family": bytearray(),
            "members": bytearray(),
            "removed": bytearray(),
        }

    def __enter__(self) -> _ShardWriter:
        self.handles = {
            name: path.open("wb", buffering=_BUFFER_BYTES) for name, path in self.paths.items()
        }
        return self

    def _append(self, name: str, line: bytes) -> None:
        buffer = self.buffers[name]
        buffer.extend(line)
        if len(buffer) >= _FLUSH_BYTES:
            self.handles[name].write(buffer)
            buffer.clear()

    def write_family(self, fields: Sequence[bytes]) -> None:
        self._append("family", b"\t".join(fields) + b"\n")

    def write_member(self, fields: Sequence[bytes]) -> None:
        self._append("members", b"\t".join(fields) + b"\n")

    def write_removed(self, fields: Sequence[bytes]) -> None:
        self._append("removed", b"\t".join(fields) + b"\n")

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        for name, buffer in self.buffers.items():
            if buffer:
                self.handles[name].write(buffer)
                buffer.clear()
        for handle in self.handles.values():
            handle.close()
        if exception_type is not None:
            for path in self.paths.values():
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass
        return False
