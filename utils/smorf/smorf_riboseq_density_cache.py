#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.007
# Function: Build reusable memory-mapped chromosome caches for Ribo-seq density.
# Input: Ordered WIG/bedGraph tracks and required chromosomes.
# Output: Sparse per-chromosome NumPy cache files and cache metadata.

"""High-performance reusable density cache for family evidence analysis."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from .smorf_riboseq_density import (
    ChromDensity,
    infer_density_format,
    iter_density_chromosomes,
    sort_bedgraph,
)


@dataclass(frozen=True, slots=True)
class CachedChromosome:
    """Store memory-mapped array paths for one chromosome."""

    starts_path: str
    ends_path: str
    values_path: str
    max_end: int
    byte_size: int


@dataclass(frozen=True, slots=True)
class CachedTrack:
    """Store cache metadata for one sample/strand density track."""

    sample: str
    strand: str
    source_path: str
    file_format: str
    chromosomes: Mapping[str, CachedChromosome]
    byte_size: int


@dataclass(frozen=True, slots=True)
class DensityCacheSet:
    """Store all cached tracks and aggregate cache statistics."""

    tracks: tuple[CachedTrack, ...]
    chromosomes: tuple[str, ...]
    byte_size: int


def _safe_token(text: str) -> str:
    """Return a filesystem-safe stable token."""
    digest = hashlib.sha1(text.encode("utf-8")).hexdigest()[:16]
    return digest


def _save_array(path: Path, values: np.ndarray) -> int:
    """Write one NumPy array without compression and return file size."""
    np.save(path, values, allow_pickle=False)
    return path.stat().st_size



def _open_binary(path: str | Path):
    """Open a plain or gzip-compressed density file in binary mode."""
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, "rb")
    return file_path.open("rb", buffering=8 * 1024 * 1024)


def _finalize_bedgraph_chromosome(
    chrom: str,
    starts: list[int],
    ends: list[int],
    values: list[float],
    ordered: bool,
) -> ChromDensity:
    """Build one sparse chromosome without sorting ordered input."""
    if not starts:
        return ChromDensity(
            chrom=chrom,
            starts=np.zeros(0, dtype=np.int64),
            ends=np.zeros(0, dtype=np.int64),
            values=np.zeros(0, dtype=np.float32),
            max_end=0,
        )
    starts_array = np.asarray(starts, dtype=np.int64)
    ends_array = np.asarray(ends, dtype=np.int64)
    values_array = np.abs(np.asarray(values, dtype=np.float32))
    if not ordered:
        order = np.lexsort((ends_array, starts_array))
        starts_array = starts_array[order]
        ends_array = ends_array[order]
        values_array = values_array[order]

    merged_starts: list[int] = []
    merged_ends: list[int] = []
    merged_values: list[float] = []
    previous_end = -1
    for start, end, value in zip(starts_array, ends_array, values_array):
        start_int = int(start)
        end_int = int(end)
        value_float = float(value)
        if start_int < 0 or end_int <= start_int:
            raise ValueError(
                f"Invalid density interval on {chrom}: "
                f"[{start_int}, {end_int})."
            )
        if start_int < previous_end:
            raise ValueError(
                f"Overlapping density intervals on {chrom}: "
                f"previous_end={previous_end}, start={start_int}."
            )
        if (
            merged_ends
            and start_int == merged_ends[-1]
            and math.isclose(
                merged_values[-1],
                value_float,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ):
            merged_ends[-1] = end_int
        else:
            merged_starts.append(start_int)
            merged_ends.append(end_int)
            merged_values.append(value_float)
        previous_end = end_int

    starts_output = np.asarray(merged_starts, dtype=np.int64)
    ends_output = np.asarray(merged_ends, dtype=np.int64)
    values_output = np.asarray(merged_values, dtype=np.float32)
    return ChromDensity(
        chrom=chrom,
        starts=starts_output,
        ends=ends_output,
        values=values_output,
        max_end=int(ends_output[-1]),
    )


def _iter_fast_bedgraph_chromosomes(
    path: str | Path,
):
    """Yield bedGraph chromosomes with a binary single-pass parser."""
    current_chrom: str | None = None
    seen: set[str] = set()
    starts: list[int] = []
    ends: list[int] = []
    values: list[float] = []
    ordered = True
    last_start = -1
    last_end = -1

    with _open_binary(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith((b"#", b"track", b"browser")):
                continue
            fields = stripped.split()
            if len(fields) < 4:
                raise ValueError(
                    f"bedGraph line {line_number} has fewer than four fields."
                )
            try:
                chrom = fields[0].decode("utf-8")
                start = int(fields[1])
                end = int(fields[2])
                value = abs(float(fields[3]))
            except (UnicodeDecodeError, ValueError) as error:
                raise ValueError(
                    f"Invalid bedGraph record at line {line_number}."
                ) from error
            if not math.isfinite(value):
                raise ValueError(
                    f"Invalid density value at line {line_number}: {value}."
                )

            if current_chrom is None:
                current_chrom = chrom
            elif chrom != current_chrom:
                if chrom in seen:
                    raise ValueError(
                        f"Density chromosome appears in multiple separated blocks: "
                        f"{chrom}."
                    )
                yield _finalize_bedgraph_chromosome(
                    current_chrom,
                    starts,
                    ends,
                    values,
                    ordered,
                )
                seen.add(current_chrom)
                current_chrom = chrom
                starts = []
                ends = []
                values = []
                ordered = True
                last_start = -1
                last_end = -1

            if start < last_start or (start == last_start and end < last_end):
                ordered = False
            starts.append(start)
            ends.append(end)
            values.append(value)
            last_start = start
            last_end = end

    if current_chrom is not None:
        yield _finalize_bedgraph_chromosome(
            current_chrom,
            starts,
            ends,
            values,
            ordered,
        )

def _cache_track_worker(
    task: tuple[int, str, str, str, str, str, tuple[str, ...]],
) -> tuple[int, CachedTrack]:
    """Parse one density track once and write chromosome caches."""
    (
        track_index,
        sample,
        strand,
        source_path,
        file_format,
        cache_root,
        required_chromosomes,
    ) = task
    required = set(required_chromosomes)
    track_directory = Path(cache_root) / f"track_{track_index:04d}"
    track_directory.mkdir(parents=True, exist_ok=True)
    resolved_format = infer_density_format(source_path, file_format)
    chromosome_cache: dict[str, CachedChromosome] = {}
    total_bytes = 0

    def populate(cache_source: str) -> None:
        nonlocal total_bytes
        density_iterator = (
            _iter_fast_bedgraph_chromosomes(cache_source)
            if resolved_format == "bedgraph"
            else iter_density_chromosomes(cache_source, resolved_format)
        )
        for density in density_iterator:
            if required and density.chrom not in required:
                continue
            token = _safe_token(density.chrom)
            starts_path = track_directory / f"{token}.starts.npy"
            ends_path = track_directory / f"{token}.ends.npy"
            values_path = track_directory / f"{token}.values.npy"
            byte_size = 0
            byte_size += _save_array(
                starts_path,
                np.asarray(density.starts, dtype=np.int64),
            )
            byte_size += _save_array(
                ends_path,
                np.asarray(density.ends, dtype=np.int64),
            )
            byte_size += _save_array(
                values_path,
                np.asarray(density.values, dtype=np.float32),
            )
            chromosome_cache[density.chrom] = CachedChromosome(
                starts_path=str(starts_path),
                ends_path=str(ends_path),
                values_path=str(values_path),
                max_end=int(density.max_end),
                byte_size=byte_size,
            )
            total_bytes += byte_size

    try:
        populate(source_path)
    except ValueError as error:
        can_sort = (
            resolved_format == "bedgraph"
            and "multiple separated blocks" in str(error)
        )
        if not can_sort:
            raise
        for cache_file in track_directory.glob("*.npy"):
            cache_file.unlink(missing_ok=True)
        chromosome_cache.clear()
        total_bytes = 0
        prepared = sort_bedgraph(source_path, track_directory)
        try:
            populate(prepared.path)
        finally:
            prepared.cleanup()

    manifest_path = track_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "sample": sample,
                "strand": strand,
                "source_path": source_path,
                "file_format": resolved_format,
                "chromosomes": sorted(chromosome_cache),
                "byte_size": total_bytes,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return (
        track_index,
        CachedTrack(
            sample=sample,
            strand=strand,
            source_path=source_path,
            file_format=resolved_format,
            chromosomes=chromosome_cache,
            byte_size=total_bytes,
        ),
    )


def _fork_context() -> mp.context.BaseContext | None:
    """Return a fork context on POSIX systems when available."""
    try:
        return mp.get_context("fork")
    except ValueError:
        return None


def build_density_cache(
    samples: Sequence[object],
    work_directory: str | Path,
    required_chromosomes: Sequence[str],
    threads: int,
    progress_callback: Callable[[int, int, CachedTrack], None] | None = None,
) -> DensityCacheSet:
    """Parse each density file once and cache sparse arrays by chromosome.

    Density files must keep each chromosome in one contiguous block. Coordinate
    order within a chromosome may be arbitrary because the existing density
    reader sorts and merges that chromosome in memory.
    """
    cache_root = Path(work_directory) / "density_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    tasks: list[tuple[int, str, str, str, str, str, tuple[str, ...]]] = []
    track_index = 0
    required = tuple(str(chrom) for chrom in required_chromosomes)
    for sample in samples:
        for track in sample.tracks:
            tasks.append(
                (
                    track_index,
                    str(sample.sample),
                    str(track.strand),
                    str(track.path),
                    str(track.file_format),
                    str(cache_root),
                    required,
                )
            )
            track_index += 1

    if not tasks:
        raise ValueError("No density tracks were provided.")

    requested = max(1, int(threads))
    worker_count = min(requested, len(tasks))
    results: list[CachedTrack | None] = [None] * len(tasks)
    completed = 0
    context = _fork_context()

    if worker_count == 1 or context is None:
        for task in tasks:
            index, cached = _cache_track_worker(task)
            results[index] = cached
            completed += 1
            if progress_callback is not None:
                progress_callback(completed, len(tasks), cached)
    else:
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=context,
        ) as executor:
            future_map = {
                executor.submit(_cache_track_worker, task): task[0]
                for task in tasks
            }
            for future in as_completed(future_map):
                index, cached = future.result()
                results[index] = cached
                completed += 1
                if progress_callback is not None:
                    progress_callback(completed, len(tasks), cached)

    tracks = tuple(track for track in results if track is not None)
    all_chromosomes = sorted(
        {
            chrom
            for track in tracks
            for chrom in track.chromosomes
        }
    )
    return DensityCacheSet(
        tracks=tracks,
        chromosomes=tuple(all_chromosomes),
        byte_size=sum(track.byte_size for track in tracks),
    )


def load_chromosome_map(
    cache: DensityCacheSet,
    chrom: str,
) -> dict[tuple[str, str], ChromDensity | None]:
    """Load one chromosome from every cached track using memory mapping."""
    output: dict[tuple[str, str], ChromDensity | None] = {}
    for track in cache.tracks:
        metadata = track.chromosomes.get(chrom)
        key = (track.sample, track.strand)
        if metadata is None:
            output[key] = None
            continue
        starts = np.load(metadata.starts_path, mmap_mode="r", allow_pickle=False)
        ends = np.load(metadata.ends_path, mmap_mode="r", allow_pickle=False)
        values = np.load(metadata.values_path, mmap_mode="r", allow_pickle=False)
        output[key] = ChromDensity(
            chrom=chrom,
            starts=starts,
            ends=ends,
            values=values,
            max_end=metadata.max_end,
        )
    return output


def chromosome_cache_bytes(cache: DensityCacheSet, chrom: str) -> int:
    """Return cached sparse-array bytes for one chromosome across tracks."""
    return sum(
        track.chromosomes[chrom].byte_size
        for track in cache.tracks
        if chrom in track.chromosomes
    )


def available_memory_bytes() -> int | None:
    """Return currently available physical memory when POSIX reports it."""
    try:
        pages = os.sysconf("SC_AVPHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
    except (AttributeError, OSError, ValueError):
        return None
    if pages <= 0 or page_size <= 0:
        return None
    return int(pages) * int(page_size)
