#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-31
# Version: 0.2.8.24-dev.016
# Function: Evaluate family evidence while separating calibration controls.
# Input: Family tables, scanner ORF table/genePred, and density design table.
# Output: Family evidence, reliable smORFs/genePred, summary, and log.

"""Bottom-up family-aware smORF evidence engine.

The engine deliberately avoids the previous representative-by-unit workflow.
It builds a SQLite family index, caches each density track once, evaluates one
family scaffold per sample, derives project reliability from independent sample support across all
sample groups, and resolves alternative starts from replicated extension evidence.

Only standard-library modules, NumPy, and the existing RiboParser density reader
are required. A malformed family is recorded as uncertain instead of aborting
an entire genome-wide run.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import shutil
import sqlite3
import tempfile
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Final, Iterable, Iterator, Mapping, Sequence, TextIO

import numpy as np

from utils.ribo.ArgsParser import message_print, progress_print, warning_print
from utils.smorf.smorf_riboseq_density import (
    ChromDensity,
    infer_density_format,
    iter_density_chromosomes,
    prepare_density,
)

ENGINE_VERSION: Final[str] = "0.2.8.24-dev.016"
SCHEMA_VERSION: Final[int] = 4
STOP_CODONS: Final[frozenset[str]] = frozenset({"TAA", "TAG", "TGA"})
ANNOTATED_CATEGORIES: Final[frozenset[str]] = frozenset(
    {"annotated_orf", "annotated_morf"}
)
FAMILY_BATCH_SIZE: Final[int] = 10_000
SQL_BATCH_SIZE: Final[int] = 50_000
MAX_WORKERS: Final[int] = 8
MEMORY_PER_WORKER_BYTES: Final[int] = 4 * 1024**3
LEVEL_RANK: Final[dict[str, int]] = {
    "NoEvidence": 0,
    "LowConfidence": 1,
    "MediumConfidence": 2,
    "HighConfidence": 3,
}

MANUAL_THRESHOLD_DEFAULTS: Final[dict[str, float | int]] = {
    "min_rpf_sum": 5.0,
    "min_rpf_per_codon": 0.10,
    "min_covered_codon": 3,
    "min_coverage_ratio": 0.10,
    "moderate_periodicity": 0.50,
    "strong_periodicity": 0.60,
    "min_window_rpf": 3.0,
    "min_window_covered": 3,
}
ADVANCED_DEFAULTS: Final[dict[str, float | int]] = {
    "short_max_codons": 30,
    "long_min_codons": 100,
    "window_codons": 20,
    "window_step_codons": 5,
    "min_supported_windows": 2,
    "min_window_gap_codons": 15,
    "min_signal_span": 0.35,
    "localized_span_max": 0.20,
    "localized_top_window_fraction": 0.70,
    "boundary_codons": 5,
}
TARGET_LENGTH_CATEGORIES: Final[frozenset[str]] = frozenset(
    {"uorf", "dorf", "lncorf"}
)


@dataclass(frozen=True, slots=True)
class DensityTrack:
    """Describe one P-site density track."""

    sample: str
    group: str
    strand: str
    path: str
    file_format: str


@dataclass(frozen=True, slots=True)
class Thresholds:
    """Store sample-specific evidence thresholds."""

    sample: str
    source: str
    control_count: int
    min_rpf_sum: float
    min_rpf_per_codon: float
    min_covered_codon: int
    min_coverage_ratio: float
    moderate_periodicity: float
    strong_periodicity: float
    min_window_rpf: float
    min_window_covered: int


@dataclass(frozen=True, slots=True)
class EngineConfig:
    """Store biological, calibration, grouping, and runtime configuration."""

    evidence_mode: str = "manual"
    group_column: str = "group"
    reliable_sample: int = 2

    min_rpf_sum: float = 5.0
    min_rpf_per_codon: float = 0.10
    min_covered_codon: int = 3
    min_coverage_ratio: float = 0.10
    moderate_periodicity: float = 0.50
    strong_periodicity: float = 0.60
    min_window_rpf: float = 3.0
    min_window_covered: int = 3

    short_max_codons: int = 30
    long_min_codons: int = 100
    window_codons: int = 20
    window_step_codons: int = 5
    min_supported_windows: int = 2
    min_window_gap_codons: int = 15
    min_signal_span: float = 0.35
    localized_span_max: float = 0.20
    localized_top_window_fraction: float = 0.70
    boundary_codons: int = 5
    hidden_overrides: frozenset[str] = frozenset()

    positive_quantile: float = 0.20
    positive_min_controls: int = 30
    positive_max_controls: int = 5000
    threads: int = 1



@dataclass(frozen=True, slots=True)
class RepresentativeGeometry:
    """Store one family representative in scaffold coordinates."""

    orf_id: str
    transcript_id: str
    category: str
    start_codon: str
    nt_length: int
    coding_nt_length: int
    aa_length: int
    starts: tuple[int, ...]
    ends: tuple[int, ...]
    start_offset: int


@dataclass(frozen=True, slots=True)
class FamilyGeometry:
    """Store one validated family scaffold."""

    family_id: str
    gene_id: str
    chrom: str
    strand: str
    category: str
    family_type: str
    family_size: int
    structural_primary: str
    scaffold_orf_id: str
    scaffold_nt_length: int
    scaffold_coding_nt_length: int
    scaffold_starts: tuple[int, ...]
    scaffold_ends: tuple[int, ...]
    common_body_orf: str
    common_body_start: int
    common_body_end: int
    representatives: tuple[RepresentativeGeometry, ...]


@dataclass(frozen=True, slots=True)
class InvalidFamily:
    """Store a family that cannot be represented by one scaffold."""

    family_id: str
    gene_id: str
    chrom: str
    strand: str
    category: str
    family_type: str
    family_size: int
    structural_primary: str
    reason: str


@dataclass(frozen=True, slots=True)
class SegmentFeatures:
    """Store sufficient statistics for one transcript segment."""

    rpf_sum: float
    rpf_per_codon: float
    covered_codon: int
    coverage_ratio: float
    frame0_density: float
    frame1_density: float
    frame2_density: float
    frame0_ratio: float
    supported_windows: int
    distributed_windows: int
    signal_span: float
    top_window_fraction: float
    start_rpf: float
    body_rpf: float
    end_rpf: float
    localized_only: bool
    score: float


@dataclass(frozen=True, slots=True)
class EngineResult:
    """Store formal output paths and project counts."""

    master_output: str
    reliable_output: str
    reliable_genepred_output: str
    summary_output: str
    total_families: int
    reliable_families: int
    uncertain_families: int
    no_evidence_families: int
    reliable_smorfs: int
    sample_count: int
    effective_workers: int


@dataclass(frozen=True, slots=True)
class ChromosomeTask:
    """Store a resumable chromosome evidence task."""

    chromosome: str
    database_path: str
    work_directory: str
    output_prefix: str
    tracks: tuple[DensityTrack, ...]
    thresholds: tuple[Thresholds, ...]
    config: EngineConfig
    group_names: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ChromosomeResult:
    """Store one chromosome shard and count summary."""

    chromosome: str
    part_path: str
    count_path: str
    total_families: int
    reliable_families: int
    uncertain_families: int
    no_evidence_families: int
    reliable_smorfs: int
    invalid_families: int


@dataclass(slots=True)
class _FamilyAccumulator:
    """Accumulate sample, replicate-group, and start-extension evidence."""

    any_signal: bool
    any_evidence: bool
    support_count: int
    high_count: int
    group_mask: int
    group_support_counts: dict[int, int]
    group_high_counts: dict[int, int]
    supporting_samples: list[str]
    best_level: str
    best_score: float
    best_sample: str
    best_features: SegmentFeatures | None
    segment_support_counts: list[int]
    segment_group_support_counts: list[dict[int, int]]
    segment_best_scores: list[float]
    segment_support_samples: list[list[str]]



class EvidenceEngineError(RuntimeError):
    """Represent a stage-specific evidence-engine failure."""


class _StageLogger:
    """Write a persistent stage log and mirror messages to stdout."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, message: str) -> None:
        text = str(message)
        message_print(text)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(text + "\n")

    def error(self, stage: str, error: BaseException) -> None:
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(f"ERROR stage={stage}: {error}\n")
            handle.write(traceback.format_exc())
            handle.write("\n")


def _smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain or gzip-compressed text."""
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8", newline="")
    return file_path.open(mode, encoding="utf-8", newline="")


def _split_ints(value: Any) -> tuple[int, ...]:
    """Parse comma-separated integer coordinates."""
    text = str(value).strip().rstrip(",")
    if not text:
        return ()
    return tuple(int(item) for item in text.split(",") if item != "")


def _coding_nt_length(record: Mapping[str, Any]) -> int:
    """Return coding length with a complete terminal stop removed."""
    length = int(record["nt_length"])
    stop = str(record.get("stop_codon", "")).upper().replace("U", "T")
    completeness = str(record.get("completeness", "complete")).lower()
    complete = completeness in {"", "complete", "cmpl", "full", "true", "yes"}
    if complete and stop in STOP_CODONS and length >= 6:
        length -= 3
    return max(0, length - length % 3)


def _file_signature(path: str | Path) -> dict[str, Any]:
    """Return a stable file signature."""
    resolved = Path(path).expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _run_signature(args: object, tracks: Sequence[DensityTrack]) -> str:
    """Hash all material inputs and analysis options."""
    files = {
        "family_table": _file_signature(args.family_table),
        "family_members": _file_signature(args.family_members),
        "orf_source": _file_signature(args.orf_source),
        "orf_genepred": _file_signature(args.orf_genepred),
        "tracks": [_file_signature(track.path) for track in tracks],
    }
    options = {
        key: value
        for key, value in sorted(vars(args).items())
        if key != "threads"
        and isinstance(value, (str, int, float, bool, type(None)))
    }
    payload = json.dumps(
        {
            "engine": ENGINE_VERSION,
            "schema": SCHEMA_VERSION,
            "files": files,
            "options": options,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _safe_name(value: str) -> str:
    """Return a deterministic filesystem-safe label."""
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:10]
    cleaned = "".join(character if character.isalnum() else "_" for character in value)
    cleaned = cleaned.strip("_")[:80] or "item"
    return f"{cleaned}.{digest}"


def _read_memory_limit() -> int | None:
    """Return the effective cgroup or host memory limit."""
    values: list[int] = []
    for path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            text = path.read_text(encoding="utf-8").strip()
        except OSError:
            continue
        if text and text != "max":
            try:
                value = int(text)
            except ValueError:
                continue
            if 0 < value < 1 << 60:
                values.append(value)
    try:
        with Path("/proc/meminfo").open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    values.append(int(line.split()[1]) * 1024)
                    break
    except OSError:
        pass
    return min(values) if values else None


def _effective_workers(requested: int, tasks: int) -> int:
    """Cap chromosome workers by memory and task count."""
    workers = max(1, min(int(requested), int(tasks), MAX_WORKERS))
    memory = _read_memory_limit()
    if memory is not None:
        workers = min(
            workers,
            max(1, int(memory * 0.70 // MEMORY_PER_WORKER_BYTES)),
        )
    return max(1, workers)


def _parse_density_list(
    path: str | Path,
    group_column: str,
) -> list[DensityTrack]:
    """Parse a precise header-based density design table.

    The design must contain sample, strand, one path column, and the biological
    replicate group column selected by ``--group``. Density format is inferred
    from each data-file path.
    """
    design_path = Path(path).expanduser().resolve()
    with _smart_open(design_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original_header = tuple(reader.fieldnames or ())
        if not original_header:
            raise ValueError(
                "Density list requires a tab-delimited header."
            )
        header_lookup = {
            str(name).strip().lower(): str(name)
            for name in original_header
        }
        path_name = next(
            (
                header_lookup[name]
                for name in ("path", "file", "density")
                if name in header_lookup
            ),
            None,
        )
        sample_name = next(
            (
                header_lookup[name]
                for name in ("sample", "name")
                if name in header_lookup
            ),
            None,
        )
        if sample_name is None or "strand" not in header_lookup:
            raise ValueError(
                "Density list requires sample/name, strand, and "
                "path/file/density columns."
            )
        group_key = str(group_column).strip().lower()
        if group_key not in {"sample", "name"} and group_key not in header_lookup:
            raise ValueError(
                f"Density-list group column was not found: {group_column}"
            )
        if path_name is None:
            raise ValueError(
                "Density list requires one path column: path, file, or density."
            )

        strand_name = header_lookup["strand"]
        group_name = (
            sample_name
            if group_key in {"sample", "name"}
            else header_lookup[group_key]
        )
        tracks: list[DensityTrack] = []
        for line_number, row in enumerate(reader, start=2):
            sample = str(row.get(sample_name, "")).strip()
            group = str(row.get(group_name, "")).strip()
            strand_text = str(row.get(strand_name, "")).strip().lower()
            density_text = str(row.get(path_name, "")).strip()
            if not sample or not group or not strand_text or not density_text:
                raise ValueError(
                    f"Density-list row {line_number} contains an empty required value."
                )
            strand_alias = {
                "+": "+",
                "plus": "+",
                "forward": "+",
                "-": "-",
                "minus": "-",
                "reverse": "-",
                ".": ".",
                "unstranded": ".",
            }
            if strand_text not in strand_alias:
                raise ValueError(
                    f"Density-list row {line_number} has invalid strand: {strand_text}"
                )
            density_path = Path(density_text).expanduser()
            if not density_path.is_absolute():
                density_path = design_path.parent / density_path
            density_path = density_path.resolve()
            if not density_path.is_file():
                raise FileNotFoundError(density_path)
            tracks.append(
                DensityTrack(
                    sample=sample,
                    group=group,
                    strand=strand_alias[strand_text],
                    path=str(density_path),
                    file_format=infer_density_format(density_path),
                )
            )
    if not tracks:
        raise ValueError("Density list contains no data rows.")
    return tracks



def read_density_tracks(args: object) -> list[DensityTrack]:
    """Read the required density design and validate sample/strand uniqueness."""
    density_list = getattr(args, "density_list", None)
    if not density_list:
        raise ValueError("--density-list is required.")
    tracks = _parse_density_list(
        density_list,
        str(getattr(args, "group_column", "group")),
    )

    seen: set[tuple[str, str]] = set()
    sample_groups: dict[str, str] = {}
    for track in tracks:
        key = (track.sample, track.strand)
        if key in seen:
            raise ValueError(
                "Only one density track per sample/strand is supported: "
                f"{track.sample} {track.strand}"
            )
        seen.add(key)
        previous_group = sample_groups.setdefault(track.sample, track.group)
        if previous_group != track.group:
            raise ValueError(
                f"Sample {track.sample} is assigned to multiple groups: "
                f"{previous_group}, {track.group}"
            )
    return tracks



def _connect_database(path: Path) -> sqlite3.Connection:
    """Open the disk index with memory-bounded pragmas."""
    connection = sqlite3.connect(path, timeout=120.0)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-262144")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    return connection


def _create_database_schema(connection: sqlite3.Connection) -> None:
    """Create the family and representative index schema."""
    connection.executescript(
        """
        CREATE TABLE families (
            family_id TEXT PRIMARY KEY,
            gene_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            strand TEXT NOT NULL,
            category TEXT NOT NULL,
            family_type TEXT NOT NULL,
            family_size INTEGER NOT NULL,
            structural_primary TEXT NOT NULL
        );
        CREATE TABLE representatives (
            family_id TEXT NOT NULL,
            orf_id TEXT PRIMARY KEY,
            is_primary INTEGER NOT NULL DEFAULT 0
        );
        CREATE INDEX representatives_family_idx
            ON representatives(family_id);
        CREATE TABLE geometry (
            orf_id TEXT PRIMARY KEY,
            transcript_id TEXT NOT NULL,
            category TEXT NOT NULL,
            start_codon TEXT NOT NULL,
            stop_codon TEXT NOT NULL,
            completeness TEXT NOT NULL,
            chrom TEXT NOT NULL,
            strand TEXT NOT NULL,
            nt_length INTEGER NOT NULL,
            aa_length INTEGER NOT NULL,
            exon_starts TEXT NOT NULL,
            exon_ends TEXT NOT NULL
        );
        CREATE INDEX geometry_chrom_idx ON geometry(chrom);
        """
    )


def _required_columns(header: Sequence[str], required: Sequence[str], name: str) -> None:
    """Validate a tabular header."""
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"{name} is missing column(s): {', '.join(missing)}")


def build_family_database(
    family_table: str | Path,
    family_members: str | Path,
    orf_source: str | Path,
    database_path: Path,
    signature: str,
    logger: _StageLogger,
) -> None:
    """Build or reuse a disk-backed family/representative index."""
    marker = database_path.with_suffix(".complete.json")
    if database_path.is_file() and marker.is_file():
        state = json.loads(marker.read_text(encoding="utf-8"))
        if state.get("signature") == signature:
            logger.write("Reuse family SQLite index.")
            return

    database_path.unlink(missing_ok=True)
    marker.unlink(missing_ok=True)
    connection = _connect_database(database_path)
    try:
        _create_database_schema(connection)

        logger.write("Index family primary records.")
        with _smart_open(family_table, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            _required_columns(
                header,
                (
                    "family_id", "gene_id", "chrom", "strand", "category",
                    "family_type", "family_size", "orf_id", "transcript_id",
                    "start_codon", "stop_codon", "completeness", "nt_length",
                    "aa_length", "exon_starts", "exon_ends",
                ),
                "family table",
            )
            family_rows: list[tuple[Any, ...]] = []
            rep_rows: list[tuple[Any, ...]] = []
            geometry_rows: list[tuple[Any, ...]] = []
            for number, row in enumerate(reader, start=1):
                family_rows.append(
                    (
                        row["family_id"], row["gene_id"], row["chrom"],
                        row["strand"], row["category"], row["family_type"],
                        int(row["family_size"]), row["orf_id"],
                    )
                )
                rep_rows.append((row["family_id"], row["orf_id"], 1))
                geometry_rows.append(
                    (
                        row["orf_id"], row["transcript_id"], row["category"],
                        row["start_codon"], row["stop_codon"],
                        row["completeness"], row["chrom"], row["strand"],
                        int(row["nt_length"]), int(row["aa_length"]),
                        row["exon_starts"], row["exon_ends"],
                    )
                )
                if len(family_rows) >= SQL_BATCH_SIZE:
                    connection.executemany(
                        "INSERT INTO families VALUES (?,?,?,?,?,?,?,?)",
                        family_rows,
                    )
                    connection.executemany(
                        "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                        rep_rows,
                    )
                    connection.executemany(
                        "INSERT OR REPLACE INTO geometry VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                        geometry_rows,
                    )
                    connection.commit()
                    family_rows.clear(); rep_rows.clear(); geometry_rows.clear()
                    progress_print(f"family index: {number:,}")
            if family_rows:
                connection.executemany(
                    "INSERT INTO families VALUES (?,?,?,?,?,?,?,?)",
                    family_rows,
                )
                connection.executemany(
                    "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                    rep_rows,
                )
                connection.executemany(
                    "INSERT OR REPLACE INTO geometry VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    geometry_rows,
                )
                connection.commit()

        logger.write("Index alternative-start representative relationships.")
        with _smart_open(family_members, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            _required_columns(
                header,
                ("family_id", "representative_orf_id"),
                "family members",
            )
            rows: list[tuple[str, str, int]] = []
            for number, row in enumerate(reader, start=1):
                representative = row["representative_orf_id"]
                if not representative:
                    continue
                rows.append(
                    (
                        row["family_id"],
                        representative,
                        int(representative == row.get("primary_orf_id", "")),
                    )
                )
                if len(rows) >= SQL_BATCH_SIZE:
                    connection.executemany(
                        "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                        rows,
                    )
                    connection.commit(); rows.clear()
                    progress_print(f"member index: {number:,}")
            if rows:
                connection.executemany(
                    "INSERT OR IGNORE INTO representatives VALUES (?,?,?)",
                    rows,
                )
                connection.commit()

        logger.write("Read source ORF table and fill representative geometry.")
        with _smart_open(orf_source, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            _required_columns(
                header,
                (
                    "orf_id", "transcript_id", "category", "start_codon",
                    "stop_codon", "completeness", "chrom", "strand",
                    "nt_length", "aa_length", "exon_starts", "exon_ends",
                ),
                "source ORF table",
            )
            connection.execute(
                """
                CREATE TEMP TABLE source_batch (
                    orf_id TEXT, transcript_id TEXT, category TEXT,
                    start_codon TEXT, stop_codon TEXT, completeness TEXT,
                    chrom TEXT, strand TEXT, nt_length INTEGER,
                    aa_length INTEGER, exon_starts TEXT, exon_ends TEXT
                )
                """
            )
            batch: list[tuple[Any, ...]] = []
            for number, row in enumerate(reader, start=1):
                batch.append(
                    (
                        row["orf_id"], row["transcript_id"], row["category"],
                        row["start_codon"], row["stop_codon"],
                        row["completeness"], row["chrom"], row["strand"],
                        int(row["nt_length"]), int(row["aa_length"]),
                        row["exon_starts"], row["exon_ends"],
                    )
                )
                if len(batch) >= SQL_BATCH_SIZE:
                    connection.executemany(
                        "INSERT INTO source_batch VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                        batch,
                    )
                    connection.execute(
                        """
                        INSERT OR REPLACE INTO geometry
                        SELECT b.* FROM source_batch b
                        INNER JOIN representatives r ON r.orf_id=b.orf_id
                        """
                    )
                    connection.execute("DELETE FROM source_batch")
                    connection.commit(); batch.clear()
                    progress_print(f"source ORFs: {number:,}")
            if batch:
                connection.executemany(
                    "INSERT INTO source_batch VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                connection.execute(
                    """
                    INSERT OR REPLACE INTO geometry
                    SELECT b.* FROM source_batch b
                    INNER JOIN representatives r ON r.orf_id=b.orf_id
                    """
                )
                connection.execute("DELETE FROM source_batch")
                connection.commit()

        missing = connection.execute(
            """
            SELECT COUNT(*) FROM representatives r
            LEFT JOIN geometry g ON g.orf_id=r.orf_id
            WHERE g.orf_id IS NULL
            """
        ).fetchone()[0]
        if missing:
            raise ValueError(
                f"Missing source coordinates for {missing:,} family representatives."
            )
        connection.execute("CREATE INDEX geometry_orf_idx ON geometry(orf_id)")
        connection.execute("ANALYZE")
        connection.commit()
        counts = {
            "families": connection.execute("SELECT COUNT(*) FROM families").fetchone()[0],
            "representatives": connection.execute(
                "SELECT COUNT(*) FROM representatives"
            ).fetchone()[0],
        }
        marker.write_text(
            json.dumps({"signature": signature, **counts}, indent=2),
            encoding="utf-8",
        )
        logger.write(
            "Family index completed: families={families:,}, representatives={representatives:,}.".format(
                **counts
            )
        )
    finally:
        connection.close()


def _density_cache_track_directory(work_directory: Path, track: DensityTrack) -> Path:
    """Return one track cache directory."""
    return work_directory / "density" / _safe_name(
        f"{track.sample}|{track.strand}|{Path(track.path).resolve()}"
    )


def build_density_cache(
    track: DensityTrack,
    work_directory: Path,
    signature: str,
) -> None:
    """Parse one density file once and save chromosome NumPy arrays."""
    directory = _density_cache_track_directory(work_directory, track)
    marker = directory / "complete.json"
    if marker.is_file():
        state = json.loads(marker.read_text(encoding="utf-8"))
        if state.get("signature") == signature:
            return
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True, exist_ok=True)
    sort_directory = directory / "sort"
    prepared = prepare_density(
        path=track.path,
        file_format=track.file_format,
        work_directory=sort_directory,
    )
    chrom_index: dict[str, str] = {}
    try:
        for number, density in enumerate(
            iter_density_chromosomes(
                prepared.path,
                "bedgraph" if prepared.was_sorted else track.file_format,
            ),
            start=1,
        ):
            file_name = _safe_name(density.chrom) + ".npy"
            output_path = directory / file_name
            array = np.empty(
                density.starts.size,
                dtype=[("start", "<i8"), ("end", "<i8"), ("value", "<f4")],
            )
            array["start"] = density.starts
            array["end"] = density.ends
            array["value"] = density.values
            np.save(output_path, array, allow_pickle=False)
            chrom_index[density.chrom] = file_name
            progress_print(
                f"density cache {track.sample}/{track.strand}: chromosome {number}"
            )
    finally:
        prepared.cleanup()
        shutil.rmtree(sort_directory, ignore_errors=True)
    marker.write_text(
        json.dumps(
            {
                "signature": signature,
                "sample": track.sample,
                "group": track.group,
                "strand": track.strand,
                "chromosomes": chrom_index,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def _load_cached_density(
    work_directory: Path,
    track: DensityTrack,
    chrom: str,
) -> ChromDensity | None:
    """Load one chromosome density as memory-mapped arrays."""
    directory = _density_cache_track_directory(work_directory, track)
    marker_path = directory / "complete.json"
    if not marker_path.is_file():
        raise FileNotFoundError(marker_path)
    state = json.loads(marker_path.read_text(encoding="utf-8"))
    file_name = state.get("chromosomes", {}).get(chrom)
    if not file_name:
        return None
    array = np.load(directory / file_name, mmap_mode="r", allow_pickle=False)
    starts = array["start"]
    ends = array["end"]
    values = array["value"]
    return ChromDensity(
        chrom=chrom,
        starts=starts,
        ends=ends,
        values=values,
        max_end=int(ends[-1]) if ends.size else 0,
    )


def _oriented_blocks(
    starts: Sequence[int],
    ends: Sequence[int],
    strand: str,
) -> tuple[tuple[int, int], ...]:
    """Return exon blocks in transcript order."""
    blocks = tuple(sorted(zip(starts, ends), key=lambda item: item[0]))
    return blocks if strand == "+" else tuple(reversed(blocks))


def _is_suffix(
    scaffold_starts: Sequence[int],
    scaffold_ends: Sequence[int],
    rep_starts: Sequence[int],
    rep_ends: Sequence[int],
    strand: str,
) -> bool:
    """Return whether a representative is an exact spliced suffix."""
    scaffold = _oriented_blocks(scaffold_starts, scaffold_ends, strand)
    representative = _oriented_blocks(rep_starts, rep_ends, strand)
    if not representative or len(representative) > len(scaffold):
        return False
    tail = scaffold[-len(representative):]
    if len(representative) > 1 and tail[1:] != representative[1:]:
        return False
    long_first = tail[0]
    short_first = representative[0]
    if strand == "+":
        return short_first[1] == long_first[1] and short_first[0] >= long_first[0]
    return short_first[0] == long_first[0] and short_first[1] <= long_first[1]


def _row_to_mapping(row: sqlite3.Row) -> dict[str, Any]:
    return {key: row[key] for key in row.keys()}


def _build_family_geometry(rows: Sequence[sqlite3.Row]) -> FamilyGeometry | InvalidFamily:
    """Compile one family and convert incompatible geometry to a record."""
    first = rows[0]
    base = {
        "family_id": first["family_id"],
        "gene_id": first["gene_id"],
        "chrom": first["chrom"],
        "strand": first["strand"],
        "category": first["family_category"],
        "family_type": first["family_type"],
        "family_size": int(first["family_size"]),
        "structural_primary": first["structural_primary"],
    }
    try:
        parsed: list[dict[str, Any]] = []
        for row in rows:
            record = _row_to_mapping(row)
            starts = _split_ints(record["exon_starts"])
            ends = _split_ints(record["exon_ends"])
            if not starts or len(starts) != len(ends):
                raise ValueError(f"invalid exon blocks for {record['orf_id']}")
            if any(end <= start for start, end in zip(starts, ends)):
                raise ValueError(f"invalid exon interval for {record['orf_id']}")
            record["starts"] = starts
            record["ends"] = ends
            record["coding_nt_length"] = _coding_nt_length(record)
            parsed.append(record)

        scaffold = max(
            parsed,
            key=lambda record: (
                int(record["nt_length"]),
                record["orf_id"] == base["structural_primary"],
            ),
        )
        representatives: list[RepresentativeGeometry] = []
        scaffold_nt_length = int(scaffold["nt_length"])
        for record in parsed:
            difference = scaffold_nt_length - int(record["nt_length"])
            if difference < 0 or difference % 3 != 0:
                raise ValueError(
                    f"phase-incompatible representative {record['orf_id']}"
                )
            if not _is_suffix(
                scaffold["starts"], scaffold["ends"],
                record["starts"], record["ends"], base["strand"],
            ):
                raise ValueError(
                    f"non-suffix representative {record['orf_id']}"
                )
            representatives.append(
                RepresentativeGeometry(
                    orf_id=record["orf_id"],
                    transcript_id=record["transcript_id"],
                    category=record["category"],
                    start_codon=record["start_codon"],
                    nt_length=int(record["nt_length"]),
                    coding_nt_length=int(record["coding_nt_length"]),
                    aa_length=int(record["aa_length"]),
                    starts=record["starts"],
                    ends=record["ends"],
                    start_offset=difference,
                )
            )
        representatives.sort(key=lambda item: (item.start_offset, item.orf_id))
        common = max(representatives, key=lambda item: item.start_offset)
        coding_end = int(scaffold["coding_nt_length"])
        if common.start_offset >= coding_end:
            raise ValueError("empty common-body coding segment")
        return FamilyGeometry(
            **base,
            scaffold_orf_id=scaffold["orf_id"],
            scaffold_nt_length=scaffold_nt_length,
            scaffold_coding_nt_length=coding_end,
            scaffold_starts=scaffold["starts"],
            scaffold_ends=scaffold["ends"],
            common_body_orf=common.orf_id,
            common_body_start=common.start_offset,
            common_body_end=coding_end,
            representatives=tuple(representatives),
        )
    except Exception as error:
        return InvalidFamily(**base, reason=str(error))


def _iter_chromosome_family_batches(
    database_path: str | Path,
    chromosome: str,
    batch_size: int = FAMILY_BATCH_SIZE,
) -> Iterator[list[FamilyGeometry | InvalidFamily]]:
    """Stream compiled family geometry from SQLite in bounded batches."""
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    query = """
        SELECT
            f.family_id, f.gene_id, f.chrom, f.strand,
            f.category AS family_category, f.family_type, f.family_size,
            f.structural_primary, r.orf_id, r.is_primary,
            g.transcript_id, g.category, g.start_codon, g.stop_codon,
            g.completeness, g.nt_length, g.aa_length,
            g.exon_starts, g.exon_ends
        FROM families f
        JOIN representatives r ON r.family_id=f.family_id
        JOIN geometry g ON g.orf_id=r.orf_id
        WHERE f.chrom=?
        ORDER BY f.family_id, g.nt_length DESC, r.orf_id
    """
    try:
        current_id: str | None = None
        current_rows: list[sqlite3.Row] = []
        batch: list[FamilyGeometry | InvalidFamily] = []
        for row in connection.execute(query, (chromosome,)):
            family_id = row["family_id"]
            if current_id is None:
                current_id = family_id
            elif family_id != current_id:
                batch.append(_build_family_geometry(current_rows))
                current_rows = []
                current_id = family_id
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
            current_rows.append(row)
        if current_rows:
            batch.append(_build_family_geometry(current_rows))
        if batch:
            yield batch
    finally:
        connection.close()


def _chromosomes(database_path: str | Path) -> list[str]:
    """Return chromosomes in deterministic lexical-natural order."""
    connection = sqlite3.connect(database_path)
    try:
        values = [row[0] for row in connection.execute(
            "SELECT DISTINCT chrom FROM families"
        )]
    finally:
        connection.close()

    def key(value: str) -> tuple[Any, ...]:
        pieces: list[Any] = []
        token = ""
        numeric = value[:1].isdigit()
        for character in value:
            now_numeric = character.isdigit()
            if token and now_numeric != numeric:
                pieces.append(int(token) if numeric else token)
                token = ""
            token += character
            numeric = now_numeric
        if token:
            pieces.append(int(token) if numeric else token)
        return tuple(pieces)

    return sorted(values, key=key)


def _track_for_family(
    tracks_by_sample: Mapping[str, Mapping[str, DensityTrack]],
    sample: str,
    strand: str,
) -> DensityTrack | None:
    """Select a strand-specific or unstranded track."""
    mapping = tracks_by_sample[sample]
    return mapping.get(strand) or mapping.get(".")


def _extract_sparse_signals(
    families: Sequence[FamilyGeometry],
    density: ChromDensity | None,
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Map one chromosome density to family scaffold transcript offsets."""
    if density is None or density.starts.size == 0 or not families:
        return {}

    block_starts: list[int] = []
    block_ends: list[int] = []
    block_offsets: list[int] = []
    block_families: list[int] = []
    block_strands: list[str] = []
    for family_index, family in enumerate(families):
        blocks = _oriented_blocks(
            family.scaffold_starts,
            family.scaffold_ends,
            family.strand,
        )
        offset = 0
        for start, end in blocks:
            block_starts.append(start)
            block_ends.append(end)
            block_offsets.append(offset)
            block_families.append(family_index)
            block_strands.append(family.strand)
            offset += end - start

    starts_array = np.asarray(block_starts, dtype=np.int64)
    ends_array = np.asarray(block_ends, dtype=np.int64)
    left = np.searchsorted(density.ends, starts_array, side="right")
    right = np.searchsorted(density.starts, ends_array, side="left")

    positions: dict[int, list[int]] = defaultdict(list)
    values: dict[int, list[float]] = defaultdict(list)
    for block_index in np.flatnonzero(right > left):
        family_index = block_families[int(block_index)]
        block_start = block_starts[int(block_index)]
        block_end = block_ends[int(block_index)]
        block_offset = block_offsets[int(block_index)]
        strand = block_strands[int(block_index)]
        for density_index in range(int(left[block_index]), int(right[block_index])):
            value = abs(float(density.values[density_index]))
            if value <= 0 or not math.isfinite(value):
                continue
            overlap_start = max(block_start, int(density.starts[density_index]))
            overlap_end = min(block_end, int(density.ends[density_index]))
            if overlap_end <= overlap_start:
                continue
            if strand == "+":
                tx_start = block_offset + overlap_start - block_start
                tx_end = block_offset + overlap_end - block_start
            else:
                tx_start = block_offset + block_end - overlap_end
                tx_end = block_offset + block_end - overlap_start
            span = tx_end - tx_start
            if span <= 0:
                continue
            positions[family_index].extend(range(tx_start, tx_end))
            values[family_index].extend([value] * span)

    output: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for family_index, family_positions in positions.items():
        order = np.argsort(np.asarray(family_positions, dtype=np.int64), kind="stable")
        output[family_index] = (
            np.asarray(family_positions, dtype=np.int64)[order],
            np.asarray(values[family_index], dtype=np.float64)[order],
        )
    return output


def segment_features(
    positions: np.ndarray,
    values: np.ndarray,
    start: int,
    end: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> SegmentFeatures:
    """Calculate sufficient statistics for one scaffold segment."""
    length = max(0, int(end) - int(start))
    codon_count = length // 3
    if codon_count <= 0 or positions.size == 0:
        return SegmentFeatures(
            0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0,
        )
    left = int(np.searchsorted(positions, start, side="left"))
    right = int(np.searchsorted(positions, end, side="left"))
    if right <= left:
        return SegmentFeatures(
            0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0,
        )

    local_positions = positions[left:right] - start
    local_values = values[left:right]
    valid = local_positions < codon_count * 3
    local_positions = local_positions[valid]
    local_values = local_values[valid]
    if local_positions.size == 0:
        return SegmentFeatures(
            0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0,
        )

    frames = local_positions % 3
    codons = local_positions // 3
    rpf_sum = float(local_values.sum())
    frame_density = np.bincount(frames, weights=local_values, minlength=3)
    total_frame = float(frame_density.sum())
    frame0_ratio = float(frame_density[0] / total_frame) if total_frame else 0.0
    codon_profile = np.bincount(
        codons,
        weights=local_values,
        minlength=codon_count,
    ).astype(np.float64, copy=False)
    covered = int(np.count_nonzero(codon_profile > 0))
    coverage_ratio = covered / codon_count
    nonzero_codons = np.flatnonzero(codon_profile > 0)
    signal_span = (
        (int(nonzero_codons[-1]) - int(nonzero_codons[0]) + 1) / codon_count
        if nonzero_codons.size else 0.0
    )

    boundary = min(config.boundary_codons, codon_count)
    start_rpf = float(codon_profile[:boundary].sum())
    end_rpf = float(codon_profile[-boundary:].sum())
    body_start = boundary
    body_end = max(body_start, codon_count - boundary)
    body_rpf = float(codon_profile[body_start:body_end].sum())

    supported_windows = 0
    supported_starts: list[int] = []
    top_window = 0.0
    if codon_count >= config.window_codons:
        prefix = np.empty(codon_count + 1, dtype=np.float64)
        prefix[0] = 0.0
        np.cumsum(codon_profile, out=prefix[1:])
        covered_prefix = np.empty(codon_count + 1, dtype=np.int64)
        covered_prefix[0] = 0
        np.cumsum((codon_profile > 0).astype(np.int64), out=covered_prefix[1:])
        frame_profiles = [
            np.bincount(
                codons[frames == frame],
                weights=local_values[frames == frame],
                minlength=codon_count,
            )
            for frame in range(3)
        ]
        frame_prefixes = []
        for profile in frame_profiles:
            current = np.empty(codon_count + 1, dtype=np.float64)
            current[0] = 0.0
            np.cumsum(profile, out=current[1:])
            frame_prefixes.append(current)
        starts = list(range(0, codon_count - config.window_codons + 1, config.window_step_codons))
        last_start = codon_count - config.window_codons
        if not starts or starts[-1] != last_start:
            starts.append(last_start)
        for window_start in starts:
            window_end = window_start + config.window_codons
            window_sum = float(prefix[window_end] - prefix[window_start])
            top_window = max(top_window, window_sum)
            window_covered = int(
                covered_prefix[window_end] - covered_prefix[window_start]
            )
            window_frames = [
                float(item[window_end] - item[window_start])
                for item in frame_prefixes
            ]
            window_total = sum(window_frames)
            window_frame0 = window_frames[0] / window_total if window_total else 0.0
            if (
                window_sum >= thresholds.min_window_rpf
                and window_covered >= thresholds.min_window_covered
                and window_frame0 >= thresholds.moderate_periodicity
            ):
                supported_windows += 1
                supported_starts.append(window_start)
    top_window_fraction = top_window / rpf_sum if rpf_sum else 0.0
    distributed = 0
    if supported_starts:
        selected = [supported_starts[0]]
        for value in supported_starts[1:]:
            if value - selected[-1] >= config.min_window_gap_codons:
                selected.append(value)
        distributed = len(selected)
    localized_only = (
        codon_count >= config.long_min_codons
        and signal_span <= config.localized_span_max
        and top_window_fraction >= config.localized_top_window_fraction
    )
    abundance_score = min(1.0, rpf_sum / max(thresholds.min_rpf_sum, 1e-9))
    coverage_score = min(1.0, coverage_ratio / max(thresholds.min_coverage_ratio, 1e-9))
    periodicity_score = max(0.0, min(1.0, (frame0_ratio - 1.0 / 3.0) / (2.0 / 3.0)))
    distribution_score = min(1.0, max(signal_span, distributed / max(config.min_supported_windows, 1)))
    score = (
        0.30 * abundance_score
        + 0.25 * coverage_score
        + 0.30 * periodicity_score
        + 0.15 * distribution_score
    )
    return SegmentFeatures(
        rpf_sum=rpf_sum,
        rpf_per_codon=rpf_sum / codon_count,
        covered_codon=covered,
        coverage_ratio=coverage_ratio,
        frame0_density=float(frame_density[0]),
        frame1_density=float(frame_density[1]),
        frame2_density=float(frame_density[2]),
        frame0_ratio=frame0_ratio,
        supported_windows=supported_windows,
        distributed_windows=distributed,
        signal_span=signal_span,
        top_window_fraction=top_window_fraction,
        start_rpf=start_rpf,
        body_rpf=body_rpf,
        end_rpf=end_rpf,
        localized_only=localized_only,
        score=score,
    )


def classify_sample(
    features: SegmentFeatures,
    codon_count: int,
    thresholds: Thresholds,
    config: EngineConfig,
) -> tuple[str, str]:
    """Classify one family common body in one independent sample."""
    if features.rpf_sum <= 0:
        return "NoEvidence", "zero_common_body_rpf"
    core_pass = (
        features.rpf_sum >= thresholds.min_rpf_sum
        and features.rpf_per_codon >= thresholds.min_rpf_per_codon
        and features.covered_codon >= thresholds.min_covered_codon
        and features.coverage_ratio >= thresholds.min_coverage_ratio
    )
    if not core_pass:
        return "LowConfidence", "abundance_or_coverage_below_threshold"
    if features.frame0_ratio < thresholds.moderate_periodicity:
        return "LowConfidence", "periodicity_below_calibrated_threshold"
    if features.localized_only:
        return "LowConfidence", "localized_only_long_orf_signal"

    if codon_count <= config.short_max_codons:
        support = True
    elif codon_count < config.long_min_codons:
        support = (
            features.supported_windows >= 1
            or features.coverage_ratio >= max(0.20, thresholds.min_coverage_ratio * 1.5)
        )
    else:
        support = (
            features.distributed_windows >= config.min_supported_windows
            or features.signal_span >= config.min_signal_span
            or (features.start_rpf > 0 and features.end_rpf > 0)
        )
    if not support:
        return "LowConfidence", "whole_orf_distribution_not_supported"

    high = (
        features.frame0_ratio >= thresholds.strong_periodicity
        and features.coverage_ratio >= max(thresholds.min_coverage_ratio, 0.20)
        and (
            codon_count <= config.short_max_codons
            or features.distributed_windows >= config.min_supported_windows
            or features.signal_span >= config.min_signal_span
        )
    )
    return (
        ("HighConfidence", "strong_replicable_common_body_evidence")
        if high
        else ("MediumConfidence", "moderate_common_body_evidence")
    )


def classify_extension(
    features: SegmentFeatures,
    thresholds: Thresholds,
) -> bool:
    """Return whether one alternative-start-specific segment is supported."""
    return (
        features.rpf_sum >= max(1.0, thresholds.min_window_rpf * 0.5)
        and features.covered_codon >= 1
        and features.coverage_ratio >= 0.10
        and features.frame0_ratio >= thresholds.moderate_periodicity
    )


def _manual_thresholds(
    sample: str,
    config: EngineConfig,
) -> Thresholds:
    """Return the eight user-defined manual evidence thresholds."""
    return Thresholds(
        sample=sample,
        source="manual",
        control_count=0,
        min_rpf_sum=config.min_rpf_sum,
        min_rpf_per_codon=config.min_rpf_per_codon,
        min_covered_codon=config.min_covered_codon,
        min_coverage_ratio=config.min_coverage_ratio,
        moderate_periodicity=config.moderate_periodicity,
        strong_periodicity=config.strong_periodicity,
        min_window_rpf=config.min_window_rpf,
        min_window_covered=config.min_window_covered,
    )




def _normalized_length_category(category: str) -> str | None:
    """Map target smORF categories to uORF, dORF, or lncORF."""
    value = str(category).strip().lower()
    for target in TARGET_LENGTH_CATEGORIES:
        if target in value:
            return target
    return None


def resolve_length_models(
    database_path: Path,
    config: EngineConfig,
    logger: _StageLogger,
) -> EngineConfig:
    """Resolve robust short and long cutoffs from target-smORF distributions.

    Equal-tertile cutoffs are inappropriate for smORFs because their length
    distribution is strongly right-skewed. The short boundary is estimated
    from the lower quartile and constrained to 10-30 codons. The long boundary
    is estimated from the upper decile and constrained to 60-100 codons.
    """
    by_category: dict[str, list[int]] = defaultdict(list)
    connection = sqlite3.connect(database_path)
    try:
        query = """
            SELECT f.category, g.nt_length, g.stop_codon, g.completeness
            FROM families f
            JOIN geometry g ON g.orf_id=f.structural_primary
        """
        for category, nt_length, stop_codon, completeness in connection.execute(query):
            normalized = _normalized_length_category(category)
            if normalized is None:
                continue
            coding = _coding_nt_length(
                {
                    "nt_length": nt_length,
                    "stop_codon": stop_codon,
                    "completeness": completeness,
                }
            )
            codons = coding // 3
            if codons > 0:
                by_category[normalized].append(codons)
    finally:
        connection.close()

    short_values: list[float] = []
    long_values: list[float] = []
    for category in sorted(TARGET_LENGTH_CATEGORIES):
        lengths = by_category.get(category, [])
        if not lengths:
            continue
        short_values.append(float(np.quantile(lengths, 0.25)))
        long_values.append(float(np.quantile(lengths, 0.90)))

    short_cutoff = config.short_max_codons
    long_cutoff = config.long_min_codons
    if "short_max_codons" not in config.hidden_overrides and short_values:
        short_cutoff = int(
            np.clip(
                round(float(np.median(short_values))),
                10,
                30,
            )
        )
    if "long_min_codons" not in config.hidden_overrides and long_values:
        long_cutoff = int(
            np.clip(
                round(float(np.median(long_values))),
                60,
                100,
            )
        )
    if long_cutoff <= short_cutoff:
        long_cutoff = max(60, short_cutoff + 1)

    logger.write(
        "Length models: short<=%d codons, medium=%d-%d codons, long>=%d codons."
        % (
            short_cutoff,
            short_cutoff + 1,
            long_cutoff - 1,
            long_cutoff,
        )
    )
    return replace(
        config,
        short_max_codons=short_cutoff,
        long_min_codons=long_cutoff,
    )



def _longest_covered_run(
    positions: np.ndarray,
    coding_nt_length: int,
) -> int:
    """Return the longest consecutive covered-codon run."""
    valid = positions[(positions >= 0) & (positions < coding_nt_length)]
    if valid.size == 0:
        return 0
    codons = np.unique(valid // 3)
    if codons.size == 0:
        return 0
    breaks = np.flatnonzero(np.diff(codons) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks + 1, [codons.size]))
    return int(np.max(ends - starts))


def _select_expression_controls(
    records: list[dict[str, float | int]],
    config: EngineConfig,
) -> list[dict[str, float | int]]:
    """Select annotated_ORFs using the lower expression quantile."""
    if len(records) < config.positive_min_controls:
        return []
    ordered = sorted(
        records,
        key=lambda item: float(item["rpf_per_codon"]),
    )
    expression = np.asarray(
        [float(item["rpf_per_codon"]) for item in ordered],
        dtype=np.float64,
    )
    cutoff = float(np.quantile(expression, config.positive_quantile))
    selected = [
        item
        for item in ordered
        if float(item["rpf_per_codon"]) >= cutoff
    ]
    if len(selected) < config.positive_min_controls:
        selected = ordered[-config.positive_min_controls:]
    if len(selected) > config.positive_max_controls:
        indices = np.linspace(
            0,
            len(selected) - 1,
            config.positive_max_controls,
            dtype=np.int64,
        )
        selected = [selected[int(index)] for index in indices]
    return selected


def _derive_canonical_model(
    config: EngineConfig,
    records: Sequence[Mapping[str, float | int]],
    logger: _StageLogger,
) -> EngineConfig:
    """Derive non-stricter hidden parameters from annotated_ORF controls.

    Canonical controls may relax the balanced hidden defaults when library
    quality is poor, but must never make the smORF whole-ORF model stricter.
    """
    if not records:
        return config

    quantile = config.positive_quantile
    longest_runs = np.asarray(
        [float(item["longest_run"]) for item in records],
        dtype=np.float64,
    )
    signal_spans = np.asarray(
        [float(item["signal_span"]) for item in records],
        dtype=np.float64,
    )
    top_fractions = np.asarray(
        [float(item["top_window_fraction"]) for item in records],
        dtype=np.float64,
    )

    default_window = int(ADVANCED_DEFAULTS["window_codons"])
    window_codons = config.window_codons
    if "window_codons" not in config.hidden_overrides:
        window_codons = int(
            np.clip(
                round(float(np.quantile(longest_runs, quantile))),
                12,
                default_window,
            )
        )

    window_step = config.window_step_codons
    if "window_step_codons" not in config.hidden_overrides:
        window_step = max(
            1,
            min(
                int(ADVANCED_DEFAULTS["window_step_codons"]),
                int(round(window_codons * 0.25)),
            ),
        )

    window_gap = config.min_window_gap_codons
    if "min_window_gap_codons" not in config.hidden_overrides:
        window_gap = max(
            window_step,
            min(
                int(ADVANCED_DEFAULTS["min_window_gap_codons"]),
                int(round(window_codons * 0.75)),
            ),
        )

    supported_windows = config.min_supported_windows
    if "min_supported_windows" not in config.hidden_overrides:
        # Keep the validated balanced requirement. Annotated mORFs are too
        # long to provide a transferable absolute window-count threshold.
        supported_windows = int(
            ADVANCED_DEFAULTS["min_supported_windows"]
        )

    signal_span = config.min_signal_span
    if "min_signal_span" not in config.hidden_overrides:
        signal_span = float(
            np.clip(
                min(
                    float(ADVANCED_DEFAULTS["min_signal_span"]),
                    float(np.quantile(signal_spans, quantile)),
                ),
                0.20,
                float(ADVANCED_DEFAULTS["min_signal_span"]),
            )
        )

    localized_span = config.localized_span_max
    if "localized_span_max" not in config.hidden_overrides:
        localized_span = float(
            np.clip(
                min(
                    float(ADVANCED_DEFAULTS["localized_span_max"]),
                    float(
                        np.quantile(
                            signal_spans,
                            max(0.05, quantile / 2.0),
                        )
                    ),
                ),
                0.08,
                float(ADVANCED_DEFAULTS["localized_span_max"]),
            )
        )

    top_fraction = config.localized_top_window_fraction
    if "localized_top_window_fraction" not in config.hidden_overrides:
        # A larger threshold makes localized-only rejection more conservative.
        top_fraction = float(
            np.clip(
                max(
                    float(
                        ADVANCED_DEFAULTS[
                            "localized_top_window_fraction"
                        ]
                    ),
                    float(np.quantile(top_fractions, 1.0 - quantile)),
                ),
                float(
                    ADVANCED_DEFAULTS[
                        "localized_top_window_fraction"
                    ]
                ),
                0.90,
            )
        )

    boundary = config.boundary_codons
    if "boundary_codons" not in config.hidden_overrides:
        boundary = max(
            3,
            min(
                int(ADVANCED_DEFAULTS["boundary_codons"]),
                int(round(window_codons * 0.25)),
            ),
        )

    resolved = replace(
        config,
        window_codons=window_codons,
        window_step_codons=window_step,
        min_supported_windows=supported_windows,
        min_window_gap_codons=window_gap,
        min_signal_span=signal_span,
        localized_span_max=min(
            localized_span,
            signal_span * 0.90,
        ),
        localized_top_window_fraction=top_fraction,
        boundary_codons=boundary,
    )
    logger.write(
        "Canonical model: window=%d, step=%d, supported_windows=%d, gap=%d, "
        "signal_span=%.3f, localized_span=%.3f, top_window_fraction=%.3f, "
        "boundary=%d."
        % (
            resolved.window_codons,
            resolved.window_step_codons,
            resolved.min_supported_windows,
            resolved.min_window_gap_codons,
            resolved.min_signal_span,
            resolved.localized_span_max,
            resolved.localized_top_window_fraction,
            resolved.boundary_codons,
        )
    )
    return resolved


def _control_families(
    database_path: str | Path,
    category: str,
    maximum: int,
) -> list[tuple[str, FamilyGeometry]]:
    """Select deterministic annotated-ORF controls."""
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    query = """
        SELECT
            f.family_id, f.gene_id, f.chrom, f.strand,
            f.category AS family_category, f.family_type, f.family_size,
            f.structural_primary, r.orf_id, r.is_primary,
            g.transcript_id, g.category, g.start_codon, g.stop_codon,
            g.completeness, g.nt_length, g.aa_length,
            g.exon_starts, g.exon_ends
        FROM families f
        JOIN representatives r ON r.family_id=f.family_id
        JOIN geometry g ON g.orf_id=r.orf_id
        WHERE f.category=?
        ORDER BY f.family_id, g.nt_length DESC
        LIMIT ?
    """
    output: list[tuple[str, FamilyGeometry]] = []
    try:
        current: str | None = None
        rows: list[sqlite3.Row] = []
        for row in connection.execute(query, (category, maximum * 4)):
            family_id = row["family_id"]
            if current is None:
                current = family_id
            elif family_id != current:
                geometry = _build_family_geometry(rows)
                if isinstance(geometry, FamilyGeometry):
                    output.append((geometry.chrom, geometry))
                    if len(output) >= maximum:
                        break
                rows = []
                current = family_id
            rows.append(row)
        if rows and len(output) < maximum:
            geometry = _build_family_geometry(rows)
            if isinstance(geometry, FamilyGeometry):
                output.append((geometry.chrom, geometry))
    finally:
        connection.close()
    return output


def calibrate_samples(
    database_path: Path,
    work_directory: Path,
    tracks: Sequence[DensityTrack],
    config: EngineConfig,
    logger: _StageLogger,
) -> tuple[list[Thresholds], EngineConfig]:
    """Resolve manual thresholds or derive canonical thresholds from annotated_ORFs."""
    samples = sorted({track.sample for track in tracks})
    output_path = work_directory / "calibration.json"

    if config.evidence_mode == "manual":
        thresholds = [
            _manual_thresholds(sample, config)
            for sample in samples
        ]
        output_path.write_text(
            json.dumps(
                {
                    "mode": "manual",
                    "thresholds": [asdict(value) for value in thresholds],
                    "model": asdict(config),
                },
                indent=2,
                default=list,
            ),
            encoding="utf-8",
        )
        logger.write(
            "Manual evidence mode: use the eight user-defined thresholds."
        )
        return thresholds, config

    controls = _control_families(
        database_path,
        "annotated_ORF",
        max(
            config.positive_min_controls,
            config.positive_max_controls * 2,
        ),
    )
    if len(controls) < config.positive_min_controls:
        raise ValueError(
            "Canonical mode requires annotated_ORF controls, but only "
            f"{len(controls):,} valid controls were found. Re-run "
            "smorf_cluster with annotated_ORF retained in --keep-categories."
        )

    tracks_by_sample: dict[str, dict[str, DensityTrack]] = defaultdict(dict)
    for track in tracks:
        tracks_by_sample[track.sample][track.strand] = track
    controls_by_chrom: dict[str, list[FamilyGeometry]] = defaultdict(list)
    for chrom, geometry in controls:
        controls_by_chrom[chrom].append(geometry)

    provisional = _manual_thresholds("canonical_provisional", config)
    records_by_sample: dict[str, list[dict[str, float | int]]] = {}
    all_selected: list[dict[str, float | int]] = []

    for sample_number, sample in enumerate(samples, start=1):
        records: list[dict[str, float | int]] = []
        mapping = tracks_by_sample[sample]
        for chrom, families in controls_by_chrom.items():
            plus_track = mapping.get("+") or mapping.get(".")
            minus_track = mapping.get("-") or mapping.get(".")
            density_by_strand = {
                "+": (
                    _load_cached_density(
                        work_directory,
                        plus_track,
                        chrom,
                    )
                    if plus_track
                    else None
                ),
                "-": (
                    _load_cached_density(
                        work_directory,
                        minus_track,
                        chrom,
                    )
                    if minus_track
                    else None
                ),
            }
            for strand in ("+", "-"):
                subset = [
                    family
                    for family in families
                    if family.strand == strand
                ]
                if not subset:
                    continue
                signals = _extract_sparse_signals(
                    subset,
                    density_by_strand[strand],
                )
                for index, family in enumerate(subset):
                    signal = signals.get(index)
                    if signal is None:
                        continue
                    features = segment_features(
                        signal[0],
                        signal[1],
                        0,
                        family.scaffold_coding_nt_length,
                        provisional,
                        config,
                    )
                    if features.rpf_sum <= 0:
                        continue
                    codon_count = max(
                        1,
                        family.scaffold_coding_nt_length // 3,
                    )
                    records.append(
                        {
                            "rpf_sum": features.rpf_sum,
                            "rpf_per_codon": features.rpf_per_codon,
                            "covered_codon": features.covered_codon,
                            "coverage_ratio": features.coverage_ratio,
                            "frame0_ratio": features.frame0_ratio,
                            "signal_span": features.signal_span,
                            "top_window_fraction": (
                                features.top_window_fraction
                            ),
                            "codon_count": codon_count,
                            "longest_run": _longest_covered_run(
                                signal[0],
                                family.scaffold_coding_nt_length,
                            ),
                        }
                    )

        selected = _select_expression_controls(records, config)
        if len(selected) < config.positive_min_controls:
            raise ValueError(
                "Canonical mode requires at least "
                f"{config.positive_min_controls} expressed annotated_ORFs "
                f"per sample. Sample {sample} provided only "
                f"{len(selected)} after expression-quantile selection. "
                "Retain annotated_ORF in smorf_cluster or use "
                "--evidence-mode manual."
            )
        records_by_sample[sample] = selected
        all_selected.extend(selected)
        logger.write(
            "Canonical controls sample=%s: raw=%d, selected=%d."
            % (sample, len(records), len(selected))
        )
        progress_print(
            f"canonical controls: {sample_number}/{len(samples)}"
        )

    config = _derive_canonical_model(
        config,
        all_selected,
        logger,
    )
    thresholds: list[Thresholds] = []
    q = config.positive_quantile
    for sample in samples:
        selected = records_by_sample[sample]

        def array(name: str) -> np.ndarray:
            return np.asarray(
                [float(item[name]) for item in selected],
                dtype=np.float64,
            )

        rpf_sum = array("rpf_sum")
        rpf_per_codon = array("rpf_per_codon")
        covered = array("covered_codon")
        coverage = array("coverage_ratio")
        frame0 = array("frame0_ratio")
        expected_window_rpf = rpf_per_codon * config.window_codons
        expected_window_covered = coverage * config.window_codons

        moderate = float(
            np.clip(np.quantile(frame0, q), 0.40, 0.85)
        )
        strong = float(
            np.clip(
                np.median(frame0),
                moderate + 0.05,
                0.95,
            )
        )
        # Canonical controls calibrate the lower data-quality boundary.
        # Thresholds may become more permissive than the balanced defaults,
        # but never stricter, because annotated mORFs are generally much more
        # abundant and longer than smORFs.
        moderate = max(
            0.40,
            min(
                float(MANUAL_THRESHOLD_DEFAULTS["moderate_periodicity"]),
                moderate,
            ),
        )
        strong = max(
            moderate + 0.05,
            min(
                float(MANUAL_THRESHOLD_DEFAULTS["strong_periodicity"]),
                strong,
            ),
        )
        threshold = Thresholds(
            sample=sample,
            source="canonical:annotated_ORF",
            control_count=len(selected),
            min_rpf_sum=max(
                1.0,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_rpf_sum"]),
                    float(np.quantile(rpf_sum, q)),
                ),
            ),
            min_rpf_per_codon=max(
                0.01,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_rpf_per_codon"]),
                    float(np.quantile(rpf_per_codon, q)),
                ),
            ),
            min_covered_codon=max(
                1,
                min(
                    int(MANUAL_THRESHOLD_DEFAULTS["min_covered_codon"]),
                    int(math.floor(float(np.quantile(covered, q)))),
                ),
            ),
            min_coverage_ratio=max(
                0.01,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_coverage_ratio"]),
                    float(np.quantile(coverage, q)),
                ),
            ),
            moderate_periodicity=moderate,
            strong_periodicity=min(0.95, strong),
            min_window_rpf=max(
                1.0,
                min(
                    float(MANUAL_THRESHOLD_DEFAULTS["min_window_rpf"]),
                    float(np.quantile(expected_window_rpf, q)),
                ),
            ),
            min_window_covered=max(
                1,
                min(
                    int(MANUAL_THRESHOLD_DEFAULTS["min_window_covered"]),
                    int(
                        math.floor(
                            float(
                                np.quantile(
                                    expected_window_covered,
                                    q,
                                )
                            )
                        )
                    ),
                ),
            ),
        )
        thresholds.append(threshold)
        logger.write(
            "Canonical thresholds sample=%s: controls=%d, RPF=%.3f, "
            "RPF/codon=%.4f, covered=%d, coverage=%.3f, "
            "moderate=%.3f, strong=%.3f, window_RPF=%.3f, "
            "window_covered=%d."
            % (
                sample,
                threshold.control_count,
                threshold.min_rpf_sum,
                threshold.min_rpf_per_codon,
                threshold.min_covered_codon,
                threshold.min_coverage_ratio,
                threshold.moderate_periodicity,
                threshold.strong_periodicity,
                threshold.min_window_rpf,
                threshold.min_window_covered,
            )
        )

    output_path.write_text(
        json.dumps(
            {
                "mode": "canonical",
                "thresholds": [asdict(value) for value in thresholds],
                "model": asdict(config),
            },
            indent=2,
            default=list,
        ),
        encoding="utf-8",
    )
    return thresholds, config



def _new_accumulator(
    family: FamilyGeometry,
) -> _FamilyAccumulator:
    """Create a bounded accumulator for one family."""
    segment_count = max(0, len(family.representatives) - 1)
    return _FamilyAccumulator(
        any_signal=False,
        any_evidence=False,
        support_count=0,
        high_count=0,
        group_mask=0,
        group_support_counts={},
        group_high_counts={},
        supporting_samples=[],
        best_level="NoEvidence",
        best_score=-1.0,
        best_sample="",
        best_features=None,
        segment_support_counts=[0] * segment_count,
        segment_group_support_counts=[
            {} for _ in range(segment_count)
        ],
        segment_best_scores=[0.0] * segment_count,
        segment_support_samples=[
            [] for _ in range(segment_count)
        ],
    )



def _output_columns() -> tuple[str, ...]:
    """Return focused family output columns."""
    return (
        "family_id", "gene_id", "chrom", "strand", "category",
        "family_type", "family_size", "structural_primary",
        "common_body_orf", "evidence_primary", "evidence_status",
        "family_translation_evidence", "reliability_reason",
        "start_site_status", "start_site_reason", "reliable_group",
        "supported_sample_count", "max_group_sample_support",
        "high_confidence_sample_count", "supporting_group_count",
        "supporting_groups", "supporting_samples", "best_sample",
        "best_sample_evidence", "best_rpf_sum", "best_rpf_per_codon",
        "best_covered_codon", "best_coverage_ratio",
        "best_frame0_ratio", "best_supported_windows",
        "best_distributed_windows", "best_signal_span",
        "best_top_window_fraction", "best_localized_only",
        "representative_count", "representative_orf_ids",
        "representative_start_codons", "invalid_family_reason",
    )



def _resolve_start(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    config: EngineConfig,
) -> tuple[str, str, str]:
    """Resolve a contiguous start chain replicated across independent samples."""
    representatives = family.representatives
    if len(representatives) == 1:
        return (
            representatives[0].orf_id,
            "Resolved",
            "singleton_family",
        )

    selected_index = len(representatives) - 1
    moved = False
    for segment_index in range(
        len(representatives) - 2,
        -1,
        -1,
    ):
        if (
            accumulator.segment_support_counts[segment_index]
            >= config.reliable_sample
        ):
            selected_index = segment_index
            moved = True
        else:
            break
    if moved:
        return (
            representatives[selected_index].orf_id,
            "Resolved",
            "replicated_contiguous_extension_support",
        )
    return (
        "",
        "Ambiguous",
        "common_body_supported_but_replicated_extensions_unresolved",
    )




def _family_output_row(
    family: FamilyGeometry,
    accumulator: _FamilyAccumulator,
    group_names: Sequence[str],
    config: EngineConfig,
) -> dict[str, Any]:
    """Convert one family accumulator into a sample-replicated final row.

    ``--group`` is annotation metadata. Reliability is determined from
    independent sample replication across the whole project, not restricted to
    one group.
    """
    supporting_group_indices = sorted(
        accumulator.group_support_counts
    )
    supporting_groups = ",".join(
        str(group_names[index])
        for index in supporting_group_indices
    )
    max_group_support = max(
        accumulator.group_support_counts.values(),
        default=0,
    )

    if accumulator.support_count >= config.reliable_sample:
        status = "Reliable"
        final_evidence = (
            "HighConfidence"
            if accumulator.high_count >= config.reliable_sample
            else "MediumConfidence"
        )
        reason = (
            f"common_body_supported_by_{accumulator.support_count}_"
            f"independent_samples_across_"
            f"{len(supporting_group_indices)}_groups"
        )
        evidence_primary, start_status, start_reason = _resolve_start(
            family,
            accumulator,
            config,
        )
        reliable_group = (
            str(group_names[supporting_group_indices[0]])
            if len(supporting_group_indices) == 1
            else "cross_group"
        )
    elif not accumulator.any_signal:
        status = "NoEvidence"
        final_evidence = "NoEvidence"
        reliable_group = ""
        reason = "no_common_body_p_site_signal_in_any_sample"
        evidence_primary = ""
        start_status = "NotEvaluated"
        start_reason = "family_translation_not_detected"
    else:
        status = "Uncertain"
        final_evidence = accumulator.best_level
        reliable_group = ""
        reason = (
            "insufficient_independent_sample_replication"
            if accumulator.any_evidence
            else "signal_present_but_no_sample_passed_translation_thresholds"
        )
        evidence_primary = ""
        start_status = "Unresolved"
        start_reason = "family_translation_not_reliably_replicated"

    features = accumulator.best_features
    return {
        "family_id": family.family_id,
        "gene_id": family.gene_id,
        "chrom": family.chrom,
        "strand": family.strand,
        "category": family.category,
        "family_type": family.family_type,
        "family_size": family.family_size,
        "structural_primary": family.structural_primary,
        "common_body_orf": family.common_body_orf,
        "evidence_primary": evidence_primary,
        "evidence_status": status,
        "family_translation_evidence": final_evidence,
        "reliability_reason": reason,
        "start_site_status": start_status,
        "start_site_reason": start_reason,
        "reliable_group": reliable_group,
        "supported_sample_count": accumulator.support_count,
        "max_group_sample_support": max_group_support,
        "high_confidence_sample_count": accumulator.high_count,
        "supporting_group_count": len(supporting_group_indices),
        "supporting_groups": supporting_groups,
        "supporting_samples": ",".join(
            accumulator.supporting_samples
        ),
        "best_sample": accumulator.best_sample,
        "best_sample_evidence": accumulator.best_level,
        "best_rpf_sum": (
            0.0 if features is None else features.rpf_sum
        ),
        "best_rpf_per_codon": (
            0.0 if features is None else features.rpf_per_codon
        ),
        "best_covered_codon": (
            0 if features is None else features.covered_codon
        ),
        "best_coverage_ratio": (
            0.0 if features is None else features.coverage_ratio
        ),
        "best_frame0_ratio": (
            0.0 if features is None else features.frame0_ratio
        ),
        "best_supported_windows": (
            0 if features is None else features.supported_windows
        ),
        "best_distributed_windows": (
            0 if features is None else features.distributed_windows
        ),
        "best_signal_span": (
            0.0 if features is None else features.signal_span
        ),
        "best_top_window_fraction": (
            0.0
            if features is None
            else features.top_window_fraction
        ),
        "best_localized_only": (
            False if features is None else features.localized_only
        ),
        "representative_count": len(family.representatives),
        "representative_orf_ids": ",".join(
            item.orf_id for item in family.representatives
        ),
        "representative_start_codons": ",".join(
            item.start_codon for item in family.representatives
        ),
        "invalid_family_reason": "",
    }




def _invalid_output_row(family: InvalidFamily) -> dict[str, Any]:
    """Convert an invalid family into an explicit uncertain row."""
    return {
        "family_id": family.family_id,
        "gene_id": family.gene_id,
        "chrom": family.chrom,
        "strand": family.strand,
        "category": family.category,
        "family_type": family.family_type,
        "family_size": family.family_size,
        "structural_primary": family.structural_primary,
        "common_body_orf": "",
        "evidence_primary": "",
        "evidence_status": "Uncertain",
        "family_translation_evidence": "InvalidFamilyGeometry",
        "reliability_reason": "family_geometry_validation_failed",
        "start_site_status": "Unresolved",
        "start_site_reason": "invalid_family_geometry",
        "reliable_group": "",
        "supported_sample_count": 0,
        "max_group_sample_support": 0,
        "high_confidence_sample_count": 0,
        "supporting_group_count": 0,
        "supporting_groups": "",
        "supporting_samples": "",
        "best_sample": "",
        "best_sample_evidence": "InvalidFamilyGeometry",
        "best_rpf_sum": 0.0,
        "best_rpf_per_codon": 0.0,
        "best_covered_codon": 0,
        "best_coverage_ratio": 0.0,
        "best_frame0_ratio": 0.0,
        "best_supported_windows": 0,
        "best_distributed_windows": 0,
        "best_signal_span": 0.0,
        "best_top_window_fraction": 0.0,
        "best_localized_only": False,
        "representative_count": 0,
        "representative_orf_ids": "",
        "representative_start_codons": "",
        "invalid_family_reason": family.reason,
    }


def _format_value(value: Any) -> str:
    if isinstance(value, float):
        return format(value, ".10g") if math.isfinite(value) else "NA"
    if isinstance(value, bool):
        return "True" if value else "False"
    return str(value)


def _is_annotated_category(value: Any) -> bool:
    """Return whether a category is reserved for canonical calibration."""
    return str(value).strip().lower() in ANNOTATED_CATEGORIES


def _write_rows(handle: TextIO, rows: Sequence[Mapping[str, Any]]) -> None:
    columns = _output_columns()
    handle.write(
        "".join(
            "\t".join(_format_value(row.get(column, "")) for column in columns)
            + "\n"
            for row in rows
        )
    )


def _process_chromosome(task: ChromosomeTask) -> ChromosomeResult:
    """Process one chromosome in bounded family batches."""
    work_directory = Path(task.work_directory)
    part_directory = work_directory / "parts"
    part_directory.mkdir(parents=True, exist_ok=True)
    safe_chrom = _safe_name(task.chromosome)
    part_path = part_directory / f"{safe_chrom}.part.tsv.gz"
    count_path = part_directory / f"{safe_chrom}.counts.json"
    done_path = part_directory / f"{safe_chrom}.done.json"
    if done_path.is_file() and part_path.is_file() and count_path.is_file():
        values = json.loads(count_path.read_text(encoding="utf-8"))
        return ChromosomeResult(
            chromosome=task.chromosome,
            part_path=str(part_path),
            count_path=str(count_path),
            **values,
        )

    thresholds_by_sample = {value.sample: value for value in task.thresholds}
    tracks_by_sample: dict[str, dict[str, DensityTrack]] = defaultdict(dict)
    sample_groups: dict[str, str] = {}
    for track in task.tracks:
        tracks_by_sample[track.sample][track.strand] = track
        sample_groups[track.sample] = track.group
    samples = sorted(tracks_by_sample)
    group_index = {group: index for index, group in enumerate(task.group_names)}
    if len(group_index) > 63:
        raise ValueError("At most 63 biological groups are supported.")

    # Open each sample/strand chromosome cache once per worker task. The arrays
    # are memory-mapped, so keeping the handles does not load all density data
    # into resident memory and avoids reopening thousands of files per batch.
    density_by_sample_strand: dict[tuple[str, str], ChromDensity | None] = {}
    for sample in samples:
        for strand in ("+", "-"):
            track = _track_for_family(tracks_by_sample, sample, strand)
            density_by_sample_strand[(sample, strand)] = (
                _load_cached_density(work_directory, track, task.chromosome)
                if track is not None
                else None
            )

    temporary = part_path.with_name(part_path.name + ".tmp")
    counts = {
        "total_families": 0,
        "reliable_families": 0,
        "uncertain_families": 0,
        "no_evidence_families": 0,
        "reliable_smorfs": 0,
        "invalid_families": 0,
    }
    with gzip.open(
        temporary,
        "wt",
        encoding="utf-8",
        compresslevel=1,
        newline="",
    ) as output_handle:
        for batch_number, batch in enumerate(
            _iter_chromosome_family_batches(
                task.database_path,
                task.chromosome,
            ),
            start=1,
        ):
            valid_families = [item for item in batch if isinstance(item, FamilyGeometry)]
            invalid_families = [item for item in batch if isinstance(item, InvalidFamily)]
            accumulators = [_new_accumulator(family) for family in valid_families]

            for sample in samples:
                threshold = thresholds_by_sample[sample]
                for strand in ("+", "-"):
                    indices = [
                        index for index, family in enumerate(valid_families)
                        if family.strand == strand
                    ]
                    if not indices:
                        continue
                    subset = [valid_families[index] for index in indices]
                    density = density_by_sample_strand[(sample, strand)]
                    signals = _extract_sparse_signals(subset, density)
                    for subset_index, signal in signals.items():
                        family_index = indices[subset_index]
                        family = valid_families[family_index]
                        accumulator = accumulators[family_index]
                        features = segment_features(
                            signal[0], signal[1],
                            family.common_body_start,
                            family.common_body_end,
                            threshold,
                            task.config,
                        )
                        if features.rpf_sum > 0:
                            accumulator.any_signal = True
                        codon_count = (
                            family.common_body_end - family.common_body_start
                        ) // 3
                        level, _ = classify_sample(
                            features,
                            codon_count,
                            threshold,
                            task.config,
                        )
                        if LEVEL_RANK[level] >= LEVEL_RANK["LowConfidence"]:
                            accumulator.any_evidence = True
                        if (
                            LEVEL_RANK[level] > LEVEL_RANK[accumulator.best_level]
                            or (
                                LEVEL_RANK[level] == LEVEL_RANK[accumulator.best_level]
                                and features.score > accumulator.best_score
                            )
                        ):
                            accumulator.best_level = level
                            accumulator.best_score = features.score
                            accumulator.best_sample = sample
                            accumulator.best_features = features
                        group = sample_groups[sample]
                        group_id = group_index[group]
                        if LEVEL_RANK[level] >= LEVEL_RANK["MediumConfidence"]:
                            accumulator.support_count += 1
                            accumulator.supporting_samples.append(sample)
                            accumulator.group_mask |= 1 << group_id
                            accumulator.group_support_counts[group_id] = (
                                accumulator.group_support_counts.get(group_id, 0)
                                + 1
                            )
                            if level == "HighConfidence":
                                accumulator.high_count += 1
                                accumulator.group_high_counts[group_id] = (
                                    accumulator.group_high_counts.get(group_id, 0)
                                    + 1
                                )

                        representatives = family.representatives
                        for segment_index in range(len(representatives) - 1):
                            segment_start = representatives[segment_index].start_offset
                            segment_end = representatives[segment_index + 1].start_offset
                            extension = segment_features(
                                signal[0], signal[1],
                                segment_start, segment_end,
                                threshold, task.config,
                            )
                            if classify_extension(extension, threshold):
                                accumulator.segment_support_counts[segment_index] += 1
                                accumulator.segment_support_samples[segment_index].append(sample)
                                group_counts = (
                                    accumulator
                                    .segment_group_support_counts[segment_index]
                                )
                                group_counts[group_id] = (
                                    group_counts.get(group_id, 0) + 1
                                )
                                accumulator.segment_best_scores[segment_index] = max(
                                    accumulator.segment_best_scores[segment_index],
                                    extension.score,
                                )

            rows: list[Mapping[str, Any]] = []
            for family, accumulator in zip(valid_families, accumulators):
                row = _family_output_row(
                    family, accumulator, task.group_names, task.config
                )
                rows.append(row)
                counts["total_families"] += 1
                if row["evidence_status"] == "Reliable":
                    counts["reliable_families"] += 1
                    if (
                        row["start_site_status"] == "Resolved"
                        and not _is_annotated_category(row["category"])
                    ):
                        counts["reliable_smorfs"] += 1
                elif row["evidence_status"] == "NoEvidence":
                    counts["no_evidence_families"] += 1
                else:
                    counts["uncertain_families"] += 1
            for invalid in invalid_families:
                rows.append(_invalid_output_row(invalid))
                counts["total_families"] += 1
                counts["uncertain_families"] += 1
                counts["invalid_families"] += 1
            _write_rows(output_handle, rows)
            progress_print(
                f"{task.chromosome}: family batches {batch_number:,}, "
                f"families {counts['total_families']:,}"
            )

    os.replace(temporary, part_path)
    count_path.write_text(json.dumps(counts, indent=2), encoding="utf-8")
    done_path.write_text(
        json.dumps({"engine": ENGINE_VERSION, "chromosome": task.chromosome}),
        encoding="utf-8",
    )
    return ChromosomeResult(
        chromosome=task.chromosome,
        part_path=str(part_path),
        count_path=str(count_path),
        **counts,
    )


def _reliable_matches(
    connection: sqlite3.Connection,
    orf_ids: Sequence[str],
) -> dict[str, int]:
    """Return reliable-index match states for one genePred input batch."""
    unique_ids = tuple(dict.fromkeys(orf_ids))
    matches: dict[str, int] = {}
    query_size = 800
    for start in range(0, len(unique_ids), query_size):
        chunk = unique_ids[start : start + query_size]
        placeholders = ",".join("?" for _ in chunk)
        query = (
            "SELECT orf_id, matched FROM reliable_export "
            f"WHERE orf_id IN ({placeholders})"
        )
        for row in connection.execute(query, chunk):
            matches[str(row["orf_id"])] = int(row["matched"])
    return matches


def _export_reliable_genepred(
    source_path: str | Path,
    temporary_path: Path,
    connection: sqlite3.Connection,
) -> int:
    """Filter scanner genePred records using the validated Reliable ID index."""
    expected = int(
        connection.execute(
            "SELECT COUNT(*) FROM reliable_export"
        ).fetchone()[0]
    )
    if expected == 0:
        temporary_path.write_text("", encoding="utf-8")
        return 0

    duplicate_count = 0
    duplicate_examples: list[str] = []
    batch: list[tuple[str, str, int]] = []
    batch_size = 25_000

    def flush_batch(output_handle: TextIO) -> None:
        nonlocal duplicate_count
        if not batch:
            return
        matches = _reliable_matches(
            connection,
            [orf_id for orf_id, _line, _number in batch],
        )
        newly_matched: set[str] = set()
        for orf_id, raw_line, line_number in batch:
            matched = matches.get(orf_id)
            if matched is None:
                continue
            if matched or orf_id in newly_matched:
                duplicate_count += 1
                if len(duplicate_examples) < 10:
                    duplicate_examples.append(
                        f"{orf_id}@line{line_number}"
                    )
                continue
            if len(raw_line.rstrip("\r\n").split("\t")) < 10:
                raise EvidenceEngineError(
                    "Reliable genePred record has fewer than 10 columns at "
                    f"line {line_number}: {orf_id}"
                )
            output_handle.write(raw_line)
            if not raw_line.endswith(("\n", "\r")):
                output_handle.write("\n")
            newly_matched.add(orf_id)
        if newly_matched:
            connection.executemany(
                "UPDATE reliable_export SET matched=1 WHERE orf_id=?",
                ((orf_id,) for orf_id in newly_matched),
            )
        batch.clear()

    with _smart_open(source_path, "rt") as source, temporary_path.open(
        "w",
        encoding="utf-8",
        buffering=8 * 1024 * 1024,
    ) as target:
        for line_number, raw_line in enumerate(source, start=1):
            if raw_line.startswith("#"):
                flush_batch(target)
                target.write(raw_line)
                continue
            text = raw_line.rstrip("\r\n")
            if not text:
                continue
            orf_id = text.split("\t", 1)[0]
            batch.append((orf_id, raw_line, line_number))
            if len(batch) >= batch_size:
                flush_batch(target)
        flush_batch(target)

    if duplicate_count:
        raise EvidenceEngineError(
            "Scanner genePred contains duplicate records for "
            f"{duplicate_count:,} Reliable ORF occurrence(s). Examples: "
            + ", ".join(duplicate_examples)
        )

    missing_count = int(
        connection.execute(
            "SELECT COUNT(*) FROM reliable_export WHERE matched=0"
        ).fetchone()[0]
    )
    if missing_count:
        examples = [
            str(row[0])
            for row in connection.execute(
                "SELECT orf_id FROM reliable_export "
                "WHERE matched=0 ORDER BY orf_id LIMIT 10"
            )
        ]
        raise EvidenceEngineError(
            f"Scanner genePred is missing {missing_count:,} Reliable ORF(s). "
            "Examples: " + ", ".join(examples)
        )
    return expected


def _merge_outputs(
    output_prefix: Path,
    chromosome_results: Sequence[ChromosomeResult],
    database_path: Path,
    config: EngineConfig,
    sample_count: int,
    effective_workers: int,
    orf_genepred: str | Path,
    logger: _StageLogger | None = None,
) -> EngineResult:
    """Merge evidence and atomically export the reliable ORF annotation."""
    master_path = Path(str(output_prefix) + ".smorf_evidence.txt.gz")
    reliable_path = Path(str(output_prefix) + ".reliable_smorf.txt")
    genepred_path = Path(str(output_prefix) + ".reliable_smorf.genepred")
    summary_path = Path(str(output_prefix) + ".evidence_summary.txt")
    master_tmp = master_path.with_name(master_path.name + ".tmp")
    reliable_tmp = reliable_path.with_name(reliable_path.name + ".tmp")
    genepred_tmp = genepred_path.with_name(genepred_path.name + ".tmp")
    summary_tmp = summary_path.with_name(summary_path.name + ".tmp")
    temporary_paths = (
        master_tmp,
        reliable_tmp,
        genepred_tmp,
        summary_tmp,
    )
    reliable_columns = (
        "family_id", "gene_id", "chrom", "strand", "category",
        "family_type", "family_size", "structural_primary",
        "evidence_primary", "transcript_id", "start_codon", "stop_codon",
        "nt_length", "aa_length", "exon_starts", "exon_ends",
        "supported_sample_count", "supporting_samples", "best_sample",
        "best_sample_evidence", "best_rpf_sum", "best_coverage_ratio",
        "best_frame0_ratio", "reliability_reason", "start_site_reason",
    )
    reliable_smorfs_written = 0
    annotated_controls_excluded = 0
    reliable_genepred_records = 0
    connection: sqlite3.Connection | None = None
    try:
        with gzip.open(
            master_tmp,
            "wt",
            encoding="utf-8",
            compresslevel=1,
        ) as handle:
            handle.write("\t".join(_output_columns()) + "\n")
        with master_tmp.open("ab") as target:
            for result in chromosome_results:
                with Path(result.part_path).open("rb") as source:
                    shutil.copyfileobj(
                        source,
                        target,
                        length=16 * 1024 * 1024,
                    )

        connection = sqlite3.connect(database_path)
        connection.row_factory = sqlite3.Row
        connection.execute(
            """
            CREATE TEMP TABLE reliable_export (
                orf_id TEXT PRIMARY KEY,
                matched INTEGER NOT NULL DEFAULT 0
            ) WITHOUT ROWID
            """
        )
        with reliable_tmp.open(
            "w", encoding="utf-8", buffering=8 * 1024 * 1024
        ) as output_handle:
            output_handle.write("\t".join(reliable_columns) + "\n")
            with gzip.open(master_tmp, "rt", encoding="utf-8") as input_handle:
                reader = csv.DictReader(input_handle, delimiter="\t")
                for row in reader:
                    if not (
                        row["evidence_status"] == "Reliable"
                        and row["start_site_status"] == "Resolved"
                        and row["evidence_primary"]
                    ):
                        continue
                    if _is_annotated_category(row.get("category", "")):
                        annotated_controls_excluded += 1
                        continue
                    geometry = connection.execute(
                        """
                        SELECT transcript_id, category, start_codon, stop_codon,
                               nt_length, aa_length, exon_starts, exon_ends
                        FROM geometry WHERE orf_id=?
                        """,
                        (row["evidence_primary"],),
                    ).fetchone()
                    if geometry is None:
                        raise EvidenceEngineError(
                            "Reliable evidence primary is absent from the "
                            f"geometry index: {row['evidence_primary']}"
                        )
                    if _is_annotated_category(geometry["category"]):
                        annotated_controls_excluded += 1
                        continue
                    try:
                        connection.execute(
                            "INSERT INTO reliable_export (orf_id) VALUES (?)",
                            (row["evidence_primary"],),
                        )
                    except sqlite3.IntegrityError as error:
                        raise EvidenceEngineError(
                            "One ORF is the evidence primary of multiple "
                            "reliable families: "
                            f"{row['evidence_primary']}"
                        ) from error
                    output = {
                        **row,
                        "transcript_id": geometry["transcript_id"],
                        "start_codon": geometry["start_codon"],
                        "stop_codon": geometry["stop_codon"],
                        "nt_length": geometry["nt_length"],
                        "aa_length": geometry["aa_length"],
                        "exon_starts": geometry["exon_starts"],
                        "exon_ends": geometry["exon_ends"],
                    }
                    output_handle.write(
                        "\t".join(
                            _format_value(output.get(column, ""))
                            for column in reliable_columns
                        )
                        + "\n"
                    )
                    reliable_smorfs_written += 1
        connection.commit()

        if logger is not None:
            logger.write(
                "Stage6: Export reliable smORF genePred annotation."
            )
        reliable_genepred_records = _export_reliable_genepred(
            source_path=orf_genepred,
            temporary_path=genepred_tmp,
            connection=connection,
        )
        if reliable_genepred_records != reliable_smorfs_written:
            raise EvidenceEngineError(
                "Reliable output and genePred record counts differ: "
                f"table={reliable_smorfs_written:,}, "
                f"genePred={reliable_genepred_records:,}."
            )

        totals = {
            "total_families": sum(
                item.total_families for item in chromosome_results
            ),
            "reliable_families": sum(
                item.reliable_families for item in chromosome_results
            ),
            "uncertain_families": sum(
                item.uncertain_families for item in chromosome_results
            ),
            "no_evidence_families": sum(
                item.no_evidence_families for item in chromosome_results
            ),
            "reliable_smorfs": reliable_smorfs_written,
            "invalid_families": sum(
                item.invalid_families for item in chromosome_results
            ),
        }
        evidence_counts: dict[str, int] = defaultdict(int)
        reason_counts: dict[str, int] = defaultdict(int)
        replicated_signal_families = 0
        with gzip.open(master_tmp, "rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                evidence_counts[row.get("best_sample_evidence", "")] += 1
                reason_counts[row.get("reliability_reason", "")] += 1
                try:
                    supported = int(
                        row.get("supported_sample_count", "0") or 0
                    )
                    if supported >= config.reliable_sample:
                        replicated_signal_families += 1
                except ValueError:
                    pass

        with summary_tmp.open("w", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write(f"engine_version\t{ENGINE_VERSION}\n")
            handle.write(f"sample_count\t{sample_count}\n")
            handle.write(f"effective_workers\t{effective_workers}\n")
            for key, value in totals.items():
                handle.write(f"{key}\t{value}\n")
            handle.write(
                "reliable_genepred_records\t"
                f"{reliable_genepred_records}\n"
            )
            handle.write(
                "annotated_controls_excluded_from_reliable_smorf\t"
                f"{annotated_controls_excluded}\n"
            )
            handle.write(
                "reliable_family_definition\t"
                "Medium/High common-body evidence in at least "
                f"{config.reliable_sample} independent samples\n"
            )
            handle.write(
                "reliable_smorf_definition\t"
                "Reliable family with a singleton start or replicated "
                "contiguous alternative-start extension support\n"
            )
            handle.write(
                "uncertain_definition\t"
                "Signal without sufficient independent-sample replication, "
                "localized long-ORF signal, or invalid family geometry\n"
            )
            handle.write(f"evidence_mode\t{config.evidence_mode}\n")
            handle.write(f"group_column\t{config.group_column}\n")
            handle.write(f"reliable_sample\t{config.reliable_sample}\n")
            handle.write(f"short_max_codons\t{config.short_max_codons}\n")
            handle.write(f"long_min_codons\t{config.long_min_codons}\n")
            handle.write(f"window_codons\t{config.window_codons}\n")
            handle.write(f"window_step_codons\t{config.window_step_codons}\n")
            handle.write(f"min_supported_windows\t{config.min_supported_windows}\n")
            handle.write(f"min_window_gap_codons\t{config.min_window_gap_codons}\n")
            handle.write(f"min_signal_span\t{config.min_signal_span}\n")
            handle.write(f"localized_span_max\t{config.localized_span_max}\n")
            handle.write(
                "localized_top_window_fraction\t"
                f"{config.localized_top_window_fraction}\n"
            )
            handle.write(f"boundary_codons\t{config.boundary_codons}\n")
            handle.write(
                f"families_with_replicated_medium_high_support\t"
                f"{replicated_signal_families}\n"
            )
            for label, count in sorted(evidence_counts.items()):
                handle.write(
                    f"best_sample_evidence_{label or 'NA'}\t{count}\n"
                )
            for reason, count in sorted(reason_counts.items()):
                handle.write(
                    f"reliability_reason_{reason or 'NA'}\t{count}\n"
                )
    except BaseException:
        for path in temporary_paths:
            path.unlink(missing_ok=True)
        raise
    finally:
        if connection is not None:
            connection.close()

    for temporary, final in (
        (master_tmp, master_path),
        (reliable_tmp, reliable_path),
        (genepred_tmp, genepred_path),
        (summary_tmp, summary_path),
    ):
        os.replace(temporary, final)

    return EngineResult(
        master_output=str(master_path),
        reliable_output=str(reliable_path),
        reliable_genepred_output=str(genepred_path),
        summary_output=str(summary_path),
        sample_count=sample_count,
        effective_workers=effective_workers,
        **{key: totals[key] for key in (
            "total_families", "reliable_families", "uncertain_families",
            "no_evidence_families", "reliable_smorfs",
        )},
    )


def _build_config(args: object) -> EngineConfig:
    """Build an engine configuration from public and hidden arguments."""
    hidden_fields = (
        "short_max_codons",
        "long_min_codons",
        "window_codons",
        "window_step_codons",
        "min_supported_windows",
        "min_window_gap_codons",
        "min_signal_span",
        "localized_span_max",
        "localized_top_window_fraction",
        "boundary_codons",
    )
    overrides = frozenset(
        field
        for field in hidden_fields
        if getattr(args, field, None) is not None
    )

    def advanced(name: str) -> float | int:
        value = getattr(args, name, None)
        return ADVANCED_DEFAULTS[name] if value is None else value

    return EngineConfig(
        evidence_mode=str(args.evidence_mode),
        group_column=str(args.group_column),
        reliable_sample=int(args.reliable_sample),
        min_rpf_sum=float(args.min_rpf_sum),
        min_rpf_per_codon=float(args.min_rpf_per_codon),
        min_covered_codon=int(args.min_covered_codon),
        min_coverage_ratio=float(args.min_codon_coverage),
        moderate_periodicity=float(args.moderate_periodicity),
        strong_periodicity=float(args.strong_periodicity),
        min_window_rpf=float(args.min_window_rpf),
        min_window_covered=int(
            args.min_window_covered_codon
        ),
        short_max_codons=int(advanced("short_max_codons")),
        long_min_codons=int(advanced("long_min_codons")),
        window_codons=int(advanced("window_codons")),
        window_step_codons=int(advanced("window_step_codons")),
        min_supported_windows=int(
            advanced("min_supported_windows")
        ),
        min_window_gap_codons=int(
            advanced("min_window_gap_codons")
        ),
        min_signal_span=float(advanced("min_signal_span")),
        localized_span_max=float(
            advanced("localized_span_max")
        ),
        localized_top_window_fraction=float(
            advanced("localized_top_window_fraction")
        ),
        boundary_codons=int(advanced("boundary_codons")),
        hidden_overrides=overrides,
        positive_quantile=float(args.positive_quantile),
        positive_min_controls=int(args.positive_min_controls),
        positive_max_controls=int(args.positive_max_controls),
        threads=int(args.threads),
    )



def run_family_evidence_engine(args: object) -> EngineResult:
    """Run the complete bottom-up family evidence engine."""
    tracks = read_density_tracks(args)
    config = _build_config(args)
    output_prefix = Path(args.output).expanduser().resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    work_directory = Path(str(output_prefix) + ".evidence_work")
    log_path = Path(str(output_prefix) + ".evidence_run.log")
    logger = _StageLogger(log_path)
    signature = _run_signature(args, tracks)

    if work_directory.exists():
        shutil.rmtree(work_directory)
    work_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = work_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {"signature": signature, "engine": ENGINE_VERSION},
            indent=2,
        ),
        encoding="utf-8",
    )

    database_path = work_directory / "family_index.sqlite"
    try:
        logger.write("Stage1: Build or reuse the family SQLite index.")
        build_family_database(
            args.family_table,
            args.family_members,
            args.orf_source,
            database_path,
            signature,
            logger,
        )

        config = resolve_length_models(
            database_path,
            config,
            logger,
        )

        logger.write("Stage2: Parse every density track once and build chromosome caches.")
        for number, track in enumerate(tracks, start=1):
            logger.write(
                f"Density track {number}/{len(tracks)}: sample={track.sample}, strand={track.strand}."
            )
            build_density_cache(track, work_directory, signature)

        logger.write("Stage3: Calibrate every sample with annotated ORFs.")
        thresholds, config = calibrate_samples(
            database_path,
            work_directory,
            tracks,
            config,
            logger,
        )

        logger.write("Stage4: Evaluate chromosome family scaffolds.")
        chromosomes = _chromosomes(database_path)
        groups = tuple(sorted({track.group for track in tracks}))
        tasks = [
            ChromosomeTask(
                chromosome=chrom,
                database_path=str(database_path),
                work_directory=str(work_directory),
                output_prefix=str(output_prefix),
                tracks=tuple(tracks),
                thresholds=tuple(thresholds),
                config=config,
                group_names=groups,
            )
            for chrom in chromosomes
        ]
        workers = _effective_workers(config.threads, len(tasks))
        logger.write(
            f"Chromosomes={len(tasks)}, requested_workers={config.threads}, effective_workers={workers}."
        )
        results_by_chrom: dict[str, ChromosomeResult] = {}
        if workers == 1:
            for number, task in enumerate(tasks, start=1):
                result = _process_chromosome(task)
                results_by_chrom[result.chromosome] = result
                progress_print(f"evidence chromosomes: {number}/{len(tasks)}")
        else:
            methods = mp.get_all_start_methods()
            context = mp.get_context("fork") if "fork" in methods else mp.get_context()
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
                futures = {executor.submit(_process_chromosome, task): task.chromosome for task in tasks}
                completed = 0
                for future in as_completed(futures):
                    chrom = futures[future]
                    try:
                        result = future.result()
                    except Exception as error:
                        raise EvidenceEngineError(
                            f"Chromosome task failed: {chrom}: {error}"
                        ) from error
                    results_by_chrom[result.chromosome] = result
                    completed += 1
                    progress_print(f"evidence chromosomes: {completed}/{len(tasks)}")
        ordered_results = [results_by_chrom[chrom] for chrom in chromosomes]

        logger.write("Stage5: Merge focused family and reliable-smORF outputs.")
        result = _merge_outputs(
            output_prefix,
            ordered_results,
            database_path,
            config,
            sample_count=len({track.sample for track in tracks}),
            effective_workers=workers,
            orf_genepred=args.orf_genepred,
            logger=logger,
        )
        logger.write(
            "Completed: total={total:,}, reliable={reliable:,}, uncertain={uncertain:,}, no_evidence={none:,}.".format(
                total=result.total_families,
                reliable=result.reliable_families,
                uncertain=result.uncertain_families,
                none=result.no_evidence_families,
            )
        )
        shutil.rmtree(work_directory, ignore_errors=True)
        return result
    except Exception as error:
        logger.error("run_family_evidence_engine", error)
        warning_print(
            f"smorf_evidence failed. The exact stage and traceback are saved in: {log_path}"
        )
        shutil.rmtree(work_directory, ignore_errors=True)
        raise
