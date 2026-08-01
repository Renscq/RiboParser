#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Read and index family-evidence inputs.
# Input: Family tables, ORF sources, and density design.
# Output: SQLite family index and cached density tracks.

"""Read and index family-evidence inputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import sqlite3
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence, TextIO

import numpy as np

from utils.ribo.ArgsParser import progress_print
from utils.smorf.density import (
    ChromDensity,
    infer_density_format,
    iter_density_chromosomes,
    prepare_density,
)

from .config import (
    FAMILY_BATCH_SIZE,
    MAX_WORKERS,
    MEMORY_PER_WORKER_BYTES,
    SCHEMA_VERSION,
    SQL_BATCH_SIZE,
    STOP_CODONS,
)
from .models import (
    DensityTrack,
    FamilyGeometry,
    InvalidFamily,
    RepresentativeGeometry,
    _StageLogger,
)


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
        if key != "threads" and isinstance(value, (str, int, float, bool, type(None)))
    }
    payload = json.dumps(
        {
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
            raise ValueError("Density list requires a tab-delimited header.")
        header_lookup = {str(name).strip().lower(): str(name) for name in original_header}
        path_name = next(
            (header_lookup[name] for name in ("path", "file", "density") if name in header_lookup),
            None,
        )
        sample_name = next(
            (header_lookup[name] for name in ("sample", "name") if name in header_lookup),
            None,
        )
        if sample_name is None or "strand" not in header_lookup:
            raise ValueError(
                "Density list requires sample/name, strand, and path/file/density columns."
            )
        group_key = str(group_column).strip().lower()
        if group_key not in {"sample", "name"} and group_key not in header_lookup:
            raise ValueError(f"Density-list group column was not found: {group_column}")
        if path_name is None:
            raise ValueError("Density list requires one path column: path, file, or density.")

        strand_name = header_lookup["strand"]
        group_name = sample_name if group_key in {"sample", "name"} else header_lookup[group_key]
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
                    "family_id",
                    "gene_id",
                    "chrom",
                    "strand",
                    "category",
                    "family_type",
                    "family_size",
                    "orf_id",
                    "transcript_id",
                    "start_codon",
                    "stop_codon",
                    "completeness",
                    "nt_length",
                    "aa_length",
                    "exon_starts",
                    "exon_ends",
                ),
                "family table",
            )
            family_rows: list[tuple[Any, ...]] = []
            rep_rows: list[tuple[Any, ...]] = []
            geometry_rows: list[tuple[Any, ...]] = []
            for number, row in enumerate(reader, start=1):
                family_rows.append(
                    (
                        row["family_id"],
                        row["gene_id"],
                        row["chrom"],
                        row["strand"],
                        row["category"],
                        row["family_type"],
                        int(row["family_size"]),
                        row["orf_id"],
                    )
                )
                rep_rows.append((row["family_id"], row["orf_id"], 1))
                geometry_rows.append(
                    (
                        row["orf_id"],
                        row["transcript_id"],
                        row["category"],
                        row["start_codon"],
                        row["stop_codon"],
                        row["completeness"],
                        row["chrom"],
                        row["strand"],
                        int(row["nt_length"]),
                        int(row["aa_length"]),
                        row["exon_starts"],
                        row["exon_ends"],
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
                    family_rows.clear()
                    rep_rows.clear()
                    geometry_rows.clear()
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
                    connection.commit()
                    rows.clear()
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
                    "orf_id",
                    "transcript_id",
                    "category",
                    "start_codon",
                    "stop_codon",
                    "completeness",
                    "chrom",
                    "strand",
                    "nt_length",
                    "aa_length",
                    "exon_starts",
                    "exon_ends",
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
                        row["orf_id"],
                        row["transcript_id"],
                        row["category"],
                        row["start_codon"],
                        row["stop_codon"],
                        row["completeness"],
                        row["chrom"],
                        row["strand"],
                        int(row["nt_length"]),
                        int(row["aa_length"]),
                        row["exon_starts"],
                        row["exon_ends"],
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
                    connection.commit()
                    batch.clear()
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
            raise ValueError(f"Missing source coordinates for {missing:,} family representatives.")
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
    return (
        work_directory
        / "density"
        / _safe_name(f"{track.sample}|{track.strand}|{Path(track.path).resolve()}")
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
            progress_print(f"density cache {track.sample}/{track.strand}: chromosome {number}")
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
    tail = scaffold[-len(representative) :]
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
                raise ValueError(f"phase-incompatible representative {record['orf_id']}")
            if not _is_suffix(
                scaffold["starts"],
                scaffold["ends"],
                record["starts"],
                record["ends"],
                base["strand"],
            ):
                raise ValueError(f"non-suffix representative {record['orf_id']}")
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
                    stop_codon=record["stop_codon"],
                    completeness=record["completeness"],
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
        values = [row[0] for row in connection.execute("SELECT DISTINCT chrom FROM families")]
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
