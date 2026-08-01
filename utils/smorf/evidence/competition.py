#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Resolve translation units and nested-ORF competition.
# Input: Reliable family evidence rows.
# Output: Competition-adjusted evidence rows.

"""Resolve translation units and nested-ORF competition."""

from __future__ import annotations

import csv
import gzip
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from .config import SQL_BATCH_SIZE
from .extent import _output_columns
from .geometry import (
    _block_overlap_length,
    _blocks_contain,
    _geometry_blocks,
    _same_overlap_frame,
    _shorter_overlap_fraction,
)
from .models import EngineConfig
from .output import _format_value


def _assign_translation_units(
    parsed: Sequence[tuple[sqlite3.Row, tuple[tuple[int, int], ...]]],
) -> dict[str, tuple[str, int, int, str]]:
    """Assign overlap-connected ORFs to independent transcript units."""
    if not parsed:
        return {}
    parents = list(range(len(parsed)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first: int, second: int) -> None:
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parents[second_root] = first_root

    interval_order = sorted(
        range(len(parsed)),
        key=lambda index: (
            min(start for start, _ in parsed[index][1]),
            max(end for _, end in parsed[index][1]),
        ),
    )
    active: list[int] = []
    for index in interval_order:
        blocks = parsed[index][1]
        left = min(start for start, _ in blocks)
        active = [other for other in active if max(end for _, end in parsed[other][1]) > left]
        for other in active:
            if _block_overlap_length(blocks, parsed[other][1]) > 0:
                union(index, other)
        active.append(index)

    components: dict[int, list[int]] = defaultdict(list)
    for index in range(len(parsed)):
        components[find(index)].append(index)
    first_row = parsed[0][0]
    strand = str(first_row["strand"])
    ordered_components = sorted(
        components.values(),
        key=lambda indices: (
            min(start for index in indices for start, _ in parsed[index][1])
            if strand == "+"
            else -max(end for index in indices for _, end in parsed[index][1])
        ),
    )
    unit_count = len(ordered_components)
    gene_id = str(first_row["gene_id"])
    transcript_id = str(first_row["transcript_id"])
    assignments: dict[str, tuple[str, int, int, str]] = {}
    for unit_index, indices in enumerate(ordered_components, start=1):
        unit_id = f"{gene_id}|{transcript_id}|TU{unit_index:03d}"
        unit_status = (
            "OverlappingCompetitionUnit" if len(indices) > 1 else "IndependentNonOverlappingUnit"
        )
        for index in indices:
            assignments[str(parsed[index][0]["family_id"])] = (
                unit_id,
                unit_index,
                unit_count,
                unit_status,
            )
    return assignments


def _apply_nested_competition(
    raw_master_path: Path,
    final_master_path: Path,
    connection: sqlite3.Connection,
    config: EngineConfig,
) -> tuple[dict[str, int], dict[str, int]]:
    """Assign transcript units and suppress explainable nested ORFs."""
    connection.execute("DROP TABLE IF EXISTS temp.candidate_calls")
    connection.execute(
        """
        CREATE TEMP TABLE candidate_calls (
            family_id TEXT PRIMARY KEY,
            orf_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            strand TEXT NOT NULL,
            category TEXT NOT NULL,
            evidence_status TEXT NOT NULL,
            frame_margin REAL NOT NULL,
            phase_rpf REAL NOT NULL,
            extent_sample_count INTEGER NOT NULL,
            extent_status TEXT NOT NULL,
            start_site_status TEXT NOT NULL,
            prior_override_status TEXT NOT NULL,
            leading_support_count INTEGER NOT NULL,
            noncanonical_support_count INTEGER NOT NULL,
            competition_status TEXT NOT NULL DEFAULT 'Independent',
            dominant_parent_orf TEXT NOT NULL DEFAULT '',
            overlap_fraction REAL NOT NULL DEFAULT 0.0,
            frame_relation TEXT NOT NULL DEFAULT 'NotEvaluated',
            competition_reason TEXT NOT NULL DEFAULT ''
        ) WITHOUT ROWID
        """
    )
    batch: list[tuple[Any, ...]] = []
    with gzip.open(raw_master_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            provisional_primary = row.get("provisional_primary") or row.get("quant_primary") or ""
            if not provisional_primary:
                continue
            if int(row.get("complete_candidate_count", "0") or 0) < 1 and not row.get(
                "quant_primary"
            ):
                continue
            batch.append(
                (
                    row["family_id"],
                    provisional_primary,
                    row["gene_id"],
                    row["chrom"],
                    row["strand"],
                    row["category"],
                    row.get("evidence_status", "Uncertain"),
                    float(row.get("frame_margin", "0") or 0.0),
                    float(
                        row.get("pooled_frame0_density", "0") or row.get("best_rpf_sum", "0") or 0.0
                    ),
                    int(
                        row.get("pooled_extent_sample_count", "0")
                        or row.get("extent_supported_sample_count", "0")
                        or 0
                    ),
                    row.get("translated_extent_status", ""),
                    row.get("start_site_status", ""),
                    row.get("prior_override_status", ""),
                    int(row.get("leading_support_sample_count", "0") or 0),
                    int(
                        row.get(
                            "noncanonical_extension_support_sample_count",
                            "0",
                        )
                        or 0
                    ),
                )
            )
            if len(batch) >= SQL_BATCH_SIZE:
                connection.executemany(
                    "INSERT INTO candidate_calls "
                    "(family_id,orf_id,gene_id,chrom,strand,category,"
                    "evidence_status,frame_margin,phase_rpf,"
                    "extent_sample_count,extent_status,start_site_status,"
                    "prior_override_status,leading_support_count,"
                    "noncanonical_support_count) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                batch.clear()
        if batch:
            connection.executemany(
                "INSERT INTO candidate_calls "
                "(family_id,orf_id,gene_id,chrom,strand,category,"
                "evidence_status,frame_margin,phase_rpf,"
                "extent_sample_count,extent_status,start_site_status,"
                "prior_override_status,leading_support_count,"
                "noncanonical_support_count) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
    connection.execute(
        "CREATE INDEX candidate_calls_group_idx ON candidate_calls(gene_id, chrom, strand)"
    )

    query = """
        SELECT c.*, g.transcript_id, g.start_codon, g.nt_length,
               g.exon_starts, g.exon_ends
        FROM candidate_calls c
        JOIN geometry g ON g.orf_id=c.orf_id
        ORDER BY c.gene_id, g.transcript_id, c.chrom, c.strand,
                 g.nt_length DESC, c.orf_id
    """
    updates: list[tuple[str, str, float, str, str, str]] = []
    group_rows: list[sqlite3.Row] = []
    group_key: tuple[str, str, str, str] | None = None
    unit_assignments: dict[str, tuple[str, int, int, str]] = {}

    def resolve_group(rows: Sequence[sqlite3.Row]) -> None:
        parsed: list[tuple[sqlite3.Row, tuple[tuple[int, int], ...]]] = [
            (row, _geometry_blocks(row)) for row in rows
        ]
        current_units = _assign_translation_units(parsed)
        unit_assignments.update(current_units)
        suppressed_parent_families: set[str] = set()
        for child_index, (child, child_blocks) in enumerate(parsed):
            competitors: list[
                tuple[
                    sqlite3.Row,
                    tuple[tuple[int, int], ...],
                    float,
                    bool,
                ]
            ] = []
            for parent, parent_blocks in parsed[:child_index]:
                if (
                    current_units[str(parent["family_id"])][0]
                    != current_units[str(child["family_id"])][0]
                ):
                    continue
                if str(parent["family_id"]) in suppressed_parent_families:
                    continue
                if str(parent["evidence_status"]) != "Reliable":
                    continue
                if int(parent["nt_length"]) <= int(child["nt_length"]):
                    continue
                contained = _blocks_contain(parent_blocks, child_blocks)
                overlap_fraction = _shorter_overlap_fraction(
                    parent_blocks,
                    child_blocks,
                )
                if contained or overlap_fraction >= config.high_overlap_fraction:
                    competitors.append(
                        (
                            parent,
                            parent_blocks,
                            overlap_fraction,
                            contained,
                        )
                    )
            if not competitors:
                continue
            parent, parent_blocks, overlap_fraction, contained = max(
                competitors,
                key=lambda value: (
                    value[2],
                    str(value[0]["start_codon"]).upper() == "ATG",
                    int(value[0]["nt_length"]),
                    str(value[0]["orf_id"]),
                ),
            )
            same_frame = _same_overlap_frame(
                parent_blocks,
                child_blocks,
                str(child["strand"]),
            )
            frame_relation = "SameFrame" if same_frame else "DifferentFrame"
            parent_noncanonical_supported = (
                str(parent["start_codon"]).upper() != "ATG"
                and str(parent["extent_status"]) == "Supported"
                and int(parent["leading_support_count"]) >= config.reliable_sample
                and int(parent["noncanonical_support_count"]) >= config.reliable_sample
                and str(parent["prior_override_status"])
                == "OverriddenByStrongNoncanonicalExtension"
            )
            canonical_child_over_parent = (
                same_frame
                and str(child["start_codon"]).upper() == "ATG"
                and str(parent["start_codon"]).upper() != "ATG"
                and not parent_noncanonical_supported
            )
            if canonical_child_over_parent:
                updates.append(
                    (
                        "ShadowedByCanonicalATG",
                        str(child["orf_id"]),
                        overlap_fraction,
                        frame_relation,
                        "canonical_ATG_preferred_without_replicated_"
                        "noncanonical_exclusive_translation",
                        str(parent["family_id"]),
                    )
                )
                suppressed_parent_families.add(str(parent["family_id"]))
                continue
            strong_overlap = overlap_fraction >= config.high_overlap_fraction
            required_margin = (
                config.overlap_min_frame_margin
                if strong_overlap
                else config.nested_min_frame_margin
            )
            required_rpf = (
                config.overlap_min_phase_rpf if strong_overlap else config.nested_min_phase_rpf
            )
            if (
                str(parent["start_codon"]).upper() == "ATG"
                and str(child["start_codon"]).upper() != "ATG"
            ):
                required_margin += 0.05
                required_rpf *= 1.25
            required_rpf = max(
                required_rpf,
                float(parent["phase_rpf"]) * 0.25,
            )
            independent_phase = (
                str(child["extent_status"]) == "Supported"
                and float(child["frame_margin"]) >= required_margin
                and float(child["phase_rpf"]) >= required_rpf
                and int(child["extent_sample_count"]) >= config.reliable_sample
            )
            start_site_resolved = str(child["start_site_status"]) in {
                "Exact",
                "UniqueCandidate",
            }
            noncanonical_start_resolved = str(child["start_codon"]).upper() == "ATG" or (
                str(child["prior_override_status"]) == "OverriddenByStrongNoncanonicalExtension"
                and int(child["leading_support_count"]) >= config.reliable_sample
                and int(child["noncanonical_support_count"]) >= config.reliable_sample
            )
            independent_start = start_site_resolved and noncanonical_start_resolved
            if not same_frame and independent_phase and independent_start:
                updates.append(
                    (
                        "IndependentOverlappingORF",
                        str(parent["orf_id"]),
                        overlap_fraction,
                        frame_relation,
                        "independent_frame_advantage_and_replicated_translation_support",
                        str(child["family_id"]),
                    )
                )
            else:
                status = "ShadowedByLongORF" if contained else "ShadowedByHighOverlapORF"
                reason = (
                    "same_frame_overlap_explained_by_longer_ORF"
                    if same_frame
                    else "insufficient_independent_start_frame_or_"
                    "abundance_evidence_in_high_overlap_region"
                )
                updates.append(
                    (
                        status,
                        str(parent["orf_id"]),
                        overlap_fraction,
                        frame_relation,
                        reason,
                        str(child["family_id"]),
                    )
                )

    for row in connection.execute(query):
        current_key = (
            str(row["gene_id"]),
            str(row["transcript_id"]),
            str(row["chrom"]),
            str(row["strand"]),
        )
        if group_key is None:
            group_key = current_key
        elif current_key != group_key:
            resolve_group(group_rows)
            group_rows = []
            group_key = current_key
        group_rows.append(row)
    if group_rows:
        resolve_group(group_rows)
    if updates:
        connection.executemany(
            "UPDATE candidate_calls SET competition_status=?, "
            "dominant_parent_orf=?, overlap_fraction=?, frame_relation=?, "
            "competition_reason=? WHERE family_id=?",
            updates,
        )

    competition = {
        str(row["family_id"]): (
            str(row["competition_status"]),
            str(row["dominant_parent_orf"]),
            float(row["overlap_fraction"]),
            str(row["frame_relation"]),
            str(row["competition_reason"]),
        )
        for row in connection.execute(
            "SELECT family_id, competition_status, dominant_parent_orf, "
            "overlap_fraction, frame_relation, competition_reason "
            "FROM candidate_calls"
        )
    }
    counts: dict[str, int] = defaultdict(int)
    unit_counts: dict[str, int] = defaultdict(int)
    for unit_id, _, _, unit_status in set(unit_assignments.values()):
        unit_counts[unit_status] += 1
    with (
        gzip.open(raw_master_path, "rt", encoding="utf-8") as source,
        final_master_path.open(
            "w",
            encoding="utf-8",
            buffering=8 * 1024 * 1024,
            newline="",
        ) as target,
    ):
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(
            target,
            fieldnames=_output_columns(),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in reader:
            status = competition.get(row["family_id"])
            if status is not None:
                row["nested_competition_status"] = status[0]
                row["dominant_parent_orf"] = status[1]
                row["competition_overlap_fraction"] = status[2]
                row["competition_frame_relation"] = status[3]
                row["competition_reason"] = status[4]
                counts[status[0]] += 1
            unit = unit_assignments.get(row["family_id"])
            if unit is not None:
                row["translation_unit_id"] = unit[0]
                row["translation_unit_index"] = unit[1]
                row["translation_unit_count"] = unit[2]
                row["translation_unit_status"] = unit[3]
            writer.writerow(
                {column: _format_value(row.get(column, "")) for column in _output_columns()}
            )
    return dict(counts), dict(unit_counts)
