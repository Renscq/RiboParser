#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev024
# Function: Orchestrate family-aware smORF evidence analysis.
# Input: Validated CLI arguments.
# Output: Family evidence, reliable smORFs, and summary.

"""Orchestrate family-aware smORF evidence analysis."""

from __future__ import annotations

import csv
import gzip
import json
import multiprocessing as mp
import os
import shutil
import sqlite3
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

from utils.ribo.ArgsParser import progress_print, warning_print

from .calibration import (
    calibrate_samples,
    resolve_length_models,
)
from .competition import _apply_nested_competition
from .config import ADVANCED_DEFAULTS
from .extent import _output_columns
from .input import (
    _chromosomes,
    _effective_workers,
    _run_signature,
    build_density_cache,
    build_family_database,
    read_density_tracks,
)
from .models import (
    ChromosomeResult,
    ChromosomeTask,
    EngineConfig,
    EngineResult,
    EvidenceEngineError,
    _StageLogger,
)
from .output import (
    _export_reliable_genepred,
    _format_value,
    _is_annotated_category,
    _process_chromosome,
)


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
    master_path = Path(str(output_prefix) + ".smorf_evidence.txt")
    reliable_path = Path(str(output_prefix) + ".reliable_smorf.txt")
    genepred_path = Path(str(output_prefix) + ".reliable_smorf.genepred")
    summary_path = Path(str(output_prefix) + ".evidence_summary.txt")
    master_raw_tmp = master_path.with_name(master_path.name + ".raw.tmp")
    master_tmp = master_path.with_name(master_path.name + ".tmp")
    reliable_tmp = reliable_path.with_name(reliable_path.name + ".tmp")
    genepred_tmp = genepred_path.with_name(genepred_path.name + ".tmp")
    summary_tmp = summary_path.with_name(summary_path.name + ".tmp")
    temporary_paths = (
        master_raw_tmp,
        master_tmp,
        reliable_tmp,
        genepred_tmp,
        summary_tmp,
    )
    reliable_columns = (
        "family_id",
        "gene_id",
        "chrom",
        "strand",
        "category",
        "family_type",
        "family_size",
        "structural_primary",
        "provisional_primary",
        "quant_primary",
        "evidence_primary",
        "transcript_id",
        "translation_unit_id",
        "translation_unit_index",
        "translation_unit_count",
        "translation_unit_status",
        "start_codon",
        "stop_codon",
        "nt_length",
        "aa_length",
        "exon_starts",
        "exon_ends",
        "family_translation_status",
        "translated_extent_status",
        "extent_supported_sample_count",
        "extent_supporting_samples",
        "extent_supported_bins",
        "start_site_status",
        "start_interval_orf_ids",
        "selection_policy",
        "selection_reason",
        "canonical_anchor_orf",
        "longest_candidate_orf",
        "prior_override_status",
        "leading_support_sample_count",
        "leading_negative_sample_count",
        "noncanonical_extension_support_sample_count",
        "noncanonical_extension_density_ratio",
        "pooled_extent_rpf_sum",
        "pooled_extent_frame0_ratio",
        "pooled_start_rpf_sum",
        "pooled_end_rpf_sum",
        "pooled_extent_sample_count",
        "nested_competition_status",
        "dominant_parent_orf",
        "competition_overlap_fraction",
        "competition_frame_relation",
        "competition_reason",
        "frame_margin",
        "supported_sample_count",
        "supporting_samples",
        "best_sample",
        "best_sample_evidence",
        "best_rpf_sum",
        "best_coverage_ratio",
        "best_frame0_ratio",
        "reliability_reason",
        "start_site_reason",
        "complete_candidate_count",
        "complete_candidate_orf_ids",
        "candidate_audit_status",
        "candidate_failure_reasons",
        "pooled_frame0_density",
        "pooled_frame1_density",
        "pooled_frame2_density",
        "pooled_dominant_frame",
        "suggested_psite_shift_nt",
        "phase_audit_status",
    )
    reliable_smorfs_written = 0
    annotated_controls_excluded = 0
    reliable_genepred_records = 0
    connection: sqlite3.Connection | None = None
    try:
        with gzip.open(
            master_raw_tmp,
            "wt",
            encoding="utf-8",
            compresslevel=1,
        ) as handle:
            handle.write("\t".join(_output_columns()) + "\n")
        with master_raw_tmp.open("ab") as target:
            for result in chromosome_results:
                with Path(result.part_path).open("rb") as source:
                    shutil.copyfileobj(
                        source,
                        target,
                        length=16 * 1024 * 1024,
                    )

        connection = sqlite3.connect(database_path)
        connection.row_factory = sqlite3.Row
        if logger is not None:
            logger.write(
                "Stage5b: Resolve transcript translation units and nested ORF competition."
            )
        competition_counts, translation_unit_counts = _apply_nested_competition(
            master_raw_tmp,
            master_tmp,
            connection,
            config,
        )
        connection.execute(
            """
            CREATE TEMP TABLE reliable_export (
                orf_id TEXT PRIMARY KEY,
                matched INTEGER NOT NULL DEFAULT 0
            ) WITHOUT ROWID
            """
        )
        with reliable_tmp.open("w", encoding="utf-8", buffering=8 * 1024 * 1024) as output_handle:
            output_handle.write("\t".join(reliable_columns) + "\n")
            with master_tmp.open("r", encoding="utf-8") as input_handle:
                reader = csv.DictReader(input_handle, delimiter="\t")
                for row in reader:
                    if not (
                        row["evidence_status"] == "Reliable"
                        and row["translated_extent_status"] in {"Supported", "Compatible"}
                        and row["quant_primary"]
                        and row["nested_competition_status"]
                        not in {
                            "ShadowedByLongORF",
                            "ShadowedByHighOverlapORF",
                            "ShadowedByCanonicalATG",
                        }
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
                        (row["quant_primary"],),
                    ).fetchone()
                    if geometry is None:
                        raise EvidenceEngineError(
                            "Reliable quantification primary is absent from "
                            f"the geometry index: {row['quant_primary']}"
                        )
                    if _is_annotated_category(geometry["category"]):
                        annotated_controls_excluded += 1
                        continue
                    try:
                        connection.execute(
                            "INSERT INTO reliable_export (orf_id) VALUES (?)",
                            (row["quant_primary"],),
                        )
                    except sqlite3.IntegrityError as error:
                        raise EvidenceEngineError(
                            "One ORF is the quantification primary of multiple "
                            "reliable families: "
                            f"{row['quant_primary']}"
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
                            _format_value(output.get(column, "")) for column in reliable_columns
                        )
                        + "\n"
                    )
                    reliable_smorfs_written += 1
        connection.commit()

        if logger is not None:
            logger.write("Stage6: Export reliable smORF genePred annotation.")
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
            "total_families": sum(item.total_families for item in chromosome_results),
            "reliable_families": sum(item.reliable_families for item in chromosome_results),
            "uncertain_families": sum(item.uncertain_families for item in chromosome_results),
            "no_evidence_families": sum(item.no_evidence_families for item in chromosome_results),
            "reliable_smorfs": reliable_smorfs_written,
            "invalid_families": sum(item.invalid_families for item in chromosome_results),
        }
        evidence_counts: dict[str, int] = defaultdict(int)
        reason_counts: dict[str, int] = defaultdict(int)
        phase_audit_counts: dict[str, int] = defaultdict(int)
        annotated_phase_shift_counts: dict[int, int] = defaultdict(int)
        replicated_signal_families = 0
        with master_tmp.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                evidence_counts[row.get("best_sample_evidence", "")] += 1
                reason_counts[row.get("reliability_reason", "")] += 1
                phase_audit_counts[row.get("phase_audit_status", "")] += 1
                if (
                    _is_annotated_category(row.get("category", ""))
                    and row.get("phase_audit_status") == "AlternativeFrameDominant"
                ):
                    annotated_phase_shift_counts[
                        int(row.get("suggested_psite_shift_nt", "0") or 0)
                    ] += 1
                try:
                    supported = int(row.get("supported_sample_count", "0") or 0)
                    if supported >= config.reliable_sample:
                        replicated_signal_families += 1
                except ValueError:
                    pass

        with summary_tmp.open("w", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write(f"sample_count\t{sample_count}\n")
            handle.write(f"effective_workers\t{effective_workers}\n")
            for key, value in totals.items():
                handle.write(f"{key}\t{value}\n")
            handle.write(f"reliable_genepred_records\t{reliable_genepred_records}\n")
            handle.write(
                f"annotated_controls_excluded_from_reliable_smorf\t{annotated_controls_excluded}\n"
            )
            handle.write(
                "reliable_family_definition\t"
                "Candidate-first distributed extent evidence or Medium/High "
                "common-body evidence in at least "
                f"{config.reliable_sample} independent samples\n"
            )
            handle.write(
                "reliable_smorf_definition\t"
                "Reliable family with a supported or prior-compatible "
                "translated extent, a quant_primary, and no dominant ORF "
                "within the same transcript translation unit\n"
            )
            handle.write(
                "uncertain_definition\t"
                "Signal without sufficient independent-sample replication, "
                "localized long-ORF signal, or invalid family geometry\n"
            )
            handle.write(f"evidence_mode\t{config.evidence_mode}\n")
            handle.write(f"group_column\t{config.group_column}\n")
            handle.write(f"reliable_sample\t{config.reliable_sample}\n")
            handle.write(f"estimated_short_max_codons\t{config.estimated_short_max_codons}\n")
            handle.write(f"estimated_long_min_codons\t{config.estimated_long_min_codons}\n")
            handle.write(f"short_max_codons\t{config.short_max_codons}\n")
            handle.write(f"long_min_codons\t{config.long_min_codons}\n")
            handle.write(f"window_codons\t{config.window_codons}\n")
            handle.write(f"window_step_codons\t{config.window_step_codons}\n")
            handle.write(f"min_supported_windows\t{config.min_supported_windows}\n")
            handle.write(f"min_window_gap_codons\t{config.min_window_gap_codons}\n")
            handle.write(f"min_signal_span\t{config.min_signal_span}\n")
            handle.write(f"localized_span_max\t{config.localized_span_max}\n")
            handle.write(f"localized_top_window_fraction\t{config.localized_top_window_fraction}\n")
            handle.write(f"boundary_codons\t{config.boundary_codons}\n")
            handle.write(f"extent_bins\t{config.extent_bins}\n")
            handle.write(f"min_extent_bins\t{config.min_extent_bins}\n")
            handle.write(f"start_resolution_codons\t{config.start_resolution_codons}\n")
            handle.write(f"leading_window_codons\t{config.leading_window_codons}\n")
            handle.write(f"min_exclusion_codons\t{config.min_exclusion_codons}\n")
            handle.write(
                f"silent_extension_density_ratio\t{config.silent_extension_density_ratio}\n"
            )
            handle.write(
                "noncanonical_override_density_ratio\t"
                f"{config.noncanonical_override_density_ratio}\n"
            )
            handle.write(
                "noncanonical_override_min_coverage_ratio\t"
                f"{config.noncanonical_override_min_coverage_ratio}\n"
            )
            handle.write(
                "noncanonical_override_min_frame_margin\t"
                f"{config.noncanonical_override_min_frame_margin}\n"
            )
            handle.write(
                f"noncanonical_min_exclusive_codons\t{config.noncanonical_min_exclusive_codons}\n"
            )
            handle.write(f"nested_min_frame_margin\t{config.nested_min_frame_margin}\n")
            handle.write(f"nested_min_phase_rpf\t{config.nested_min_phase_rpf}\n")
            handle.write(f"high_overlap_fraction\t{config.high_overlap_fraction}\n")
            handle.write(f"overlap_min_frame_margin\t{config.overlap_min_frame_margin}\n")
            handle.write(f"overlap_min_phase_rpf\t{config.overlap_min_phase_rpf}\n")
            for label, count in sorted(competition_counts.items()):
                handle.write(f"nested_competition_{label}\t{count}\n")
            handle.write(f"translation_unit_total\t{sum(translation_unit_counts.values())}\n")
            for label, count in sorted(translation_unit_counts.items()):
                handle.write(f"translation_unit_{label}\t{count}\n")
            handle.write(
                f"families_with_replicated_medium_high_support\t{replicated_signal_families}\n"
            )
            for label, count in sorted(evidence_counts.items()):
                handle.write(f"best_sample_evidence_{label or 'NA'}\t{count}\n")
            for label, count in sorted(phase_audit_counts.items()):
                handle.write(f"phase_audit_{label or 'NA'}\t{count}\n")
            for shift, count in sorted(annotated_phase_shift_counts.items()):
                handle.write(f"annotated_control_suggested_psite_shift_{shift:+d}_nt\t{count}\n")
            for reason, count in sorted(reason_counts.items()):
                handle.write(f"reliability_reason_{reason or 'NA'}\t{count}\n")
        if logger is not None and annotated_phase_shift_counts:
            dominant_shift, dominant_count = max(
                annotated_phase_shift_counts.items(),
                key=lambda item: item[1],
            )
            if dominant_count >= config.positive_min_controls:
                warning_print(
                    "Annotated ORF phase audit detected a systematic "
                    f"alternative frame: suggested P-site coordinate shift "
                    f"{dominant_shift:+d} nt in {dominant_count:,} controls. "
                    "Verify the upstream P-site offset before interpreting "
                    "candidate reading frames."
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
    master_raw_tmp.unlink(missing_ok=True)

    return EngineResult(
        master_output=str(master_path),
        reliable_output=str(reliable_path),
        reliable_genepred_output=str(genepred_path),
        summary_output=str(summary_path),
        sample_count=sample_count,
        effective_workers=effective_workers,
        **{
            key: totals[key]
            for key in (
                "total_families",
                "reliable_families",
                "uncertain_families",
                "no_evidence_families",
                "reliable_smorfs",
            )
        },
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
        "extent_bins",
        "min_extent_bins",
        "start_resolution_codons",
        "min_frame_margin",
        "leading_window_codons",
        "min_exclusion_codons",
        "silent_extension_density_ratio",
        "noncanonical_min_exclusive_codons",
        "noncanonical_override_density_ratio",
        "noncanonical_override_min_coverage_ratio",
        "noncanonical_override_min_frame_margin",
        "nested_min_frame_margin",
        "nested_min_phase_rpf",
        "high_overlap_fraction",
        "overlap_min_frame_margin",
        "overlap_min_phase_rpf",
    )
    overrides = frozenset(
        field for field in hidden_fields if getattr(args, field, None) is not None
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
        min_window_covered=int(args.min_window_covered_codon),
        short_max_codons=int(advanced("short_max_codons")),
        long_min_codons=int(advanced("long_min_codons")),
        window_codons=int(advanced("window_codons")),
        window_step_codons=int(advanced("window_step_codons")),
        min_supported_windows=int(advanced("min_supported_windows")),
        min_window_gap_codons=int(advanced("min_window_gap_codons")),
        min_signal_span=float(advanced("min_signal_span")),
        localized_span_max=float(advanced("localized_span_max")),
        localized_top_window_fraction=float(advanced("localized_top_window_fraction")),
        boundary_codons=int(advanced("boundary_codons")),
        extent_bins=int(advanced("extent_bins")),
        min_extent_bins=int(advanced("min_extent_bins")),
        start_resolution_codons=int(advanced("start_resolution_codons")),
        min_frame_margin=float(advanced("min_frame_margin")),
        leading_window_codons=int(advanced("leading_window_codons")),
        min_exclusion_codons=int(advanced("min_exclusion_codons")),
        silent_extension_density_ratio=float(advanced("silent_extension_density_ratio")),
        noncanonical_min_exclusive_codons=int(advanced("noncanonical_min_exclusive_codons")),
        noncanonical_override_density_ratio=float(advanced("noncanonical_override_density_ratio")),
        noncanonical_override_min_coverage_ratio=float(
            advanced("noncanonical_override_min_coverage_ratio")
        ),
        noncanonical_override_min_frame_margin=float(
            advanced("noncanonical_override_min_frame_margin")
        ),
        nested_min_frame_margin=float(advanced("nested_min_frame_margin")),
        nested_min_phase_rpf=float(advanced("nested_min_phase_rpf")),
        high_overlap_fraction=float(advanced("high_overlap_fraction")),
        overlap_min_frame_margin=float(advanced("overlap_min_frame_margin")),
        overlap_min_phase_rpf=float(advanced("overlap_min_phase_rpf")),
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
    logger = _StageLogger()
    signature = _run_signature(args, tracks)

    if work_directory.exists():
        shutil.rmtree(work_directory)
    work_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = work_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {"signature": signature},
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
                futures = {
                    executor.submit(_process_chromosome, task): task.chromosome for task in tasks
                }
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
            f"Completed: total={result.total_families:,}, reliable={result.reliable_families:,}, uncertain={result.uncertain_families:,}, no_evidence={result.no_evidence_families:,}."
        )
        shutil.rmtree(work_directory, ignore_errors=True)
        return result
    except Exception as error:
        logger.error("run_family_evidence_engine", error)
        shutil.rmtree(work_directory, ignore_errors=True)
        raise
