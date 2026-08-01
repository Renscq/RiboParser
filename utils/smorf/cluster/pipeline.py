#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev008
# Function: Orchestrate parallel smORF family clustering.
# Input: Scanner table, genePred annotation, and cluster settings.
# Output: Family tables, removed candidates, and summary.

"""Orchestrate parallel smORF family clustering."""

from __future__ import annotations

import bisect
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import BinaryIO

from .annotation import AnnotationIndex
from .config import _BUFFER_BYTES, _FLUSH_BYTES, _MIN_PARALLEL_BYTES, MEMBER_COLUMNS, STOP_CODONS
from .family import _FamilyMixin
from .filtering import _FilteringMixin
from .input import _ColumnPlan
from .kozak import KozakModel
from .models import (
    ClusterSummary,
    _ClusterConfig,
    _FastKozakScorer,
    _KozakSpec,
    _RangeResult,
    _RangeTask,
)
from .output import _ShardWriter

_WORKER_ANNOTATION: AnnotationIndex | None = None
_WORKER_PLAN: _ColumnPlan | None = None
_WORKER_CONFIG: _ClusterConfig | None = None
_WORKER_INPUT: str | None = None


def _worker_initialize(
    annotation_path: str,
    input_path: str,
    header: tuple[bytes, ...],
    config: _ClusterConfig,
) -> None:
    """Initialize one long-lived clustering worker."""
    global _WORKER_ANNOTATION, _WORKER_PLAN, _WORKER_CONFIG, _WORKER_INPUT
    _WORKER_ANNOTATION = AnnotationIndex.from_genepred(annotation_path)
    _WORKER_PLAN = _ColumnPlan.from_header(header)
    _WORKER_CONFIG = config
    _WORKER_INPUT = input_path


def _process_range_worker(task: _RangeTask) -> _RangeResult:
    """Process one independent gene-safe byte range."""
    if (
        _WORKER_ANNOTATION is None
        or _WORKER_PLAN is None
        or _WORKER_CONFIG is None
        or _WORKER_INPUT is None
    ):
        raise RuntimeError("smorf_cluster worker was not initialized.")
    writer = _ShardWriter(Path(task.temp_directory), task.shard_id)
    engine = _ClusterEngine(
        annotation=_WORKER_ANNOTATION,
        plan=_WORKER_PLAN,
        config=_WORKER_CONFIG,
    )
    with writer:
        summary = engine.process_range(
            input_path=_WORKER_INPUT,
            start=task.start,
            end=task.end,
            writer=writer,
        )
    return _RangeResult(
        shard_id=task.shard_id,
        paths=(
            str(writer.paths["family"]),
            str(writer.paths["members"]),
            str(writer.paths["removed"]),
        ),
        summary=summary,
    )


class _ClusterEngine(_FilteringMixin, _FamilyMixin):
    """Coordinate filtering and family construction for one input range."""

    def __init__(
        self,
        annotation: AnnotationIndex,
        plan: _ColumnPlan,
        config: _ClusterConfig,
    ) -> None:
        self.annotation = annotation
        self.plan = plan
        self.config = config
        self.start_rank = dict(config.start_rank)
        self.kozak_scorer = (
            _FastKozakScorer(config.kozak_spec) if config.kozak_spec is not None else None
        )

        lookup = plan.lookup
        self.i_orf_id = lookup[b"orf_id"]
        self.i_gene_id = lookup[b"gene_id"]
        self.i_transcript_id = lookup[b"transcript_id"]
        self.i_chrom = lookup[b"chrom"]
        self.i_strand = lookup[b"strand"]
        self.i_source_strand = lookup[b"source_strand"]
        self.i_category = lookup[b"category"]
        self.i_priority = lookup[b"priority"]
        self.i_start_codon = lookup[b"start_codon"]
        self.i_stop_codon = lookup[b"stop_codon"]
        self.i_nt_length = lookup[b"nt_length"]
        self.i_aa_length = lookup[b"aa_length"]
        self.i_completeness = lookup[b"completeness"]
        self.i_exon_count = lookup[b"exon_count"]
        self.i_exon_starts = lookup[b"exon_starts"]
        self.i_exon_ends = lookup[b"exon_ends"]
        self.i_ambiguous = lookup.get(b"ambiguous_codon_count", -1)
        self.i_kozak_seq = lookup.get(b"kozak_seq", -1)
        self.i_kozak_start = lookup.get(b"kozak_start_index", -1)
        self.i_kozak_name = lookup[b"kozak_pwm_name"]
        self.i_kozak_score = lookup[b"kozak_pwm_score"]
        self.i_kozak_level = lookup[b"kozak_pwm_level"]
        self.i_kozak_valid = lookup[b"kozak_valid_ratio"]
        self.i_structure = lookup[b"structure_status"]
        self.i_filter_status = lookup[b"filter_status"]
        self.i_filter_reason = lookup[b"filter_reason"]
        self.i_matched_annotated = lookup[b"matched_annotated_orf_id"]
        self.i_annotated_relation = lookup[b"annotated_overlap_relation"]
        self.i_annotated_overlap_nt = lookup[b"annotated_overlap_nt"]
        self.i_annotated_overlap_codon = lookup[b"annotated_overlap_codon"]

    def process_range(
        self,
        input_path: str | Path,
        start: int,
        end: int,
        writer: _ShardWriter,
    ) -> ClusterSummary:
        """Filter and cluster one gene-safe input range."""
        summary = ClusterSummary()
        transcript_groups = self._iter_transcript_groups(
            input_path=input_path,
            start=start,
            end=end,
            summary=summary,
            writer=writer,
        )
        local_family_number = 0
        for _gene_id, candidates in self._iter_gene_groups(transcript_groups):
            summary.genes_processed += 1
            candidates = self._remove_annotated_overlaps(
                candidates,
                summary,
                writer,
            )
            representatives = self._deduplicate_exact(candidates, summary)
            families = self._cluster_representatives(representatives)
            for family in families:
                local_family_number += 1
                self._write_family(
                    family,
                    local_family_number,
                    summary,
                    writer,
                )
        return summary


class SmORFCluster:
    """Filter and cluster smORFs with deterministic high-performance I/O."""

    def __init__(
        self,
        annotation: str | Path,
        keep_start_codons: str = "",
        min_aa: int = 8,
        max_aa: int = 10000,
        keep_categories: str = (
            "uORF,dORF,lncORF,iORF,same_frame_iORF,emORF,"
            "overlap_uORF,overlap_dORF,other_ORF,"
            "annotated_ORF,annotated_mORF"
        ),
        remove_categories: str = "antisense_ORF",
        require_sense: bool = True,
        max_ambiguous_codons: int = 0,
        kozak_model: KozakModel | None = None,
        min_kozak_score: float = 0.0,
        start_priority: str = "ATG",
        threads: int = 0,
    ) -> None:
        self.annotation_path = str(Path(annotation).resolve())
        self.annotation = AnnotationIndex.from_genepred(annotation)
        self.threads = int(threads)
        if self.threads < 0:
            raise ValueError("threads must be >= 0.")
        if min_aa < 1:
            raise ValueError("min_aa must be >= 1.")
        if max_aa < min_aa:
            raise ValueError("max_aa must be >= min_aa.")
        if max_ambiguous_codons < 0:
            raise ValueError("max_ambiguous_codons must be >= 0.")
        if not 0 <= min_kozak_score <= 1:
            raise ValueError("min_kozak_score must be in [0, 1].")
        if min_kozak_score > 0 and kozak_model is None:
            raise ValueError("A Kozak model is required when min_kozak_score > 0.")

        keep_start_set = self._parse_bytes_set(keep_start_codons, normalize_codon=True)
        keep_category_set = self._parse_bytes_set(keep_categories)
        remove_category_set = self._parse_bytes_set(remove_categories)
        conflict = keep_category_set.intersection(remove_category_set)
        if conflict:
            raise ValueError(
                "Category cannot be both retained and removed: "
                + ", ".join(item.decode("utf-8", errors="replace") for item in sorted(conflict))
            )
        start_order = self._parse_start_order(start_priority)
        self.config = _ClusterConfig(
            keep_start_codons=frozenset(keep_start_set),
            min_aa=int(min_aa),
            max_aa=int(max_aa),
            keep_categories=frozenset(keep_category_set),
            remove_categories=frozenset(remove_category_set),
            require_sense=bool(require_sense),
            max_ambiguous_codons=int(max_ambiguous_codons),
            min_kozak_score=float(min_kozak_score),
            start_rank=tuple((codon, rank) for rank, codon in enumerate(start_order)),
            kozak_spec=_KozakSpec.from_model(kozak_model),
        )

    @staticmethod
    def _parse_bytes_set(
        value: str,
        normalize_codon: bool = False,
    ) -> set[bytes]:
        output: set[bytes] = set()
        for item in str(value).split(","):
            text = item.strip()
            if not text:
                continue
            if normalize_codon:
                text = text.upper().replace("U", "T")
            output.add(text.encode("utf-8"))
        return output

    @staticmethod
    def _parse_start_order(value: str) -> list[bytes]:
        codons: list[bytes] = []
        seen: set[bytes] = set()
        for item in str(value).split(","):
            codon_text = item.strip().upper().replace("U", "T")
            if not codon_text:
                continue
            if len(codon_text) != 3 or set(codon_text).difference("ACGT"):
                raise ValueError(f"Invalid start-priority codon: {codon_text}")
            codon = codon_text.encode("ascii")
            if codon in STOP_CODONS:
                raise ValueError(f"Stop codon cannot be a start-priority codon: {codon_text}")
            if codon not in seen:
                seen.add(codon)
                codons.append(codon)
        if not codons:
            raise ValueError("Start-codon priority cannot be empty.")
        return codons

    @staticmethod
    def read_header(path: str | Path) -> tuple[bytes, ...]:
        """Read the first non-empty table header line."""
        with Path(path).open("rb", buffering=_BUFFER_BYTES) as handle:
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError(f"Empty ORF message table: {path}")
                stripped = line.strip()
                if stripped:
                    return tuple(stripped.split(b"\t"))

    @staticmethod
    def _data_start(path: str | Path) -> int:
        with Path(path).open("rb", buffering=_BUFFER_BYTES) as handle:
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError(f"Empty ORF message table: {path}")
                if line.strip():
                    return handle.tell()

    @staticmethod
    def _extract_field(line: bytes, index: int) -> bytes:
        """Extract one tab field without splitting the complete row."""
        start = 0
        for current in range(index):
            position = line.find(b"\t", start)
            if position < 0:
                return b""
            start = position + 1
        end = line.find(b"\t", start)
        if end < 0:
            end = len(line)
        return line[start:end].strip()

    def _effective_workers(self, file_size: int) -> int:
        if self.threads > 0:
            requested = self.threads
        else:
            slurm = os.environ.get("SLURM_CPUS_PER_TASK")
            requested = int(slurm) if slurm and slurm.isdigit() else (os.cpu_count() or 1)
        requested = max(1, requested)
        if file_size < _MIN_PARALLEL_BYTES:
            return 1
        return requested

    def _discover_ranges(
        self,
        input_path: str | Path,
        plan: _ColumnPlan,
        workers: int,
    ) -> list[tuple[int, int]]:
        """Find approximately balanced cuts that never split one gene."""
        input_file = Path(input_path)
        file_size = input_file.stat().st_size
        data_start = self._data_start(input_file)
        if workers <= 1 or file_size <= data_start:
            return [(data_start, file_size)]

        transcript_index = plan.lookup[b"transcript_id"]
        safe_offsets = [data_start]
        current_transcript: bytes | None = None
        current_order: int | None = None
        unknown_transcript = False

        with input_file.open("rb", buffering=_BUFFER_BYTES) as handle:
            handle.seek(data_start)
            while True:
                line_start = handle.tell()
                line = handle.readline()
                if not line:
                    break
                transcript_id = self._extract_field(
                    line.rstrip(b"\r\n"),
                    transcript_index,
                )
                if transcript_id != current_transcript:
                    if (
                        current_transcript is not None
                        and current_order is not None
                        and self.annotation.is_safe_cut(current_order)
                    ):
                        safe_offsets.append(line_start)
                    current_transcript = transcript_id
                    meta = self.annotation.transcripts.get(transcript_id)
                    if meta is None:
                        unknown_transcript = True
                        current_order = None
                    else:
                        current_order = meta.annotation_order

        if unknown_transcript:
            return [(data_start, file_size)]
        if safe_offsets[-1] != file_size:
            safe_offsets.append(file_size)
        safe_offsets = sorted(set(safe_offsets))
        if len(safe_offsets) <= 2:
            return [(data_start, file_size)]

        selected = [data_start]
        span = file_size - data_start
        for part in range(1, workers):
            target = data_start + span * part // workers
            position = bisect.bisect_left(safe_offsets, target)
            candidates = []
            if position < len(safe_offsets):
                candidates.append(safe_offsets[position])
            if position > 0:
                candidates.append(safe_offsets[position - 1])
            boundary = min(candidates, key=lambda value: abs(value - target))
            if boundary > selected[-1] and boundary < file_size:
                selected.append(boundary)
        selected.append(file_size)
        return [(left, right) for left, right in zip(selected, selected[1:]) if right > left]

    @staticmethod
    def _copy_file(source: Path, target: BinaryIO) -> None:
        with source.open("rb", buffering=_BUFFER_BYTES) as handle:
            shutil.copyfileobj(handle, target, length=_BUFFER_BYTES)

    @staticmethod
    def _format_family_id(number: int) -> bytes:
        return f"SMORF_FAM{number:09d}".encode("ascii")

    def _merge_family_shard(
        self,
        source: Path,
        target: BinaryIO,
        family_id_index: int,
        offset: int,
    ) -> None:
        with source.open("rb", buffering=_BUFFER_BYTES) as handle:
            buffer = bytearray()
            for line in handle:
                fields = line.rstrip(b"\r\n").split(b"\t")
                local_number = int(fields[family_id_index])
                fields[family_id_index] = self._format_family_id(offset + local_number)
                buffer.extend(b"\t".join(fields) + b"\n")
                if len(buffer) >= _FLUSH_BYTES:
                    target.write(buffer)
                    buffer.clear()
            if buffer:
                target.write(buffer)

    def _merge_member_shard(
        self,
        source: Path,
        target: BinaryIO,
        offset: int,
    ) -> None:
        with source.open("rb", buffering=_BUFFER_BYTES) as handle:
            buffer = bytearray()
            for line in handle:
                fields = line.rstrip(b"\r\n").split(b"\t")
                fields[0] = self._format_family_id(offset + int(fields[0]))
                buffer.extend(b"\t".join(fields) + b"\n")
                if len(buffer) >= _FLUSH_BYTES:
                    target.write(buffer)
                    buffer.clear()
            if buffer:
                target.write(buffer)

    @staticmethod
    def _write_summary(path: Path, summary: ClusterSummary) -> None:
        with path.open("wb", buffering=_BUFFER_BYTES) as handle:
            handle.write(b"section\titem\tcount\n")
            overall = {
                "input_orfs": summary.input_orfs,
                "basic_passed": summary.basic_passed,
                "basic_removed": summary.basic_removed,
                "annotated_overlap_removed": (summary.annotated_overlap_removed),
                "lnc_morf_removed": summary.lnc_morf_removed,
                "post_annotation_retained": (
                    summary.basic_passed - summary.annotated_overlap_removed
                ),
                "exact_duplicates_collapsed": summary.exact_duplicates_collapsed,
                "alt_starts_collapsed": summary.alt_starts_collapsed,
                "total_collapsed": (
                    summary.exact_duplicates_collapsed + summary.alt_starts_collapsed
                ),
                "family_count": summary.family_count,
                "singleton_families": summary.singleton_families,
                "multi_member_families": summary.multi_member_families,
                "genes_processed": summary.genes_processed,
                "effective_workers": summary.effective_workers,
            }
            for item, count in overall.items():
                handle.write(f"overall\t{item}\t{count}\n".encode())
            for reason, count in sorted(summary.removal_reasons.items()):
                handle.write(f"removal_reason\t{reason}\t{count}\n".encode())
            for category, count in sorted(summary.primary_categories.items()):
                handle.write(f"primary_category\t{category}\t{count}\n".encode())

    def cluster_file(
        self,
        input_path: str | Path,
        output_prefix: str | Path,
    ) -> ClusterSummary:
        """Cluster a scanner table with optional gene-safe multiprocessing."""
        input_file = Path(input_path).resolve()
        header = self.read_header(input_file)
        plan = _ColumnPlan.from_header(header)
        file_size = input_file.stat().st_size
        requested_workers = self._effective_workers(file_size)
        ranges = self._discover_ranges(input_file, plan, requested_workers)
        effective_workers = len(ranges)

        output_prefix_path = Path(output_prefix)
        output_prefix_path.parent.mkdir(parents=True, exist_ok=True)
        final_paths = {
            "family": Path(f"{output_prefix_path}.family.message.txt"),
            "members": Path(f"{output_prefix_path}.family.members.txt"),
            "removed": Path(f"{output_prefix_path}.removed.message.txt"),
            "summary": Path(f"{output_prefix_path}.cluster_summary.txt"),
        }
        temp_final = {
            name: path.with_name(path.name + ".tmp") for name, path in final_paths.items()
        }

        with tempfile.TemporaryDirectory(
            prefix=f".{output_prefix_path.name}.cluster.",
            dir=str(output_prefix_path.parent),
        ) as temp_directory:
            tasks = [
                _RangeTask(
                    shard_id=index,
                    start=start,
                    end=end,
                    temp_directory=temp_directory,
                )
                for index, (start, end) in enumerate(ranges)
            ]
            if effective_workers == 1:
                writer = _ShardWriter(Path(temp_directory), 0)
                engine = _ClusterEngine(
                    annotation=self.annotation,
                    plan=plan,
                    config=self.config,
                )
                with writer:
                    shard_summary = engine.process_range(
                        input_path=input_file,
                        start=ranges[0][0],
                        end=ranges[0][1],
                        writer=writer,
                    )
                results = [
                    _RangeResult(
                        shard_id=0,
                        paths=(
                            str(writer.paths["family"]),
                            str(writer.paths["members"]),
                            str(writer.paths["removed"]),
                        ),
                        summary=shard_summary,
                    )
                ]
            else:
                with ProcessPoolExecutor(
                    max_workers=effective_workers,
                    initializer=_worker_initialize,
                    initargs=(
                        self.annotation_path,
                        str(input_file),
                        header,
                        self.config,
                    ),
                ) as executor:
                    results = list(executor.map(_process_range_worker, tasks))

            results.sort(key=lambda item: item.shard_id)
            summary = ClusterSummary(effective_workers=effective_workers)
            family_id_index = plan.family_lookup[b"family_id"]
            family_offset = 0
            try:
                with (
                    temp_final["family"].open("wb", buffering=_BUFFER_BYTES) as family_handle,
                    temp_final["members"].open("wb", buffering=_BUFFER_BYTES) as member_handle,
                    temp_final["removed"].open("wb", buffering=_BUFFER_BYTES) as removed_handle,
                ):
                    family_handle.write(b"\t".join(plan.family_header) + b"\n")
                    member_handle.write(b"\t".join(MEMBER_COLUMNS) + b"\n")
                    removed_handle.write(b"\t".join(plan.filter_header) + b"\n")
                    for result in results:
                        family_path, member_path, removed_path = map(Path, result.paths)
                        self._merge_family_shard(
                            family_path,
                            family_handle,
                            family_id_index,
                            family_offset,
                        )
                        self._merge_member_shard(
                            member_path,
                            member_handle,
                            family_offset,
                        )
                        self._copy_file(removed_path, removed_handle)
                        family_offset += result.summary.family_count
                        summary.merge(result.summary)
                summary.effective_workers = effective_workers
                self._write_summary(temp_final["summary"], summary)
                for name, final_path in final_paths.items():
                    os.replace(temp_final[name], final_path)
            except Exception:
                for path in temp_final.values():
                    try:
                        path.unlink()
                    except FileNotFoundError:
                        pass
                raise

        return summary
