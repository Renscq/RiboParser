#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-24
# Version: 0.2.8.24-dev.002
# Function: Run high-performance transcript-centric smORF scanning and streaming output.
# Input: Genome FASTA, genePred annotation, and ORF scanning parameters.
# Output: genePred, metadata, nucleotide FASTA, and peptide FASTA files.

"""High-performance transcript-centric smORF scanning pipeline.

The command-line workflow uses two execution paths:

1. Single-process or compatibility mode scans each transcript once and streams
   records directly to the final atomic writer.
2. Large multi-process jobs use a count-and-shard workflow. A lightweight first
   pass determines deterministic ORF identifier ranges and output workload.
   Workers then scan balanced contiguous transcript partitions and write local
   shards directly. The parent process concatenates those shards without
   transferring sequence-heavy ORF objects through multiprocessing pipes.

The shard mode preserves transcript input order, stable ORF identifiers, output
schemas, and atomic final-file replacement.
"""

from __future__ import annotations

import math
import multiprocessing
import shutil
import sys
import tempfile
import threading
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from utils.ribo.ArgsParser import message_print, progress_print, result_print

from .smorf_classifier import ORFClassifier
from .smorf_coordinate import CoordinateMapper
from .smorf_fasta import FastaParser
from .smorf_genepred import GenePredParser
from .smorf_models import ORFRecord, Transcript
from .smorf_overlap import ORFOverlapMarker
from .smorf_scanner import ORFScanner
from .smorf_writer import (
    FastaWriter,
    GenePredWriter,
    MessageWriter,
    ORFOutputWriter,
)

_WORKER_GENOME: dict[str, str] | None = None
_WORKER_CONFIG: dict[str, Any] | None = None
_WORKER_SCANNER: ORFScanner | None = None


TranscriptResult = tuple[int, str, str, list[ORFRecord]]
TranscriptSummary = tuple[int, str, str, int, int, int]
PartEntry = tuple[int, Transcript, int, int]
PartTask = tuple[int, str, list[PartEntry]]
PartResult = tuple[int, str, int, str, str]


def _build_scanner(config: dict[str, Any]) -> ORFScanner:
    """Build one scanner from immutable pipeline configuration."""
    return ORFScanner(
        start_codons=config["start_codons"],
        min_aa=config["min_aa"],
        max_aa=config["max_aa"],
        scan_strand=config["scan_strand"],
        kozak_up=config["kozak_up"],
        kozak_down=config["kozak_down"],
        include_stop=config["include_stop"],
        keep_partial=config["keep_partial"],
        allow_ambiguous=config["allow_ambiguous"],
    )


def _init_worker(
    genome: dict[str, str],
    config: dict[str, Any],
) -> None:
    """Initialize process-local immutable resources once."""
    global _WORKER_GENOME
    global _WORKER_CONFIG
    global _WORKER_SCANNER

    _WORKER_GENOME = genome
    _WORKER_CONFIG = config
    _WORKER_SCANNER = _build_scanner(config)


def _require_worker_resources() -> tuple[
    dict[str, str],
    dict[str, Any],
    ORFScanner,
]:
    """Return initialized worker resources."""
    if (
        _WORKER_GENOME is None
        or _WORKER_CONFIG is None
        or _WORKER_SCANNER is None
    ):
        raise RuntimeError("smORF worker was not initialized.")
    return _WORKER_GENOME, _WORKER_CONFIG, _WORKER_SCANNER


def _scan_transcript(
    tx: Transcript,
    genome: dict[str, str],
    config: dict[str, Any],
    scanner: ORFScanner,
) -> list[ORFRecord]:
    """Run transcript reconstruction, scanning, and classification."""
    try:
        CoordinateMapper.build_transcript_sequence(tx, genome)
        records = scanner.scan_transcript(tx)
        ORFClassifier.classify(tx, records)

        if config["mark_overlap"]:
            ORFOverlapMarker.mark(records)
        if config["remove_discarded"]:
            records = [
                record
                for record in records
                if record.priority != "discarded"
            ]
        return records
    except Exception as error:
        raise RuntimeError(
            "smORF scanning failed for "
            f"gene={tx.gene_id}, transcript={tx.transcript_id}, "
            f"chrom={tx.chrom}, strand={tx.strand}: {error}"
        ) from error


def _scan_transcript_worker(
    task: tuple[int, Transcript],
) -> TranscriptResult:
    """Scan one transcript in an initialized worker process."""
    genome, config, scanner = _require_worker_resources()
    index, transcript = task
    records = _scan_transcript(
        tx=transcript,
        genome=genome,
        config=config,
        scanner=scanner,
    )
    return (
        index,
        transcript.transcript_id,
        transcript.gene_id,
        records,
    )


def _summarize_transcript_worker(
    task: tuple[int, Transcript],
) -> TranscriptSummary:
    """Count retained ORFs without constructing sequence-heavy records."""
    genome, config, scanner = _require_worker_resources()
    index, transcript = task

    try:
        CoordinateMapper.build_transcript_sequence(transcript, genome)
        orf_count, nucleotide_length, peptide_length = (
            scanner.summarize_transcript(
                transcript,
                remove_discarded=config["remove_discarded"],
            )
        )
    except Exception as error:
        raise RuntimeError(
            "smORF summary failed for "
            f"gene={transcript.gene_id}, "
            f"transcript={transcript.transcript_id}, "
            f"chrom={transcript.chrom}, strand={transcript.strand}: {error}"
        ) from error

    return (
        index,
        transcript.transcript_id,
        transcript.gene_id,
        orf_count,
        nucleotide_length,
        peptide_length,
    )


def _write_part_worker(task: PartTask) -> PartResult:
    """Scan one balanced transcript partition and write local shards."""
    genome, config, scanner = _require_worker_resources()
    part_index, part_prefix, entries = task
    written_count = 0
    last_transcript = "NA"
    last_gene = "NA"

    with ORFOutputWriter(part_prefix) as output_writer:
        for _input_index, transcript, first_orf_index, expected_count in entries:
            records = _scan_transcript(
                tx=transcript,
                genome=genome,
                config=config,
                scanner=scanner,
            )
            if len(records) != expected_count:
                raise RuntimeError(
                    "ORF count changed between summary and output passes for "
                    f"transcript={transcript.transcript_id}: "
                    f"expected={expected_count}, observed={len(records)}."
                )

            for local_index, record in enumerate(records):
                record.orf_id = (
                    f"{config['orf_prefix']}"
                    f"{first_orf_index + local_index:08d}"
                )

            output_writer.write_records(records)
            written_count += len(records)
            last_transcript = transcript.transcript_id
            last_gene = transcript.gene_id

    return (
        part_index,
        part_prefix,
        written_count,
        last_transcript,
        last_gene,
    )


class SmORFPipeline:
    """Run transcript-centric smORF discovery."""

    def __init__(
        self,
        genome: str,
        annotation: str,
        out_prefix: str = "ORF",
        orf_prefix: str = "ORF",
        start_codons: str = "ATG",
        min_aa: int = 8,
        max_aa: int = 10000,
        scan_strand: str = "sense",
        kozak_up: int = 6,
        kozak_down: int = 6,
        mark_overlap: bool = False,
        remove_discarded: bool = False,
        include_stop: bool = False,
        keep_partial: bool = False,
        allow_ambiguous: bool = False,
        threads: int = 1,
        retain_records: bool = True,
        process_start_method: str = "auto",
    ) -> None:
        """Initialize and validate the scanning pipeline."""
        self.genome_path = str(genome)
        self.annotation_path = str(annotation)
        self.out_prefix = str(out_prefix)
        self.orf_prefix = str(orf_prefix)
        self.start_codons = [
            codon.strip().upper().replace("U", "T")
            for codon in str(start_codons).split(",")
            if codon.strip()
        ]
        self.min_aa = int(min_aa)
        self.max_aa = int(max_aa)
        self.scan_strand = str(scan_strand)
        self.kozak_up = int(kozak_up)
        self.kozak_down = int(kozak_down)
        self.mark_overlap = bool(mark_overlap)
        self.remove_discarded = bool(remove_discarded)
        self.include_stop = bool(include_stop)
        self.keep_partial = bool(keep_partial)
        self.allow_ambiguous = bool(allow_ambiguous)
        self.threads = max(1, int(threads))
        self.retain_records = bool(retain_records)
        self.process_start_method = str(process_start_method)
        self.records: list[ORFRecord] = []
        self.transcript_count = 0
        self.orf_count = 0
        self._outputs_written = False

        _build_scanner(self._config())
        if not self.orf_prefix:
            raise ValueError("orf_prefix must not be empty.")

    def _config(self) -> dict[str, Any]:
        """Return worker-safe immutable scanner configuration."""
        return {
            "start_codons": tuple(self.start_codons),
            "min_aa": self.min_aa,
            "max_aa": self.max_aa,
            "scan_strand": self.scan_strand,
            "kozak_up": self.kozak_up,
            "kozak_down": self.kozak_down,
            "include_stop": self.include_stop,
            "keep_partial": self.keep_partial,
            "allow_ambiguous": self.allow_ambiguous,
            "mark_overlap": self.mark_overlap,
            "remove_discarded": self.remove_discarded,
            "orf_prefix": self.orf_prefix,
        }

    @staticmethod
    def _progress_interval(total_items: int) -> int:
        """Return an adaptive progress interval."""
        if total_items <= 0:
            return 1
        return max(
            1,
            min(
                1000,
                int(math.ceil(total_items / 100.0)),
            ),
        )

    def _process_context(self) -> multiprocessing.context.BaseContext:
        """Return the configured multiprocessing context."""
        available_methods = multiprocessing.get_all_start_methods()
        if self.process_start_method == "auto":
            if (
                sys.platform.startswith("linux")
                and "fork" in available_methods
                and threading.active_count() == 1
            ):
                method = "fork"
            elif "forkserver" in available_methods:
                method = "forkserver"
            else:
                method = multiprocessing.get_start_method()
        else:
            method = self.process_start_method

        if method not in available_methods:
            raise ValueError(
                f"Unavailable multiprocessing start method: {method}. "
                f"Available methods: {', '.join(available_methods)}"
            )
        return multiprocessing.get_context(method)

    def _chunksize(self, total_transcripts: int) -> int:
        """Return a bounded process-map chunk size."""
        target_chunks = max(1, self.threads * 32)
        return max(
            1,
            min(
                32,
                int(math.ceil(total_transcripts / target_chunks)),
            ),
        )

    def _single_process_results(
        self,
        genome: dict[str, str],
        transcripts: list[Transcript],
    ) -> Iterator[TranscriptResult]:
        """Yield results while reusing one scanner."""
        config = self._config()
        scanner = _build_scanner(config)
        for index, transcript in enumerate(transcripts, start=1):
            records = _scan_transcript(
                tx=transcript,
                genome=genome,
                config=config,
                scanner=scanner,
            )
            yield (
                index,
                transcript.transcript_id,
                transcript.gene_id,
                records,
            )

    def _multi_process_results(
        self,
        genome: dict[str, str],
        transcripts: list[Transcript],
    ) -> Iterator[TranscriptResult]:
        """Yield deterministic compatibility-mode process results."""
        config = self._config()
        context = self._process_context()
        chunksize = self._chunksize(len(transcripts))
        message_print(
            "Multiprocessing: workers={workers}, start_method={method}, "
            "chunksize={chunksize}.".format(
                workers=self.threads,
                method=context.get_start_method(),
                chunksize=chunksize,
            )
        )

        tasks = enumerate(transcripts, start=1)
        with ProcessPoolExecutor(
            max_workers=self.threads,
            mp_context=context,
            initializer=_init_worker,
            initargs=(genome, config),
        ) as executor:
            yield from executor.map(
                _scan_transcript_worker,
                tasks,
                chunksize=chunksize,
            )

    def _use_sharded_mode(self, total_transcripts: int) -> bool:
        """Return whether sequence-heavy IPC should be bypassed."""
        return (
            self.threads > 1
            and not self.retain_records
            and total_transcripts >= max(512, self.threads * 32)
        )

    def _build_part_tasks(
        self,
        transcripts: list[Transcript],
        summaries: list[TranscriptSummary],
        part_directory: Path,
    ) -> tuple[list[PartTask], int]:
        """Build contiguous output partitions balanced by estimated work."""
        if len(transcripts) != len(summaries):
            raise RuntimeError("Transcript summary count is inconsistent.")

        entries_with_weight: list[tuple[PartEntry, int]] = []
        next_orf_index = 1
        for expected_index, (transcript, summary) in enumerate(
            zip(transcripts, summaries),
            start=1,
        ):
            (
                input_index,
                transcript_id,
                gene_id,
                orf_count,
                nucleotide_length,
                peptide_length,
            ) = summary
            if input_index != expected_index:
                raise RuntimeError("Transcript summaries are out of order.")
            if (
                transcript.transcript_id != transcript_id
                or transcript.gene_id != gene_id
            ):
                raise RuntimeError(
                    "Transcript summary identity is inconsistent at index "
                    f"{expected_index}."
                )

            entry: PartEntry = (
                input_index,
                transcript,
                next_orf_index,
                orf_count,
            )
            weight = (
                transcript.transcript_length()
                + nucleotide_length
                + peptide_length
                + orf_count * 256
            )
            entries_with_weight.append((entry, max(1, weight)))
            next_orf_index += orf_count

        total_orfs = next_orf_index - 1
        desired_parts = min(
            len(entries_with_weight),
            max(self.threads, self.threads * 4),
        )
        total_weight = sum(weight for _entry, weight in entries_with_weight)
        remaining_weight = total_weight
        remaining_parts = desired_parts
        current_entries: list[PartEntry] = []
        current_weight = 0
        partitions: list[list[PartEntry]] = []

        for position, (entry, weight) in enumerate(entries_with_weight):
            current_entries.append(entry)
            current_weight += weight
            remaining_entries = len(entries_with_weight) - position - 1
            target_weight = remaining_weight / max(1, remaining_parts)

            should_close = (
                remaining_parts > 1
                and remaining_entries >= remaining_parts - 1
                and current_weight >= target_weight
            )
            if should_close:
                partitions.append(current_entries)
                remaining_weight -= current_weight
                remaining_parts -= 1
                current_entries = []
                current_weight = 0

        if current_entries:
            partitions.append(current_entries)

        tasks: list[PartTask] = []
        for part_index, entries in enumerate(partitions, start=1):
            part_prefix = part_directory / f"part{part_index:05d}"
            tasks.append((part_index, str(part_prefix), entries))
        return tasks, total_orfs

    def _run_sharded_multiprocessing(
        self,
        genome: dict[str, str],
        transcripts: list[Transcript],
    ) -> None:
        """Run balanced worker-side output without ORF-record IPC."""
        config = self._config()
        context = self._process_context()
        chunksize = self._chunksize(len(transcripts))
        output_prefix = Path(self.out_prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        part_directory = Path(
            tempfile.mkdtemp(
                prefix=f".{output_prefix.name}.smorf_parts.",
                dir=str(output_prefix.parent),
            )
        )

        message_print(
            "Fast multiprocessing: workers={workers}, start_method={method}, "
            "summary_chunksize={chunksize}.".format(
                workers=self.threads,
                method=context.get_start_method(),
                chunksize=chunksize,
            )
        )

        try:
            with ProcessPoolExecutor(
                max_workers=self.threads,
                mp_context=context,
                initializer=_init_worker,
                initargs=(genome, config),
            ) as executor:
                summary_iterator = executor.map(
                    _summarize_transcript_worker,
                    enumerate(transcripts, start=1),
                    chunksize=chunksize,
                )
                summaries: list[TranscriptSummary] = []
                summary_interval = self._progress_interval(
                    self.transcript_count
                )
                for completed, summary in enumerate(
                    summary_iterator,
                    start=1,
                ):
                    summaries.append(summary)
                    if (
                        completed == self.transcript_count
                        or completed % summary_interval == 0
                    ):
                        progress_print(
                            "indexed_transcripts={done:,}/{total:,}.".format(
                                done=completed,
                                total=self.transcript_count,
                            )
                        )

                part_tasks, expected_orf_count = self._build_part_tasks(
                    transcripts=transcripts,
                    summaries=summaries,
                    part_directory=part_directory,
                )
                message_print(
                    "Balanced output partitions: {parts:,}.".format(
                        parts=len(part_tasks)
                    )
                )

                part_prefixes: list[str] = []
                observed_orf_count = 0
                for completed, result in enumerate(
                    executor.map(
                        _write_part_worker,
                        part_tasks,
                        chunksize=1,
                    ),
                    start=1,
                ):
                    (
                        part_index,
                        part_prefix,
                        written_count,
                        last_transcript,
                        last_gene,
                    ) = result
                    if part_index != completed:
                        raise RuntimeError("Output partitions are out of order.")
                    part_prefixes.append(part_prefix)
                    observed_orf_count += written_count
                    progress_print(
                        "output_parts={done:,}/{total:,}, "
                        "last_gene={gene}, last_transcript={transcript}, "
                        "total_ORFs={orfs:,}.".format(
                            done=completed,
                            total=len(part_tasks),
                            gene=last_gene,
                            transcript=last_transcript,
                            orfs=observed_orf_count,
                        )
                    )

            if observed_orf_count != expected_orf_count:
                raise RuntimeError(
                    "Final ORF count differs from the summary pass: "
                    f"expected={expected_orf_count}, "
                    f"observed={observed_orf_count}."
                )

            ORFOutputWriter.merge_parts(
                output_prefix=self.out_prefix,
                part_prefixes=part_prefixes,
            )
            self.orf_count = observed_orf_count
        finally:
            shutil.rmtree(part_directory, ignore_errors=True)

    def _run_streaming(
        self,
        genome: dict[str, str],
        transcripts: list[Transcript],
    ) -> None:
        """Run the standard one-pass streaming workflow."""
        progress_interval = self._progress_interval(self.transcript_count)
        if self.threads == 1:
            result_iterator = self._single_process_results(
                genome=genome,
                transcripts=transcripts,
            )
        else:
            result_iterator = self._multi_process_results(
                genome=genome,
                transcripts=transcripts,
            )

        next_orf_index = 1
        with ORFOutputWriter(self.out_prefix) as output_writer:
            for completed, (
                _input_index,
                transcript_id,
                gene_id,
                transcript_records,
            ) in enumerate(result_iterator, start=1):
                for record in transcript_records:
                    record.orf_id = (
                        f"{self.orf_prefix}{next_orf_index:08d}"
                    )
                    next_orf_index += 1

                output_writer.write_records(transcript_records)
                self.orf_count += len(transcript_records)
                if self.retain_records:
                    self.records.extend(transcript_records)

                if (
                    completed == 1
                    or completed == self.transcript_count
                    or completed % progress_interval == 0
                ):
                    progress_print(
                        "transcripts={done:,}/{total:,}, "
                        "last_gene={gene}, last_transcript={transcript}, "
                        "total_ORFs={orfs:,}.".format(
                            done=completed,
                            total=self.transcript_count,
                            gene=gene_id,
                            transcript=transcript_id,
                            orfs=self.orf_count,
                        )
                    )

    def run(self) -> None:
        """Run scanning and write all output formats."""
        self.records.clear()
        self.transcript_count = 0
        self.orf_count = 0
        self._outputs_written = False

        genome = FastaParser.read_fasta(self.genome_path)
        transcripts = GenePredParser.read_genepred(self.annotation_path)
        if not transcripts:
            raise ValueError("No transcript is available for smORF scanning.")

        self.transcript_count = len(transcripts)
        if self._use_sharded_mode(self.transcript_count):
            self._run_sharded_multiprocessing(
                genome=genome,
                transcripts=transcripts,
            )
        else:
            self._run_streaming(
                genome=genome,
                transcripts=transcripts,
            )

        self._outputs_written = True
        result_print(
            [
                ("Transcripts", f"{self.transcript_count:,}"),
                ("ORFs", f"{self.orf_count:,}"),
                ("Retained in memory", f"{len(self.records):,}"),
            ]
        )

    def write_outputs(self) -> None:
        """Write retained records through the legacy API."""
        if self._outputs_written:
            return
        if not self.records:
            raise ValueError("No ORF record is available for output.")

        GenePredWriter.write(
            f"{self.out_prefix}.genePred",
            self.records,
        )
        MessageWriter.write(
            f"{self.out_prefix}.message.txt",
            self.records,
        )
        FastaWriter.write_nt(
            f"{self.out_prefix}.nt.fa",
            self.records,
        )
        FastaWriter.write_pep(
            f"{self.out_prefix}.pep.fa",
            self.records,
        )
        self._outputs_written = True
