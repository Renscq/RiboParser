#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Run transcript-centric smORF scanning and stream output files.
# Input: Genome FASTA, genePred annotation, and ORF scanning parameters.
# Output: genePred, message table, nucleotide FASTA, and peptide FASTA files.

"""High-level transcript-centric smORF scanning pipeline.

The pipeline performs the following operations:

1. Read and validate genome FASTA records.
2. Read and validate genePred transcript annotation.
3. Reconstruct spliced transcript sequences.
4. Scan complete and optional 3-prime partial ORFs.
5. Classify ORFs relative to the annotated CDS.
6. Mark overlap relationships when requested.
7. Assign stable ORF identifiers in transcript input order.
8. Stream four output formats through atomic temporary files.

Multiprocessing uses ordered ``ProcessPoolExecutor.map`` with an internally
calculated chunk size. It avoids creating one Future object per transcript and
preserves deterministic output ordering. On Linux, ``fork`` is selected
explicitly in automatic mode so the large genome dictionary can be shared
copy-on-write. Other platforms use their available default process context.
"""

from __future__ import annotations

import math
import multiprocessing
import sys
import threading
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
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


def _init_worker(
    genome: dict[str, str],
    config: dict[str, Any],
) -> None:
    """Initialize process-local genome and scanner configuration.

    Args:
        genome: Genome sequence dictionary.
        config: Immutable scanner and post-processing configuration.
    """
    global _WORKER_GENOME
    global _WORKER_CONFIG

    _WORKER_GENOME = genome
    _WORKER_CONFIG = config


def _scan_transcript(
    tx: Transcript,
    genome: dict[str, str],
    config: dict[str, Any],
) -> list[ORFRecord]:
    """Run the complete transcript-local scanning workflow.

    Args:
        tx: Transcript to scan.
        genome: Genome sequence dictionary.
        config: Scanner and post-processing configuration.

    Returns:
        Classified and optionally overlap-filtered ORF records.

    Raises:
        RuntimeError: If transcript reconstruction or scanning fails.
    """
    try:
        scanner = ORFScanner(
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
) -> tuple[int, str, str, list[ORFRecord]]:
    """Scan one transcript in a worker process.

    Args:
        task: Transcript input index and Transcript object.

    Returns:
        Input index, transcript ID, gene ID, and ORF records.

    Raises:
        RuntimeError: If worker initialization is incomplete.
    """
    if _WORKER_GENOME is None or _WORKER_CONFIG is None:
        raise RuntimeError("smORF worker was not initialized.")

    index, transcript = task
    records = _scan_transcript(
        tx=transcript,
        genome=_WORKER_GENOME,
        config=_WORKER_CONFIG,
    )
    return (
        index,
        transcript.transcript_id,
        transcript.gene_id,
        records,
    )


class SmORFPipeline:
    """Run transcript-centric smORF discovery.

    Args:
        genome: Genome FASTA or FASTA.GZ path.
        annotation: Basic genePred or genePredExt annotation path.
        out_prefix: Prefix of output files.
        orf_prefix: Prefix of generated ORF identifiers.
        start_codons: Comma-separated candidate start codons.
        min_aa: Minimum peptide length.
        max_aa: Maximum peptide length.
        scan_strand: ``sense``, ``antisense``, or ``both``.
        kozak_up: Upstream Kozak-context length.
        kozak_down: Downstream Kozak-context length.
        mark_overlap: Mark same-frame and different-frame overlaps.
        remove_discarded: Remove same-frame internal ORFs classified as
            discarded.
        include_stop: Retain the terminal ``*`` in peptide output.
        keep_partial: Retain 3-prime partial ORFs without an in-frame stop.
        allow_ambiguous: Retain ORFs containing ambiguous codons.
        threads: Number of worker processes.
        retain_records: Retain all ORF objects in ``records`` after writing.
            Disable this for command-line runs with millions of ORFs.
        process_start_method: Multiprocessing start method. ``auto`` selects
            ``fork`` on Linux when available and the platform default
            elsewhere.

    Raises:
        ValueError: If a pipeline parameter is invalid.
    """

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
        keep_partial: bool = True,
        allow_ambiguous: bool = False,
        threads: int = 1,
        retain_records: bool = True,
        process_start_method: str = "auto",
    ) -> None:
        """Initialize and validate the scanning pipeline.

        Args:
            genome: Genome FASTA path.
            annotation: genePred annotation path.
            out_prefix: Output file prefix.
            orf_prefix: Generated ORF identifier prefix.
            start_codons: Comma-separated candidate start codons.
            min_aa: Minimum peptide length.
            max_aa: Maximum peptide length.
            scan_strand: Transcript orientation to scan.
            kozak_up: Upstream Kozak-context length.
            kozak_down: Downstream Kozak-context length.
            mark_overlap: Annotate ORF overlaps.
            remove_discarded: Remove discarded same-frame internal ORFs.
            include_stop: Retain terminal stop symbols in peptides.
            keep_partial: Retain 3-prime partial ORFs.
            allow_ambiguous: Retain ORFs with ambiguous codons.
            threads: Worker process count.
            retain_records: Retain all ORF objects in memory.
            process_start_method: Multiprocessing start method.

        Raises:
            ValueError: If any pipeline parameter is invalid.
        """
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

        # Validate all scanner-level parameters once before reading large files.
        ORFScanner(
            start_codons=self.start_codons,
            min_aa=self.min_aa,
            max_aa=self.max_aa,
            scan_strand=self.scan_strand,
            kozak_up=self.kozak_up,
            kozak_down=self.kozak_down,
            include_stop=self.include_stop,
            keep_partial=self.keep_partial,
            allow_ambiguous=self.allow_ambiguous,
        )
        if not self.orf_prefix:
            raise ValueError("orf_prefix must not be empty.")

    def _config(self) -> dict[str, Any]:
        """Return worker-safe scanner configuration.

        Returns:
            Plain dictionary suitable for process initialization.
        """
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
        }

    @staticmethod
    def _progress_interval(total_transcripts: int) -> int:
        """Return an adaptive progress interval.

        Args:
            total_transcripts: Total number of transcripts.

        Returns:
            Number of completed transcripts between progress messages.
        """
        if total_transcripts <= 0:
            return 1
        return max(
            1,
            min(
                1000,
                int(math.ceil(total_transcripts / 100.0)),
            ),
        )

    def _process_context(self) -> multiprocessing.context.BaseContext:
        """Return the configured multiprocessing context.

        Returns:
            Multiprocessing context.

        Raises:
            ValueError: If the requested start method is unavailable.
        """
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
        """Return an internal executor map chunk size.

        Args:
            total_transcripts: Total number of transcripts.

        Returns:
            Bounded positive chunk size.
        """
        target_chunks = max(1, self.threads * 8)
        return max(
            1,
            min(
                64,
                int(math.ceil(total_transcripts / target_chunks)),
            ),
        )

    def _single_process_results(
        self,
        genome: dict[str, str],
        transcripts: list[Transcript],
    ) -> Iterator[tuple[int, str, str, list[ORFRecord]]]:
        """Yield transcript results in single-process mode.

        Args:
            genome: Genome sequence dictionary.
            transcripts: Transcripts in stable input order.

        Yields:
            Input index, transcript ID, gene ID, and ORF records.
        """
        config = self._config()
        for index, transcript in enumerate(transcripts, start=1):
            records = _scan_transcript(
                tx=transcript,
                genome=genome,
                config=config,
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
    ) -> Iterator[tuple[int, str, str, list[ORFRecord]]]:
        """Yield transcript results in deterministic multiprocessing order.

        Args:
            genome: Genome sequence dictionary.
            transcripts: Transcripts in stable input order.

        Yields:
            Input index, transcript ID, gene ID, and ORF records.
        """
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

    def run(self) -> None:
        """Run the complete scanning and output workflow.

        The method is safe to call more than once on the same object. Previous
        in-memory records and counters are cleared before every run.

        Raises:
            ValueError: If no transcript is available.
            RuntimeError: If transcript scanning fails.
        """
        self.records.clear()
        self.transcript_count = 0
        self.orf_count = 0
        self._outputs_written = False

        genome = FastaParser.read_fasta(self.genome_path)
        transcripts = GenePredParser.read_genepred(self.annotation_path)
        if not transcripts:
            raise ValueError("No transcript is available for smORF scanning.")

        self.transcript_count = len(transcripts)
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
                        "transcripts={done:,}/{total:,}, last_gene={gene}, "
                        "last_transcript={transcript}, total_ORFs={orfs:,}.".format(
                            done=completed,
                            total=self.transcript_count,
                            gene=gene_id,
                            transcript=transcript_id,
                            orfs=self.orf_count,
                        )
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
        """Write retained records through the legacy non-streaming API.

        This method is retained for external callers that construct or modify
        ``records`` manually. ``run`` already writes outputs and therefore this
        method becomes a no-op after a successful pipeline run.

        Raises:
            ValueError: If outputs have not been written and no records exist.
        """
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
