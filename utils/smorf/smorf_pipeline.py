# Author: Rensc
# date: 2026-05-21

"""
Main smORF pipeline.

This pipeline connects all functional modules:
1. Read genome FASTA.
2. Read genePred annotation.
3. Reconstruct transcript sequences.
4. Scan ORFs.
5. Classify ORFs.
6. Mark ORF overlaps.
7. Write output files.

Parallel mode uses multiprocessing instead of threading because ORF scanning
is CPU-intensive and Python threads are limited by the GIL.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed

from .smorf_fasta import FastaParser
from .smorf_genepred import GenePredParser
from .smorf_coordinate import CoordinateMapper
from .smorf_scanner import ORFScanner
from .smorf_classifier import ORFClassifier
from .smorf_overlap import ORFOverlapMarker
from .smorf_writer import GenePredWriter, MessageWriter, FastaWriter


_WORKER_GENOME = None
_WORKER_CONFIG = None


def _init_worker(genome, config):
    """
    Initialize worker-level global objects.

    This avoids sending the genome dictionary to every single transcript task.
    """

    global _WORKER_GENOME
    global _WORKER_CONFIG

    _WORKER_GENOME = genome
    _WORKER_CONFIG = config


def _scan_transcript_worker(task):
    """
    Worker function for scanning one transcript.

    Parameters
    ----------
    task : tuple
        Tuple of transcript index and Transcript object.

    Returns
    -------
    tuple
        Transcript index, transcript ID, gene ID, and ORF records.
    """

    idx, tx = task

    scanner = ORFScanner(
        start_codons=_WORKER_CONFIG["start_codons"],
        min_aa=_WORKER_CONFIG["min_aa"],
        max_aa=_WORKER_CONFIG["max_aa"],
        scan_strand=_WORKER_CONFIG["scan_strand"],
        kozak_up=_WORKER_CONFIG["kozak_up"],
        kozak_down=_WORKER_CONFIG["kozak_down"],
        include_stop=_WORKER_CONFIG["include_stop"],
    )

    # Reconstruct spliced transcript sequence before ORF scanning.
    CoordinateMapper.build_transcript_sequence(tx, _WORKER_GENOME)

    # Scan candidate ORFs from the transcript sequence.
    tx_records = scanner.scan_transcript(tx)

    # Assign ORF category labels.
    ORFClassifier.classify(tx, tx_records)

    # Mark nested or overlapping ORFs if requested.
    if _WORKER_CONFIG["mark_overlap"]:
        ORFOverlapMarker.mark(tx_records)

    # Remove same-frame internal ORFs if requested.
    if _WORKER_CONFIG["remove_discarded"]:
        tx_records = [x for x in tx_records if x.priority != "discarded"]

    return idx, tx.transcript_id, tx.gene_id, tx_records


class SmORFPipeline:
    """
    High-level pipeline for transcript-centric smORF detection.
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
        threads: int = 1,
    ):
        """
        Initialize smORF pipeline.
        """

        self.genome_path = genome
        self.annotation_path = annotation
        self.out_prefix = out_prefix
        self.orf_prefix = orf_prefix
        self.start_codons = [x.strip().upper() for x in start_codons.split(",")]
        self.min_aa = min_aa
        self.max_aa = max_aa
        self.scan_strand = scan_strand
        self.kozak_up = kozak_up
        self.kozak_down = kozak_down
        self.mark_overlap = mark_overlap
        self.remove_discarded = remove_discarded
        self.include_stop = include_stop
        self.threads = max(1, int(threads))
        self.records = []

    def run(self) -> None:
        """
        Run the complete smORF scanning pipeline.
        """

        genome = FastaParser.read_fasta(self.genome_path)
        transcripts = GenePredParser.read_genepred(self.annotation_path)

        if self.threads == 1:
            self._run_single_process(genome, transcripts)
        else:
            self._run_multi_process(genome, transcripts)

        self.write_outputs()

    def _run_single_process(self, genome, transcripts) -> None:
        """
        Run smORF scanning in single-process mode.
        """

        scanner = ORFScanner(
            start_codons=self.start_codons,
            min_aa=self.min_aa,
            max_aa=self.max_aa,
            scan_strand=self.scan_strand,
            kozak_up=self.kozak_up,
            kozak_down=self.kozak_down,
            include_stop=self.include_stop,
        )

        total_tx = len(transcripts)
        orf_index = 1

        for idx, tx in enumerate(transcripts, start=1):
            # Print scanning progress.
            print(
                "[smORFScanner] [{}/{}] Scanning gene={}, transcript={}, chrom={}, strand={}".format(
                    idx,
                    total_tx,
                    tx.gene_id,
                    tx.transcript_id,
                    tx.chrom,
                    tx.strand,
                ),
                flush=True,
            )

            # Reconstruct spliced transcript sequence before ORF scanning.
            CoordinateMapper.build_transcript_sequence(tx, genome)

            # Scan candidate ORFs from the transcript sequence.
            tx_records = scanner.scan_transcript(tx)

            # Assign ORF category labels.
            ORFClassifier.classify(tx, tx_records)

            # Mark nested or overlapping ORFs if requested.
            if self.mark_overlap:
                ORFOverlapMarker.mark(tx_records)

            # Remove same-frame internal ORFs if requested.
            if self.remove_discarded:
                tx_records = [x for x in tx_records if x.priority != "discarded"]

            # Assign stable ORF IDs.
            for rec in tx_records:
                rec.orf_id = "{}{:08d}".format(self.orf_prefix, orf_index)
                orf_index += 1

            self.records.extend(tx_records)

    def _run_multi_process(self, genome, transcripts) -> None:
        """
        Run smORF scanning in multiprocessing mode.
        """

        total_tx = len(transcripts)

        config = {
            "start_codons": self.start_codons,
            "min_aa": self.min_aa,
            "max_aa": self.max_aa,
            "scan_strand": self.scan_strand,
            "kozak_up": self.kozak_up,
            "kozak_down": self.kozak_down,
            "include_stop": self.include_stop,
            "mark_overlap": self.mark_overlap,
            "remove_discarded": self.remove_discarded,
        }

        print(
            "[smORFScanner] Running in multiprocessing mode with {} workers.".format(
                self.threads
            ),
            flush=True,
        )

        results_by_index = {}

        with ProcessPoolExecutor(
            max_workers=self.threads,
            initializer=_init_worker,
            initargs=(genome, config),
        ) as executor:
            future_to_index = {
                executor.submit(_scan_transcript_worker, (idx, tx)): idx
                for idx, tx in enumerate(transcripts, start=1)
            }

            finished = 0

            for future in as_completed(future_to_index):
                idx, transcript_id, gene_id, tx_records = future.result()
                results_by_index[idx] = tx_records

                finished += 1

                # Print completed transcript progress.
                print(
                    "[smORFScanner] [{}/{}] Finished gene={}, transcript={}, ORFs={}".format(
                        finished,
                        total_tx,
                        gene_id,
                        transcript_id,
                        len(tx_records),
                    ),
                    flush=True,
                )

        # Rebuild records in original transcript order and assign stable ORF IDs.
        orf_index = 1

        for idx in range(1, total_tx + 1):
            tx_records = results_by_index.get(idx, [])

            for rec in tx_records:
                rec.orf_id = "{}{:08d}".format(self.orf_prefix, orf_index)
                orf_index += 1

            self.records.extend(tx_records)

    def write_outputs(self) -> None:
        """
        Write all output files.
        """

        GenePredWriter.write("{}.genePred".format(self.out_prefix), self.records)
        MessageWriter.write("{}.message.txt".format(self.out_prefix), self.records)
        FastaWriter.write_nt("{}.nt.fa".format(self.out_prefix), self.records)
        FastaWriter.write_pep("{}.pep.fa".format(self.out_prefix), self.records)
