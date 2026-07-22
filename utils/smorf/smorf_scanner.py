#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Scan complete and 3-prime partial ORFs from transcript sequences.
# Input: Reconstructed Transcript objects and ORF scanning parameters.
# Output: Transcript-coordinate and genome-coordinate ORFRecord objects.

"""Linear-time ORF scanner for reconstructed transcript sequences.

The historical implementation searched downstream for a stop codon every time
a start codon was encountered. A start-rich transcript could therefore require
quadratic work. This implementation scans each reading frame from right to left
once, records the nearest downstream stop, and then builds candidates in
ascending coordinate order. The scanning complexity is linear in transcript
length plus the number of emitted ORFs.

Kozak contexts are padded with ``N`` at transcript boundaries so every record
has a stable start-codon index. This avoids incorrect overlap-priority scoring
for boundary ORFs and for non-default upstream context lengths.
"""

from __future__ import annotations

from collections.abc import Sequence

from .smorf_coordinate import CoordinateMapper
from .smorf_models import ORFRecord, Transcript
from .smorf_sequence import STOP_CODONS, SeqUtils


DNA_BASES = frozenset("ACGT")
VALID_SCAN_STRANDS = frozenset({"sense", "antisense", "both"})


class ORFScanner:
    """Scan candidate ORFs from reconstructed transcript sequences.

    Args:
        start_codons: Candidate start codons.
        min_aa: Minimum peptide length, excluding a terminal stop.
        max_aa: Maximum peptide length, excluding a terminal stop.
        scan_strand: ``sense``, ``antisense``, or ``both``.
        kozak_up: Number of upstream nucleotides in the Kozak context.
        kozak_down: Number of downstream nucleotides after the start codon.
        include_stop: Retain the terminal ``*`` in peptide output.
        keep_partial: Emit 3-prime partial ORFs without an in-frame stop.
        allow_ambiguous: Retain ORFs containing non-ACGT codons. Such records
            store the number of ambiguous codons in
            ``ambiguous_codon_count``.

    Raises:
        ValueError: If a scanner parameter or start codon is invalid.
    """

    def __init__(
        self,
        start_codons: Sequence[str],
        min_aa: int = 8,
        max_aa: int = 10000,
        scan_strand: str = "sense",
        kozak_up: int = 6,
        kozak_down: int = 6,
        include_stop: bool = False,
        keep_partial: bool = True,
        allow_ambiguous: bool = False,
    ) -> None:
        """Initialize and validate scanner configuration.

        Args:
            start_codons: Candidate start codons.
            min_aa: Minimum peptide length.
            max_aa: Maximum peptide length.
            scan_strand: Transcript orientation to scan.
            kozak_up: Upstream Kozak-context length.
            kozak_down: Downstream Kozak-context length.
            include_stop: Retain the terminal stop symbol in peptide output.
            keep_partial: Retain 3-prime partial ORFs.
            allow_ambiguous: Retain ORFs containing ambiguous codons.

        Raises:
            ValueError: If any scanner parameter is invalid.
        """
        normalized_start_codons = {
            str(codon).strip().upper().replace("U", "T")
            for codon in start_codons
            if str(codon).strip()
        }
        if not normalized_start_codons:
            raise ValueError("At least one start codon is required.")

        invalid_start_codons = sorted(
            codon
            for codon in normalized_start_codons
            if len(codon) != 3 or set(codon).difference(DNA_BASES)
        )
        if invalid_start_codons:
            raise ValueError(
                "Invalid start codon(s): "
                + ", ".join(invalid_start_codons)
            )

        stop_as_start = sorted(
            normalized_start_codons.intersection(STOP_CODONS)
        )
        if stop_as_start:
            raise ValueError(
                "Stop codons cannot be used as start codons: "
                + ", ".join(stop_as_start)
            )
        if min_aa < 1:
            raise ValueError("min_aa must be >= 1.")
        if max_aa < min_aa:
            raise ValueError("max_aa must be >= min_aa.")
        if scan_strand not in VALID_SCAN_STRANDS:
            raise ValueError(
                "scan_strand must be one of: sense, antisense, both."
            )
        if kozak_up < 0 or kozak_down < 0:
            raise ValueError("Kozak flank lengths must be >= 0.")

        self.start_codons = frozenset(normalized_start_codons)
        self.min_aa = int(min_aa)
        self.max_aa = int(max_aa)
        self.scan_strand = scan_strand
        self.kozak_up = int(kozak_up)
        self.kozak_down = int(kozak_down)
        self.include_stop = bool(include_stop)
        self.keep_partial = bool(keep_partial)
        self.allow_ambiguous = bool(allow_ambiguous)

    def scan_transcript(self, tx: Transcript) -> list[ORFRecord]:
        """Scan one reconstructed transcript.

        Args:
            tx: Transcript with a non-empty ``tx_seq``.

        Returns:
            ORF records in deterministic orientation, frame, and start order.

        Raises:
            ValueError: If the transcript sequence is missing or has a length
                inconsistent with its exon structure.
        """
        if not tx.tx_seq:
            raise ValueError(
                f"Transcript sequence has not been reconstructed: "
                f"{tx.transcript_id}"
            )
        if len(tx.tx_seq) != tx.transcript_length():
            raise ValueError(
                f"Transcript sequence length mismatch for {tx.transcript_id}."
            )

        records: list[ORFRecord] = []
        if self.scan_strand in {"sense", "both"}:
            records.extend(
                self._scan_sequence(
                    tx=tx,
                    seq=tx.tx_seq,
                    source_strand="sense",
                )
            )

        if self.scan_strand in {"antisense", "both"}:
            records.extend(
                self._scan_sequence(
                    tx=tx,
                    seq=SeqUtils.reverse_complement(tx.tx_seq),
                    source_strand="antisense",
                )
            )
        return records

    def _scan_sequence(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
    ) -> list[ORFRecord]:
        """Scan all three reading frames of one oriented sequence.

        Args:
            tx: Source transcript.
            seq: Sequence in the orientation being scanned.
            source_strand: ``sense`` or ``antisense`` relative to ``tx``.

        Returns:
            Predicted ORF records.
        """
        sequence = seq.upper().replace("U", "T")
        sequence_length = len(sequence)
        records: list[ORFRecord] = []

        for frame in range(3):
            for start_position, stop_position in self._frame_candidates(
                seq=sequence,
                frame=frame,
            ):
                record = self._build_record(
                    tx=tx,
                    seq=sequence,
                    source_strand=source_strand,
                    frame=frame,
                    start_position=start_position,
                    stop_position=stop_position,
                )
                if record is not None:
                    records.append(record)

        return records

    def _frame_candidates(
        self,
        seq: str,
        frame: int,
    ) -> list[tuple[int, int | None]]:
        """Return start positions and their nearest downstream stops.

        Each reading frame is traversed once from right to left. The nearest
        stop encountered so far is associated with every upstream start codon.

        Args:
            seq: Oriented transcript sequence.
            frame: Reading frame, 0, 1, or 2.

        Returns:
            ``(start_position, stop_position)`` tuples in ascending start order.
            ``stop_position`` is ``None`` for a 3-prime partial candidate.
        """
        positions = range(frame, len(seq) - 2, 3)
        nearest_stop: int | None = None
        candidates: list[tuple[int, int | None]] = []

        for position in reversed(positions):
            codon = seq[position:position + 3]
            if codon in STOP_CODONS:
                nearest_stop = position
                continue
            if codon in self.start_codons:
                if nearest_stop is not None or self.keep_partial:
                    candidates.append((position, nearest_stop))

        candidates.reverse()
        return candidates

    def _build_record(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
        frame: int,
        start_position: int,
        stop_position: int | None,
    ) -> ORFRecord | None:
        """Build and validate one candidate ORF record.

        Args:
            tx: Source transcript.
            seq: Oriented sequence being scanned.
            source_strand: Relative scanning orientation.
            frame: Reading frame in ``seq``.
            start_position: Start-codon position in ``seq``.
            stop_position: Nearest downstream in-frame stop, if present.

        Returns:
            A validated ORF record, or ``None`` when the candidate fails length
            or ambiguity filters.
        """
        if stop_position is None:
            orf_end_scan = frame + ((len(seq) - frame) // 3) * 3
            stop_codon = "NA"
            completeness = "3prime_partial"
            aa_length = (orf_end_scan - start_position) // 3
        else:
            orf_end_scan = stop_position + 3
            stop_codon = seq[stop_position:stop_position + 3]
            completeness = "complete"
            aa_length = (orf_end_scan - start_position) // 3 - 1

        if aa_length < self.min_aa or aa_length > self.max_aa:
            return None

        nucleotide_sequence = seq[start_position:orf_end_scan]
        ambiguous_codon_count = self._ambiguous_codon_count(
            nucleotide_sequence
        )
        if ambiguous_codon_count and not self.allow_ambiguous:
            return None

        translated_sequence = SeqUtils.translate(nucleotide_sequence)
        if completeness == "complete":
            if not translated_sequence.endswith("*"):
                raise ValueError(
                    f"Terminal stop translation failed for transcript "
                    f"{tx.transcript_id} at scan position {start_position}."
                )
            peptide_sequence = (
                translated_sequence
                if self.include_stop
                else translated_sequence[:-1]
            )
        else:
            peptide_sequence = translated_sequence

        tx_orf_start, tx_orf_end, orf_strand = (
            self._scan_to_transcript_interval(
                tx=tx,
                sequence_length=len(seq),
                source_strand=source_strand,
                scan_start=start_position,
                scan_end=orf_end_scan,
            )
        )

        exon_starts, exon_ends = (
            CoordinateMapper.tx_interval_to_genomic_blocks(
                tx=tx,
                t_start=tx_orf_start,
                t_end=tx_orf_end,
            )
        )
        if not exon_starts:
            return None

        kozak_sequence = self._extract_kozak(
            seq=seq,
            start=start_position,
            up=self.kozak_up,
            down=self.kozak_down,
        )

        return ORFRecord(
            orf_id="TEMP",
            transcript_id=tx.transcript_id,
            gene_id=tx.gene_id,
            chrom=tx.chrom,
            strand=orf_strand,
            source_strand=source_strand,
            frame=frame,
            tx_orf_start=tx_orf_start,
            tx_orf_end=tx_orf_end,
            genomic_start=min(exon_starts),
            genomic_end=max(exon_ends),
            exon_starts=exon_starts,
            exon_ends=exon_ends,
            start_codon=seq[start_position:start_position + 3],
            stop_codon=stop_codon,
            nt_length=len(nucleotide_sequence),
            aa_length=aa_length,
            nt_seq=nucleotide_sequence,
            pep_seq=peptide_sequence,
            kozak_seq=kozak_sequence,
            completeness=completeness,
            kozak_start_index=self.kozak_up,
            ambiguous_codon_count=ambiguous_codon_count,
        )

    @staticmethod
    def _scan_to_transcript_interval(
        tx: Transcript,
        sequence_length: int,
        source_strand: str,
        scan_start: int,
        scan_end: int,
    ) -> tuple[int, int, str]:
        """Convert a scanned interval to transcript coordinates.

        Args:
            tx: Source transcript.
            sequence_length: Length of the oriented scan sequence.
            source_strand: ``sense`` or ``antisense``.
            scan_start: Interval start in the scanned sequence.
            scan_end: Interval end in the scanned sequence.

        Returns:
            Transcript start, transcript end, and genomic ORF strand.
        """
        if source_strand == "sense":
            return scan_start, scan_end, tx.strand

        tx_start = sequence_length - scan_end
        tx_end = sequence_length - scan_start
        orf_strand = "-" if tx.strand == "+" else "+"
        return tx_start, tx_end, orf_strand

    @staticmethod
    def _ambiguous_codon_count(seq: str) -> int:
        """Count codons containing at least one non-ACGT base.

        Args:
            seq: In-frame nucleotide sequence.

        Returns:
            Number of ambiguous codons.
        """
        return sum(
            bool(set(seq[position:position + 3]).difference(DNA_BASES))
            for position in range(0, len(seq), 3)
        )

    @staticmethod
    def _find_stop(seq: str, start: int) -> int | None:
        """Find the first downstream in-frame stop codon.

        This compatibility helper is retained for external callers. The main
        scanner no longer calls it repeatedly.

        Args:
            seq: Sequence to search.
            start: In-frame search start.

        Returns:
            Stop-codon position, or ``None`` when absent.
        """
        for position in range(start, len(seq) - 2, 3):
            if seq[position:position + 3] in STOP_CODONS:
                return position
        return None

    @staticmethod
    def _extract_kozak(
        seq: str,
        start: int,
        up: int,
        down: int,
    ) -> str:
        """Extract a fixed-width, boundary-padded Kozak context.

        Args:
            seq: Oriented scan sequence.
            start: Start-codon position.
            up: Required upstream context length.
            down: Required downstream length after the start codon.

        Returns:
            Sequence of length ``up + 3 + down``. Missing transcript-boundary
            bases are represented by ``N``.
        """
        left_start = max(0, start - up)
        left_sequence = seq[left_start:start]
        left_padding = "N" * (up - len(left_sequence))

        right_end = min(len(seq), start + 3 + down)
        right_sequence = seq[start:right_end]
        right_padding = "N" * (3 + down - len(right_sequence))

        return left_padding + left_sequence + right_sequence + right_padding
