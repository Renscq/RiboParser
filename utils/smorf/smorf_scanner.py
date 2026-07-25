#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.003
# Function: Scan complete smORFs with annotation-boundary mORF rescue.
# Input: Reconstructed Transcript objects and ORF scanning parameters.
# Output: Transcript-coordinate and genome-coordinate ORFRecord objects.

"""High-performance transcript-centric ORF scanner.

The scanner translates every reading frame exactly once. Candidate ORFs then
reuse slices of the cached frame peptide instead of translating each nested ORF
independently. Ambiguous-codon counts are accumulated while scanning each
stop-delimited frame segment, avoiding repeated codon traversal.
"""

from __future__ import annotations

from collections.abc import Sequence

from .smorf_coordinate import CoordinateMapper
from .smorf_models import ORFRecord, Transcript
from .smorf_sequence import GENETIC_CODE, STOP_CODONS, SeqUtils

DNA_BASES = frozenset("ACGT")
VALID_SCAN_STRANDS = frozenset({"sense", "antisense", "both"})


class ORFScanner:
    """Scan candidate ORFs from reconstructed transcript sequences."""

    def __init__(
        self,
        start_codons: Sequence[str],
        min_aa: int = 8,
        max_aa: int = 10000,
        scan_strand: str = "sense",
        kozak_up: int = 6,
        kozak_down: int = 6,
        include_stop: bool = False,
        keep_partial: bool = False,
        allow_ambiguous: bool = False,
    ) -> None:
        """Initialize and validate scanner configuration."""
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
                "Invalid start codon(s): " + ", ".join(invalid_start_codons)
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
        """Scan one reconstructed transcript."""
        if not tx.tx_seq:
            raise ValueError(
                "Transcript sequence has not been reconstructed: "
                f"{tx.transcript_id}"
            )
        if len(tx.tx_seq) != tx.transcript_length():
            raise ValueError(
                f"Transcript sequence length mismatch for {tx.transcript_id}."
            )

        records: list[ORFRecord] = []
        if self.scan_strand in {"sense", "both"}:
            sense_records = self._scan_sequence(
                tx=tx,
                seq=tx.tx_seq,
                source_strand="sense",
            )
            if self._ensure_annotated_morf(
                tx=tx,
                seq=tx.tx_seq,
                records=sense_records,
            ):
                sense_records.sort(
                    key=lambda record: (
                        record.frame,
                        record.tx_orf_start,
                        record.tx_orf_end,
                    )
                )
            records.extend(sense_records)

        if self.scan_strand in {"antisense", "both"}:
            records.extend(
                self._scan_sequence(
                    tx=tx,
                    seq=SeqUtils.reverse_complement(tx.tx_seq),
                    source_strand="antisense",
                )
            )
        return records

    @staticmethod
    def _normalize_sequence(seq: str) -> str:
        """Normalize a nucleotide sequence without copying clean input."""
        if seq.isupper() and "U" not in seq:
            return seq
        return seq.upper().replace("U", "T")

    def summarize_transcript(
        self,
        tx: Transcript,
        remove_discarded: bool = False,
    ) -> tuple[int, int, int]:
        """Summarize retained ORFs without constructing output objects.

        Parameters
        ----------
        tx : Transcript
            Reconstructed transcript.
        remove_discarded : bool, optional
            Exclude complete same-frame internal ORFs that the classifier
            labels as discarded.

        Returns
        -------
        tuple of int
            Retained ORF count, total nucleotide output length, and total
            peptide output length.
        """
        if not tx.tx_seq:
            raise ValueError(
                "Transcript sequence has not been reconstructed: "
                f"{tx.transcript_id}"
            )
        if len(tx.tx_seq) != tx.transcript_length():
            raise ValueError(
                f"Transcript sequence length mismatch for {tx.transcript_id}."
            )

        summary = [0, 0, 0, 0]
        if self.scan_strand in {"sense", "both"}:
            self._summarize_sequence(
                tx=tx,
                seq=tx.tx_seq,
                source_strand="sense",
                remove_discarded=remove_discarded,
                summary=summary,
            )
        if self.scan_strand in {"antisense", "both"}:
            self._summarize_sequence(
                tx=tx,
                seq=SeqUtils.reverse_complement(tx.tx_seq),
                source_strand="antisense",
                remove_discarded=remove_discarded,
                summary=summary,
            )
        if self.scan_strand in {"sense", "both"}:
            self._summarize_missing_annotated_morf(
                tx=tx,
                seq=tx.tx_seq,
                summary=summary,
            )

        return summary[0], summary[1], summary[2]

    def _summarize_sequence(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
        remove_discarded: bool,
        summary: list[int],
    ) -> None:
        """Accumulate output-size statistics for one oriented sequence."""
        sequence = self._normalize_sequence(seq)
        sequence_length = len(sequence)

        for frame in range(3):
            frame_end = frame + ((sequence_length - frame) // 3) * 3
            if frame_end - frame < 3:
                continue

            segment_starts: list[tuple[int, int, int]] = []
            ambiguous_count = 0
            codon_index = 0

            for position in range(frame, frame_end, 3):
                codon = sequence[position:position + 3]
                if codon not in GENETIC_CODE:
                    ambiguous_count += 1

                if codon in STOP_CODONS:
                    self._summarize_segment(
                        tx=tx,
                        source_strand=source_strand,
                        sequence_length=sequence_length,
                        segment_starts=segment_starts,
                        stop_codon_index=codon_index,
                        scan_end=position + 3,
                        ambiguous_count=ambiguous_count,
                        completeness="complete",
                        remove_discarded=remove_discarded,
                        summary=summary,
                    )
                    segment_starts = []
                    ambiguous_count = 0
                elif codon in self.start_codons:
                    segment_starts.append(
                        (position, codon_index, ambiguous_count)
                    )
                codon_index += 1

            if self.keep_partial and segment_starts:
                self._summarize_segment(
                    tx=tx,
                    source_strand=source_strand,
                    sequence_length=sequence_length,
                    segment_starts=segment_starts,
                    stop_codon_index=(frame_end - frame) // 3,
                    scan_end=frame_end,
                    ambiguous_count=ambiguous_count,
                    completeness="3prime_partial",
                    remove_discarded=remove_discarded,
                    summary=summary,
                )

    def _summarize_segment(
        self,
        tx: Transcript,
        source_strand: str,
        sequence_length: int,
        segment_starts: list[tuple[int, int, int]],
        stop_codon_index: int,
        scan_end: int,
        ambiguous_count: int,
        completeness: str,
        remove_discarded: bool,
        summary: list[int],
    ) -> None:
        """Accumulate valid candidates from one stop-delimited segment."""
        for start_position, start_codon_index, ambiguous_at_start in segment_starts:
            aa_length = stop_codon_index - start_codon_index
            if aa_length < self.min_aa or aa_length > self.max_aa:
                continue

            candidate_ambiguous_count = (
                ambiguous_count - ambiguous_at_start
            )
            if candidate_ambiguous_count and not self.allow_ambiguous:
                continue

            tx_start, tx_end, _orf_strand = (
                self._scan_to_transcript_interval(
                    tx=tx,
                    sequence_length=sequence_length,
                    source_strand=source_strand,
                    scan_start=start_position,
                    scan_end=scan_end,
                )
            )
            if remove_discarded and self._is_discarded_internal_orf(
                tx=tx,
                source_strand=source_strand,
                tx_start=tx_start,
                tx_end=tx_end,
            ):
                continue

            summary[0] += 1
            summary[1] += scan_end - start_position
            summary[2] += aa_length + (
                1 if completeness == "complete" and self.include_stop else 0
            )
            if (
                source_strand == "sense"
                and self._matches_annotated_interval(
                    tx=tx,
                    tx_start=tx_start,
                    tx_end=tx_end,
                )
            ):
                summary[3] = 1

    def _annotated_morf_spec(
        self,
        tx: Transcript,
        seq: str,
    ) -> tuple[int, int] | None:
        """Return the annotation-supported CDS start and stop positions.

        The genePred ecosystem contains both conventions in which ``cdsEnd``
        includes the terminal stop codon and conventions in which ``cdsEnd``
        points immediately before that stop codon. Both are accepted.

        Parameters
        ----------
        tx : Transcript
            Reconstructed coding transcript.
        seq : str
            Sense transcript sequence.

        Returns
        -------
        tuple of int or None
            Start-codon position and stop-codon position in transcript
            coordinates, or ``None`` when no complete annotation-supported ORF
            can be constructed.
        """
        if (
            not tx.is_coding()
            or tx.cds_tx_start is None
            or tx.cds_tx_end is None
        ):
            return None

        sequence = self._normalize_sequence(seq)
        cds_start = int(tx.cds_tx_start)
        cds_end = int(tx.cds_tx_end)

        if cds_start < 0 or cds_end <= cds_start:
            return None
        if cds_start + 3 > len(sequence):
            return None
        if sequence[cds_start:cds_start + 3] not in self.start_codons:
            return None

        stop_position: int | None = None

        if (
            cds_end <= len(sequence)
            and cds_end - cds_start >= 6
            and (cds_end - cds_start) % 3 == 0
            and sequence[cds_end - 3:cds_end] in STOP_CODONS
        ):
            stop_position = cds_end - 3
        elif (
            cds_end + 3 <= len(sequence)
            and cds_end - cds_start >= 3
            and (cds_end - cds_start) % 3 == 0
            and sequence[cds_end:cds_end + 3] in STOP_CODONS
        ):
            stop_position = cds_end

        if stop_position is None:
            return None

        first_stop = self._find_stop(sequence, cds_start + 3)
        if first_stop != stop_position:
            return None

        aa_length = (stop_position - cds_start) // 3
        if aa_length < self.min_aa or aa_length > self.max_aa:
            return None

        nucleotide_sequence = sequence[cds_start:stop_position + 3]
        ambiguous_codon_count = self._ambiguous_codon_count(
            nucleotide_sequence
        )
        if ambiguous_codon_count and not self.allow_ambiguous:
            return None

        return cds_start, stop_position

    @staticmethod
    def _matches_annotated_interval(
        tx: Transcript,
        tx_start: int,
        tx_end: int,
    ) -> bool:
        """Return whether an ORF interval matches an annotated CDS boundary."""
        if tx.cds_tx_start is None or tx.cds_tx_end is None:
            return False
        if tx_start != tx.cds_tx_start:
            return False
        return tx_end in {
            tx.cds_tx_end,
            tx.cds_tx_end + 3,
        }

    def _ensure_annotated_morf(
        self,
        tx: Transcript,
        seq: str,
        records: list[ORFRecord],
    ) -> bool:
        """Rescue one complete annotated mORF missing from de novo output.

        Returns
        -------
        bool
            ``True`` when a record was appended.
        """
        specification = self._annotated_morf_spec(tx=tx, seq=seq)
        if specification is None:
            return False

        start_position, stop_position = specification
        expected_end = stop_position + 3
        if any(
            record.source_strand == "sense"
            and record.tx_orf_start == start_position
            and record.tx_orf_end == expected_end
            for record in records
        ):
            return False

        record = self._build_record(
            tx=tx,
            seq=self._normalize_sequence(seq),
            source_strand="sense",
            frame=start_position % 3,
            start_position=start_position,
            stop_position=stop_position,
        )
        if record is None:
            return False

        records.append(record)
        return True

    def _summarize_missing_annotated_morf(
        self,
        tx: Transcript,
        seq: str,
        summary: list[int],
    ) -> None:
        """Add a rescued annotated mORF to the lightweight summary pass."""
        if summary[3]:
            return

        specification = self._annotated_morf_spec(tx=tx, seq=seq)
        if specification is None:
            return

        start_position, stop_position = specification
        aa_length = (stop_position - start_position) // 3
        summary[0] += 1
        summary[1] += stop_position + 3 - start_position
        summary[2] += aa_length + int(self.include_stop)
        summary[3] = 1

    @staticmethod
    def _is_discarded_internal_orf(
        tx: Transcript,
        source_strand: str,
        tx_start: int,
        tx_end: int,
    ) -> bool:
        """Return whether classification marks a candidate as discarded."""
        if source_strand == "antisense":
            return False
        if (
            not tx.is_coding()
            or tx.cds_tx_start is None
            or tx.cds_tx_end is None
        ):
            return False

        cds_start = tx.cds_tx_start
        cds_end = tx.cds_tx_end
        if tx_start == cds_start and tx_end == cds_end:
            return False
        return (
            tx_start >= cds_start
            and tx_end <= cds_end
            and (tx_start - cds_start) % 3 == 0
        )

    def _scan_sequence(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
    ) -> list[ORFRecord]:
        """Scan all three reading frames of one oriented sequence."""
        sequence = self._normalize_sequence(seq)
        records: list[ORFRecord] = []
        for frame in range(3):
            records.extend(
                self._scan_frame(
                    tx=tx,
                    seq=sequence,
                    source_strand=source_strand,
                    frame=frame,
                )
            )
        return records

    def _scan_frame(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
        frame: int,
    ) -> list[ORFRecord]:
        """Scan one frame using a single translation pass."""
        frame_end = frame + ((len(seq) - frame) // 3) * 3
        if frame_end - frame < 3:
            return []

        positions = range(frame, frame_end, 3)
        codons = [seq[position:position + 3] for position in positions]
        frame_peptide = "".join(
            GENETIC_CODE.get(codon, "X") for codon in codons
        )

        records: list[ORFRecord] = []
        segment_starts: list[tuple[int, int, int]] = []
        ambiguous_count = 0

        for codon_index, codon in enumerate(codons):
            position = frame + codon_index * 3
            if codon not in GENETIC_CODE:
                ambiguous_count += 1

            if codon in STOP_CODONS:
                if segment_starts:
                    records.extend(
                        self._emit_complete_segment(
                            tx=tx,
                            seq=seq,
                            source_strand=source_strand,
                            frame=frame,
                            stop_position=position,
                            stop_codon_index=codon_index,
                            frame_peptide=frame_peptide,
                            segment_starts=segment_starts,
                            ambiguous_count=ambiguous_count,
                        )
                    )
                segment_starts = []
                ambiguous_count = 0
                continue

            if codon in self.start_codons:
                segment_starts.append(
                    (position, codon_index, ambiguous_count)
                )

        if self.keep_partial and segment_starts:
            records.extend(
                self._emit_partial_segment(
                    tx=tx,
                    seq=seq,
                    source_strand=source_strand,
                    frame=frame,
                    frame_end=frame_end,
                    frame_peptide=frame_peptide,
                    segment_starts=segment_starts,
                    ambiguous_count=ambiguous_count,
                )
            )

        return records

    def _emit_complete_segment(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
        frame: int,
        stop_position: int,
        stop_codon_index: int,
        frame_peptide: str,
        segment_starts: list[tuple[int, int, int]],
        ambiguous_count: int,
    ) -> list[ORFRecord]:
        """Build complete ORFs ending at one shared stop codon."""
        records: list[ORFRecord] = []
        orf_end_scan = stop_position + 3
        stop_codon = seq[stop_position:orf_end_scan]

        for start_position, start_codon_index, ambiguous_at_start in segment_starts:
            aa_length = stop_codon_index - start_codon_index
            if aa_length < self.min_aa or aa_length > self.max_aa:
                continue

            ambiguous_codon_count = ambiguous_count - ambiguous_at_start
            if ambiguous_codon_count and not self.allow_ambiguous:
                continue

            peptide_end = stop_codon_index + int(self.include_stop)
            peptide_sequence = frame_peptide[
                start_codon_index:peptide_end
            ]
            record = self._create_record(
                tx=tx,
                seq=seq,
                source_strand=source_strand,
                frame=frame,
                start_position=start_position,
                orf_end_scan=orf_end_scan,
                stop_codon=stop_codon,
                aa_length=aa_length,
                peptide_sequence=peptide_sequence,
                completeness="complete",
                ambiguous_codon_count=ambiguous_codon_count,
            )
            if record is not None:
                records.append(record)

        return records

    def _emit_partial_segment(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
        frame: int,
        frame_end: int,
        frame_peptide: str,
        segment_starts: list[tuple[int, int, int]],
        ambiguous_count: int,
    ) -> list[ORFRecord]:
        """Build optional 3-prime partial ORFs from the final frame segment."""
        records: list[ORFRecord] = []
        frame_codon_count = len(frame_peptide)

        for start_position, start_codon_index, ambiguous_at_start in segment_starts:
            aa_length = frame_codon_count - start_codon_index
            if aa_length < self.min_aa or aa_length > self.max_aa:
                continue

            ambiguous_codon_count = ambiguous_count - ambiguous_at_start
            if ambiguous_codon_count and not self.allow_ambiguous:
                continue

            record = self._create_record(
                tx=tx,
                seq=seq,
                source_strand=source_strand,
                frame=frame,
                start_position=start_position,
                orf_end_scan=frame_end,
                stop_codon="NA",
                aa_length=aa_length,
                peptide_sequence=frame_peptide[start_codon_index:],
                completeness="3prime_partial",
                ambiguous_codon_count=ambiguous_codon_count,
            )
            if record is not None:
                records.append(record)

        return records

    def _create_record(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
        frame: int,
        start_position: int,
        orf_end_scan: int,
        stop_codon: str,
        aa_length: int,
        peptide_sequence: str,
        completeness: str,
        ambiguous_codon_count: int,
    ) -> ORFRecord | None:
        """Create one output-ready ORF record."""
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

        nucleotide_sequence = seq[start_position:orf_end_scan]
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
            kozak_seq=self._extract_kozak(
                seq=seq,
                start=start_position,
                up=self.kozak_up,
                down=self.kozak_down,
            ),
            completeness=completeness,
            kozak_start_index=self.kozak_up,
            ambiguous_codon_count=ambiguous_codon_count,
        )

    def _frame_candidates(
        self,
        seq: str,
        frame: int,
    ) -> list[tuple[int, int | None]]:
        """Return start positions and their nearest downstream stops."""
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
        """Compatibility path for callers building one candidate directly."""
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
            peptide_sequence = (
                translated_sequence
                if self.include_stop
                else translated_sequence[:-1]
            )
        else:
            peptide_sequence = translated_sequence

        return self._create_record(
            tx=tx,
            seq=seq,
            source_strand=source_strand,
            frame=frame,
            start_position=start_position,
            orf_end_scan=orf_end_scan,
            stop_codon=stop_codon,
            aa_length=aa_length,
            peptide_sequence=peptide_sequence,
            completeness=completeness,
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
        """Convert a scanned interval to transcript coordinates."""
        if source_strand == "sense":
            return scan_start, scan_end, tx.strand
        tx_start = sequence_length - scan_end
        tx_end = sequence_length - scan_start
        orf_strand = "-" if tx.strand == "+" else "+"
        return tx_start, tx_end, orf_strand

    @staticmethod
    def _ambiguous_codon_count(seq: str) -> int:
        """Count codons containing at least one non-ACGT base."""
        return sum(
            seq[position:position + 3] not in GENETIC_CODE
            for position in range(0, len(seq), 3)
        )

    @staticmethod
    def _find_stop(seq: str, start: int) -> int | None:
        """Find the first downstream in-frame stop codon."""
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
        """Extract a fixed-width, boundary-padded Kozak context."""
        left_start = max(0, start - up)
        left_sequence = seq[left_start:start]
        left_padding = "N" * (up - len(left_sequence))
        right_end = min(len(seq), start + 3 + down)
        right_sequence = seq[start:right_end]
        right_padding = "N" * (3 + down - len(right_sequence))
        return left_padding + left_sequence + right_sequence + right_padding
