# Author: Rensc
# date: 2026-05-21

"""
ORF overlap marker.

This module marks ORF overlap types and assigns priority labels.

Priority rule:
1. annotated_mORF is preferred.
2. complete ORF is preferred over partial ORF.
3. ATG start codon is preferred over non-ATG start codons.
4. Stronger Kozak context is preferred.
5. Longer ORF is preferred.
6. More upstream start site is preferred.

Important:
Different-frame overlapping ORFs are not suppressed by default.
"""

from typing import List, Dict, Tuple
from .smorf_models import ORFRecord


START_CODON_RANK = {
    "ATG": 1,
    "CTG": 2,
    "GTG": 3,
    "TTG": 4,
    "ACG": 5,
    "ATA": 6,
    "ATT": 7,
    "ATC": 8,
}


class ORFOverlapMarker:
    """
    Mark ORF overlap type and priority.
    """

    @staticmethod
    def mark(records: List[ORFRecord]) -> None:
        """
        Mark ORF overlap relationships.

        Parameters
        ----------
        records : list
            List of ORFRecord objects.
        """

        ORFOverlapMarker._mark_different_frame_overlap(records)
        ORFOverlapMarker._mark_same_frame_overlap(records)

    @staticmethod
    def _mark_same_frame_overlap(records: List[ORFRecord]) -> None:
        """
        Mark overlaps among ORFs with the same transcript, strand, and frame.
        """

        grouped: Dict[Tuple[str, str, int], List[ORFRecord]] = {}

        for rec in records:
            key = (rec.transcript_id, rec.source_strand, rec.frame)
            grouped.setdefault(key, []).append(rec)

        for _, items in grouped.items():
            items.sort(key=lambda x: (x.tx_orf_start, x.tx_orf_end))

            for rec in items:
                competitors = [
                    other for other in items
                    if other is not rec
                    and ORFOverlapMarker._is_overlap(rec, other)
                ]

                if not competitors:
                    continue

                best = ORFOverlapMarker._select_best_orf([rec] + competitors)

                if rec is best:
                    if rec.overlap_type == "none":
                        rec.overlap_type = "same_frame_overlap"
                    rec.priority = "primary"
                else:
                    ORFOverlapMarker._downgrade_orf(rec, best)

    @staticmethod
    def _mark_different_frame_overlap(records: List[ORFRecord]) -> None:
        """
        Mark different-frame overlaps without suppressing either ORF.
        """

        grouped: Dict[Tuple[str, str], List[ORFRecord]] = {}

        for rec in records:
            key = (rec.transcript_id, rec.source_strand)
            grouped.setdefault(key, []).append(rec)

        for _, items in grouped.items():
            for i, rec in enumerate(items):
                for j, other in enumerate(items):
                    if i >= j:
                        continue

                    if rec.frame == other.frame:
                        continue

                    if ORFOverlapMarker._is_overlap(rec, other):
                        if rec.overlap_type == "none":
                            rec.overlap_type = "overlap_different_frame"
                        if other.overlap_type == "none":
                            other.overlap_type = "overlap_different_frame"

    @staticmethod
    def _downgrade_orf(rec: ORFRecord, best: ORFRecord) -> None:
        """
        Downgrade an ORF according to its relationship with the selected best ORF.
        """

        rec.priority = "secondary"

        if ORFOverlapMarker._is_identical(rec, best):
            rec.overlap_type = "identical_ORF"
        elif ORFOverlapMarker._is_nested(rec, best):
            if rec.start_codon != "ATG" and best.start_codon == "ATG":
                rec.overlap_type = "secondary_noncanonical_start"
            else:
                rec.overlap_type = "nested"
        elif rec.tx_orf_end == best.tx_orf_end:
            if rec.start_codon != best.start_codon:
                rec.overlap_type = "alternative_start_same_stop"
            else:
                rec.overlap_type = "same_stop_overlap"
        elif rec.tx_orf_start == best.tx_orf_start:
            rec.overlap_type = "same_start_different_stop"
        else:
            rec.overlap_type = "same_frame_overlap_different_stop"

    @staticmethod
    def _select_best_orf(records: List[ORFRecord]) -> ORFRecord:
        """
        Select the most reliable ORF from overlapping ORFs.
        """

        return sorted(records, key=ORFOverlapMarker._priority_key)[0]

    @staticmethod
    def _priority_key(rec: ORFRecord):
        """
        Build sorting key for ORF priority.

        Lower value means higher priority.
        """

        annotated_rank = 0 if rec.category in {"annotated_ORF", "annotated_mORF"} else 1
        completeness_rank = 0 if rec.completeness == "complete" else 1
        start_rank = START_CODON_RANK.get(rec.start_codon, 99)
        kozak_rank = -ORFOverlapMarker._kozak_score(rec.kozak_seq)

        # Longer ORFs are preferred after biological confidence rules.
        length_rank = -rec.aa_length

        # More upstream start site is preferred if all other ranks are equal.
        start_position_rank = rec.tx_orf_start

        return (
            annotated_rank,
            completeness_rank,
            start_rank,
            kozak_rank,
            length_rank,
            start_position_rank,
        )

    @staticmethod
    def _kozak_score(kozak_seq: str) -> int:
        """
        Calculate a simple Kozak score.

        Rule:
        - Position -3 is A/G: +1
        - Position +4 is G: +1

        The input sequence is expected to contain:
        upstream sequence + start codon + downstream sequence.
        """

        if not kozak_seq:
            return 0

        seq = kozak_seq.upper()
        score = 0

        # Default scanner extracts 6 nt upstream + 3 nt start codon + downstream.
        # Therefore start codon begins at index 6 if full Kozak sequence exists.
        start_index = 6 if len(seq) >= 9 else max(0, len(seq) // 2 - 1)

        minus3_index = start_index - 3
        plus4_index = start_index + 3

        if 0 <= minus3_index < len(seq) and seq[minus3_index] in {"A", "G"}:
            score += 1

        if 0 <= plus4_index < len(seq) and seq[plus4_index] == "G":
            score += 1

        return score

    @staticmethod
    def _is_overlap(a: ORFRecord, b: ORFRecord) -> bool:
        """
        Check whether two ORFs overlap in transcript coordinates.
        """

        return a.tx_orf_start < b.tx_orf_end and a.tx_orf_end > b.tx_orf_start

    @staticmethod
    def _is_nested(a: ORFRecord, b: ORFRecord) -> bool:
        """
        Check whether ORF a is fully contained within ORF b.
        """

        return b.tx_orf_start <= a.tx_orf_start and b.tx_orf_end >= a.tx_orf_end

    @staticmethod
    def _is_identical(a: ORFRecord, b: ORFRecord) -> bool:
        """
        Check whether two ORFs have identical transcript coordinates.
        """

        return a.tx_orf_start == b.tx_orf_start and a.tx_orf_end == b.tx_orf_end
