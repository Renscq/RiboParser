#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Mark ORF overlap relationships and assign deterministic priorities.
# Input: Classified ORFRecord objects from one or more transcripts.
# Output: Updated overlap_type and priority fields.

"""Efficient overlap annotation for transcript-centric ORF records.

Same-frame ORFs are ranked greedily by biological priority. A record is
downgraded only when it directly overlaps an already selected higher-priority
primary ORF. This avoids both local-ranking inconsistencies and over-suppression
across transitive overlap chains.

Different-frame overlaps are marked without suppressing either ORF.
"""

from __future__ import annotations

from collections import defaultdict

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
    """Mark overlap relationships and preserve classifier-level discards."""

    @staticmethod
    def mark(records: list[ORFRecord]) -> None:
        """Mark all same-frame and different-frame overlaps.

        Args:
            records: Classified ORF records.

        Notes:
            Calling this method repeatedly is deterministic. Existing overlap
            labels are reset, while ``same_frame_iORF`` records retain the
            classifier-level ``discarded`` priority.
        """
        for record in records:
            record.overlap_type = "none"
            if record.category == "same_frame_iORF":
                record.priority = "discarded"
            elif record.priority != "discarded":
                record.priority = "primary"

        ORFOverlapMarker._mark_different_frame_overlap(records)
        ORFOverlapMarker._mark_same_frame_overlap(records)

    @staticmethod
    def _mark_same_frame_overlap(records: list[ORFRecord]) -> None:
        """Resolve directly overlapping ORFs within each reading frame.

        Args:
            records: ORF records to update.
        """
        grouped: dict[
            tuple[str, str, int],
            list[ORFRecord],
        ] = defaultdict(list)

        for record in records:
            grouped[
                (
                    record.transcript_id,
                    record.source_strand,
                    record.frame,
                )
            ].append(record)

        for items in grouped.values():
            if len(items) < 2:
                continue

            ranked = sorted(items, key=ORFOverlapMarker._priority_key)
            primaries: list[ORFRecord] = []

            for record in ranked:
                overlapping_primaries = [
                    primary
                    for primary in primaries
                    if ORFOverlapMarker._is_overlap(record, primary)
                ]

                if overlapping_primaries:
                    best_primary = min(
                        overlapping_primaries,
                        key=ORFOverlapMarker._priority_key,
                    )
                    ORFOverlapMarker._downgrade_orf(
                        record,
                        best_primary,
                    )
                    continue

                if record.priority != "discarded":
                    record.priority = "primary"
                    primaries.append(record)

            for primary in primaries:
                if any(
                    other is not primary
                    and ORFOverlapMarker._is_overlap(primary, other)
                    for other in items
                ):
                    if primary.overlap_type == "none":
                        primary.overlap_type = "same_frame_overlap"

    @staticmethod
    def _mark_different_frame_overlap(
        records: list[ORFRecord],
    ) -> None:
        """Mark different-frame overlaps using an interval sweep.

        Args:
            records: ORF records to update.
        """
        grouped: dict[
            tuple[str, str],
            list[ORFRecord],
        ] = defaultdict(list)

        for record in records:
            grouped[
                (record.transcript_id, record.source_strand)
            ].append(record)

        for items in grouped.values():
            ordered = sorted(
                items,
                key=lambda record: (
                    record.tx_orf_start,
                    record.tx_orf_end,
                    record.frame,
                ),
            )
            for index, record in enumerate(ordered):
                for other in ordered[index + 1:]:
                    if other.tx_orf_start >= record.tx_orf_end:
                        break
                    if record.frame == other.frame:
                        continue
                    if not ORFOverlapMarker._is_overlap(record, other):
                        continue
                    if record.overlap_type == "none":
                        record.overlap_type = "overlap_different_frame"
                    if other.overlap_type == "none":
                        other.overlap_type = "overlap_different_frame"

    @staticmethod
    def _downgrade_orf(
        record: ORFRecord,
        best: ORFRecord,
    ) -> None:
        """Downgrade an ORF relative to a directly overlapping primary ORF.

        Args:
            record: Lower-priority ORF.
            best: Directly overlapping higher-priority ORF.
        """
        if record.priority != "discarded":
            record.priority = "secondary"

        if ORFOverlapMarker._is_identical(record, best):
            record.overlap_type = "identical_ORF"
        elif ORFOverlapMarker._is_nested(record, best):
            if (
                record.start_codon != "ATG"
                and best.start_codon == "ATG"
            ):
                record.overlap_type = "secondary_noncanonical_start"
            else:
                record.overlap_type = "nested"
        elif record.tx_orf_end == best.tx_orf_end:
            if record.start_codon != best.start_codon:
                record.overlap_type = "alternative_start_same_stop"
            else:
                record.overlap_type = "same_stop_overlap"
        elif record.tx_orf_start == best.tx_orf_start:
            record.overlap_type = "same_start_different_stop"
        else:
            record.overlap_type = "same_frame_overlap_different_stop"

    @staticmethod
    def _select_best_orf(records: list[ORFRecord]) -> ORFRecord:
        """Return the highest-priority ORF.

        Args:
            records: Candidate ORF records.

        Returns:
            Highest-priority ORF.
        """
        if not records:
            raise ValueError("Cannot select an ORF from an empty collection.")
        return min(records, key=ORFOverlapMarker._priority_key)

    @staticmethod
    def _priority_key(record: ORFRecord) -> tuple:
        """Build a deterministic biological priority key.

        Priority order:
            1. Annotated CDS-matching ORF.
            2. Non-discarded positional category.
            3. Complete ORF.
            4. Canonical or stronger start codon.
            5. Stronger simple Kozak context.
            6. Longer peptide.
            7. More upstream start in transcript orientation.

        Args:
            record: ORF record.

        Returns:
            Tuple ordered from highest to lowest priority.
        """
        annotated_rank = (
            0
            if record.category in {"annotated_ORF", "annotated_mORF"}
            else 1
        )
        discarded_rank = (
            1
            if record.priority == "discarded"
            or record.category == "same_frame_iORF"
            else 0
        )
        completeness_rank = (
            0 if record.completeness == "complete" else 1
        )
        start_rank = START_CODON_RANK.get(record.start_codon, 99)
        kozak_rank = -ORFOverlapMarker._kozak_score(record)
        length_rank = -record.aa_length
        start_position_rank = record.tx_orf_start

        return (
            annotated_rank,
            discarded_rank,
            completeness_rank,
            start_rank,
            kozak_rank,
            length_rank,
            start_position_rank,
            record.tx_orf_end,
        )

    @staticmethod
    def _kozak_score(record: ORFRecord) -> int:
        """Calculate the simple -3/+4 Kozak score.

        Args:
            record: ORF record containing a padded Kozak context and the
                start-codon index.

        Returns:
            Integer score from 0 to 2.
        """
        sequence = record.kozak_seq.upper()
        if not sequence:
            return 0

        start_index = int(record.kozak_start_index)
        minus3_index = start_index - 3
        plus4_index = start_index + 3
        score = 0

        if (
            0 <= minus3_index < len(sequence)
            and sequence[minus3_index] in {"A", "G"}
        ):
            score += 1
        if (
            0 <= plus4_index < len(sequence)
            and sequence[plus4_index] == "G"
        ):
            score += 1
        return score

    @staticmethod
    def _is_overlap(first: ORFRecord, second: ORFRecord) -> bool:
        """Return whether two transcript intervals overlap."""
        return (
            first.tx_orf_start < second.tx_orf_end
            and first.tx_orf_end > second.tx_orf_start
        )

    @staticmethod
    def _is_nested(first: ORFRecord, second: ORFRecord) -> bool:
        """Return whether ``first`` is contained in ``second``."""
        return (
            second.tx_orf_start <= first.tx_orf_start
            and second.tx_orf_end >= first.tx_orf_end
        )

    @staticmethod
    def _is_identical(first: ORFRecord, second: ORFRecord) -> bool:
        """Return whether two transcript intervals are identical."""
        return (
            first.tx_orf_start == second.tx_orf_start
            and first.tx_orf_end == second.tx_orf_end
        )
