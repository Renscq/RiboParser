#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-25
# Version: 0.2.8.24-dev.003
# Function: Classify ORFs with compatible annotated-CDS terminal boundaries.
# Input: Transcript objects and predicted ORFRecord objects.
# Output: ORF category and priority annotations.

"""Classify predicted ORFs relative to the annotated CDS."""

from __future__ import annotations

from collections.abc import Sequence

from .smorf_models import ORFRecord, Transcript
from .smorf_sequence import STOP_CODONS


class ORFClassifier:
    """Classify ORFs by their position relative to the annotated CDS."""

    @staticmethod
    def _is_annotated_orf(
        tx: Transcript,
        record: ORFRecord,
    ) -> bool:
        """Return whether a record represents the annotated complete mORF.

        Both common genePred terminal-boundary conventions are supported:

        1. ``cds_tx_end`` includes the terminal stop codon.
        2. ``cds_tx_end`` points immediately before the terminal stop codon.

        Parameters
        ----------
        tx : Transcript
            Source transcript.
        record : ORFRecord
            Complete sense-strand ORF candidate.

        Returns
        -------
        bool
            Whether the candidate matches the annotated coding ORF.
        """
        if tx.cds_tx_start is None or tx.cds_tx_end is None:
            return False
        if record.source_strand != "sense":
            return False
        if record.completeness != "complete":
            return False
        if record.tx_orf_start != tx.cds_tx_start:
            return False

        if record.tx_orf_end == tx.cds_tx_end:
            return True

        return (
            record.stop_codon in STOP_CODONS
            and record.tx_orf_end == tx.cds_tx_end + 3
            and (tx.cds_tx_end - tx.cds_tx_start) % 3 == 0
        )

    @staticmethod
    def classify(
        tx: Transcript,
        records: Sequence[ORFRecord],
    ) -> None:
        """Assign category labels to ORF records.

        Parameters
        ----------
        tx : Transcript
            Source transcript.
        records : sequence of ORFRecord
            Predicted ORF records modified in place.
        """
        for record in records:
            if record.source_strand == "antisense":
                record.category = "antisense_ORF"
                continue

            if (
                not tx.is_coding()
                or tx.cds_tx_start is None
                or tx.cds_tx_end is None
            ):
                record.category = "lncORF"
                continue

            cds_start = tx.cds_tx_start
            cds_end = tx.cds_tx_end
            orf_start = record.tx_orf_start
            orf_end = record.tx_orf_end

            if ORFClassifier._is_annotated_orf(tx, record):
                record.category = "annotated_ORF"
                record.priority = "primary"
                continue

            if orf_end <= cds_start:
                record.category = "uORF"
                continue

            if orf_start >= cds_end:
                record.category = "dORF"
                continue

            if orf_start < cds_start and orf_end in {
                cds_end,
                cds_end + 3,
            }:
                record.category = "emORF"
                continue

            if orf_start >= cds_start and orf_end <= cds_end:
                cds_relative_frame = (orf_start - cds_start) % 3
                if cds_relative_frame == 0:
                    record.category = "same_frame_iORF"
                    record.priority = "discarded"
                else:
                    record.category = "iORF"
                continue

            if orf_start < cds_start and orf_end > cds_start:
                record.category = "overlap_uORF"
                continue

            if orf_start < cds_end and orf_end > cds_end:
                record.category = "overlap_dORF"
                continue

            record.category = "other_ORF"
