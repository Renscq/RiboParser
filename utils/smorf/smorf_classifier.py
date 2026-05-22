# Author: Rensc
# date: 2026-05-21

"""
ORF category classifier.

Categories:
1. annotated_ORF: ORF exactly matches the annotated CDS.
2. uORF: ORF fully located in 5' UTR.
3. dORF: ORF fully located in 3' UTR.
4. lncORF: ORF from non-coding transcript.
5. emORF: ORF extends the annotated CDS upstream and shares CDS end.
6. iORF: ORF fully inside annotated CDS but in a different reading frame.
7. same_frame_iORF: ORF fully inside annotated CDS and in the same frame.
8. overlap_uORF: ORF overlaps 5' UTR and CDS.
9. overlap_dORF: ORF overlaps CDS and 3' UTR.
10. antisense_ORF: ORF detected from antisense transcript sequence.
"""

from typing import List
from .smorf_models import Transcript, ORFRecord


class ORFClassifier:
    """
    Classify ORFs by relative position to annotated CDS.
    """

    @staticmethod
    def classify(tx: Transcript, records: List[ORFRecord]) -> None:
        """
        Assign category labels to ORF records.

        Parameters
        ----------
        tx : Transcript
            Source transcript.
        records : list
            List of ORFRecord objects.
        """

        for rec in records:
            # Antisense ORFs are not classified by sense CDS structure.
            if rec.source_strand == "antisense":
                rec.category = "antisense_ORF"
                continue

            # Transcripts without valid CDS are treated as non-coding transcripts.
            if not tx.is_coding() or tx.cds_tx_start is None or tx.cds_tx_end is None:
                rec.category = "lncORF"
                continue

            cds_start = tx.cds_tx_start
            cds_end = tx.cds_tx_end
            orf_start = rec.tx_orf_start
            orf_end = rec.tx_orf_end

            # The annotated coding ORF must be judged before internal ORF logic.
            if orf_start == cds_start and orf_end == cds_end:
                rec.category = "annotated_ORF"
                rec.priority = "primary"
                continue

            # Upstream ORF located completely in the 5' UTR.
            if orf_end <= cds_start:
                rec.category = "uORF"
                continue

            # Downstream ORF located completely in the 3' UTR.
            if orf_start >= cds_end:
                rec.category = "dORF"
                continue

            # Extended mORF: upstream start codon with the same CDS stop site.
            if orf_start < cds_start and orf_end == cds_end:
                rec.category = "emORF"
                continue

            # Internal ORF located completely inside annotated CDS.
            if orf_start >= cds_start and orf_end <= cds_end:
                cds_relative_frame = (orf_start - cds_start) % 3

                if cds_relative_frame == 0:
                    rec.category = "same_frame_iORF"
                    rec.priority = "discarded"
                else:
                    rec.category = "iORF"

                continue

            # ORF starts in 5' UTR and overlaps CDS.
            if orf_start < cds_start and orf_end > cds_start:
                rec.category = "overlap_uORF"
                continue

            # ORF starts in CDS and extends into 3' UTR.
            if orf_start < cds_end and orf_end > cds_end:
                rec.category = "overlap_dORF"
                continue

            rec.category = "other_ORF"
