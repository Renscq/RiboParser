# Author: Rensc
# date: 2026-05-21

"""
ORF scanner for transcript sequences.

This scanner supports:
1. Sense strand scanning.
2. Antisense strand scanning.
3. Custom start codons.
4. Standard stop-codon terminated ORFs.
5. 3' partial ORFs without downstream stop codons.
6. Kozak sequence extraction.
"""

from typing import List, Optional
from .smorf_models import Transcript, ORFRecord
from .smorf_sequence import SeqUtils, STOP_CODONS
from .smorf_coordinate import CoordinateMapper


class ORFScanner:
    """
    Scan ORFs from reconstructed transcript sequences.
    """

    def __init__(
        self,
        start_codons: List[str],
        min_aa: int = 8,
        max_aa: int = 10000,
        scan_strand: str = "sense",
        kozak_up: int = 6,
        kozak_down: int = 6,
        include_stop: bool = False,
    ):
        """
        Initialize ORFScanner.

        Parameters
        ----------
        start_codons : list
            Candidate start codons.
        min_aa : int
            Minimum peptide length.
        max_aa : int
            Maximum peptide length.
        scan_strand : str
            sense, antisense, or both.
        kozak_up : int
            Upstream nucleotide length for Kozak sequence.
        kozak_down : int
            Downstream nucleotide length after start codon.
        include_stop : bool
            Whether to keep stop symbol in peptide sequence.
        """

        self.start_codons = set(x.upper() for x in start_codons)
        self.min_aa = min_aa
        self.max_aa = max_aa
        self.scan_strand = scan_strand
        self.kozak_up = kozak_up
        self.kozak_down = kozak_down
        self.include_stop = include_stop

    def scan_transcript(self, tx: Transcript) -> List[ORFRecord]:
        """
        Scan ORFs from one transcript.

        Parameters
        ----------
        tx : Transcript
            Transcript object with reconstructed sequence.

        Returns
        -------
        list
            List of ORFRecord objects.
        """

        records = []

        if self.scan_strand in {"sense", "both"}:
            records.extend(self._scan_sequence(tx, tx.tx_seq, "sense"))

        if self.scan_strand in {"antisense", "both"}:
            rc_seq = SeqUtils.reverse_complement(tx.tx_seq)
            records.extend(self._scan_sequence(tx, rc_seq, "antisense"))

        return records

    def _scan_sequence(
        self,
        tx: Transcript,
        seq: str,
        source_strand: str,
    ) -> List[ORFRecord]:
        """
        Scan all ORFs from a sequence in three reading frames.

        Parameters
        ----------
        tx : Transcript
            Source transcript.
        seq : str
            Sequence to scan.
        source_strand : str
            sense or antisense.

        Returns
        -------
        list
            List of ORFRecord objects.
        """

        records = []
        seq_len = len(seq)

        # Scan three reading frames.
        for frame in range(3):
            i = frame

            while i <= seq_len - 3:
                codon = seq[i:i + 3]

                # Start a candidate ORF once a valid start codon is found.
                if codon in self.start_codons:
                    stop_pos = self._find_stop(seq, i + 3)

                    if stop_pos is not None:
                        orf_start_scan = i
                        orf_end_scan = stop_pos + 3
                        stop_codon = seq[stop_pos:stop_pos + 3]
                        completeness = "complete"
                    else:
                        # Keep 3' partial ORFs if no in-frame stop codon is found.
                        orf_start_scan = i
                        orf_end_scan = seq_len - ((seq_len - i) % 3)
                        stop_codon = "NA"
                        completeness = "3prime_partial"

                    nt_seq = seq[orf_start_scan:orf_end_scan]
                    aa_seq = SeqUtils.translate(nt_seq)

                    if aa_seq.endswith("*") and not self.include_stop:
                        pep_seq = aa_seq[:-1]
                    else:
                        pep_seq = aa_seq

                    aa_len = len(pep_seq.replace("*", ""))

                    # Filter ORFs by peptide length.
                    if aa_len < self.min_aa or aa_len > self.max_aa:
                        i += 3
                        continue

                    # Convert antisense coordinates back to transcript coordinates.
                    if source_strand == "sense":
                        tx_orf_start = orf_start_scan
                        tx_orf_end = orf_end_scan
                        orf_strand = tx.strand
                    else:
                        tx_orf_start = seq_len - orf_end_scan
                        tx_orf_end = seq_len - orf_start_scan
                        orf_strand = "-" if tx.strand == "+" else "+"

                    exon_starts, exon_ends = CoordinateMapper.tx_interval_to_genomic_blocks(
                        tx,
                        tx_orf_start,
                        tx_orf_end,
                    )

                    if not exon_starts:
                        i += 3
                        continue

                    genomic_start = min(exon_starts)
                    genomic_end = max(exon_ends)

                    kozak_seq = self._extract_kozak(
                        seq,
                        orf_start_scan,
                        self.kozak_up,
                        self.kozak_down,
                    )

                    records.append(
                        ORFRecord(
                            orf_id="TEMP",
                            transcript_id=tx.transcript_id,
                            gene_id=tx.gene_id,
                            chrom=tx.chrom,
                            strand=orf_strand,
                            source_strand=source_strand,
                            frame=frame,
                            tx_orf_start=tx_orf_start,
                            tx_orf_end=tx_orf_end,
                            genomic_start=genomic_start,
                            genomic_end=genomic_end,
                            exon_starts=exon_starts,
                            exon_ends=exon_ends,
                            start_codon=codon,
                            stop_codon=stop_codon,
                            nt_length=len(nt_seq),
                            aa_length=aa_len,
                            nt_seq=nt_seq,
                            pep_seq=pep_seq,
                            kozak_seq=kozak_seq,
                            completeness=completeness,
                        )
                    )

                i += 3

        return records

    @staticmethod
    def _find_stop(seq: str, start: int) -> Optional[int]:
        """
        Find the first in-frame stop codon downstream of a start codon.

        Parameters
        ----------
        seq : str
            Sequence to search.
        start : int
            Search start position.

        Returns
        -------
        int or None
            Stop codon position if found, otherwise None.
        """

        for i in range(start, len(seq) - 2, 3):
            if seq[i:i + 3] in STOP_CODONS:
                return i

        return None

    @staticmethod
    def _extract_kozak(seq: str, start: int, up: int, down: int) -> str:
        """
        Extract the Kozak-like sequence around the start codon.

        Parameters
        ----------
        seq : str
            Scanned sequence.
        start : int
            Start codon position.
        up : int
            Upstream length.
        down : int
            Downstream length after the start codon.

        Returns
        -------
        str
            Extracted Kozak-like sequence.
        """

        left = max(0, start - up)
        right = min(len(seq), start + 3 + down)

        return seq[left:right]
