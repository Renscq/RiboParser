# Author: Rensc
# date: 2026-05-21

"""
Sequence utility functions for smORF scanning.

This file provides:
1. Standard stop codons.
2. Standard nuclear genetic code.
3. Reverse-complement and translation functions.
"""

# Standard stop codons.
STOP_CODONS = {"TAA", "TAG", "TGA"}


# Standard genetic code.
GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class SeqUtils:
    """
    Basic sequence operation utilities.
    """

    # Translation table for reverse-complement conversion.
    COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """
        Return the reverse-complement sequence.

        Parameters
        ----------
        seq : str
            Input nucleotide sequence.

        Returns
        -------
        str
            Reverse-complement sequence in uppercase.
        """

        return seq.translate(SeqUtils.COMP)[::-1].upper()

    @staticmethod
    def translate(seq: str) -> str:
        """
        Translate nucleotide sequence into peptide sequence.

        Unknown or ambiguous codons are translated as X.

        Parameters
        ----------
        seq : str
            Input nucleotide sequence.

        Returns
        -------
        str
            Translated peptide sequence.
        """

        aa = []
        seq = seq.upper()

        # Translate only complete codons.
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            aa.append(GENETIC_CODE.get(codon, "X"))

        return "".join(aa)
