# Author: Rensc
# date: 2026-05-21

"""
FASTA parser for genome sequence loading.

This parser supports plain FASTA and gzip-compressed FASTA files.
"""

import gzip
from typing import Dict


class FastaParser:
    """
    Read genome FASTA files into an in-memory dictionary.
    """

    @staticmethod
    def open_file(path: str):
        """
        Open a plain text or gzip-compressed FASTA file.

        Parameters
        ----------
        path : str
            Input FASTA file path.

        Returns
        -------
        file object
            Opened file handle.
        """

        if path.endswith(".gz"):
            return gzip.open(path, "rt")

        return open(path, "r")

    @staticmethod
    def read_fasta(path: str) -> Dict[str, str]:
        """
        Read FASTA file into a dictionary.

        The sequence ID is extracted from the first token after '>'.

        Parameters
        ----------
        path : str
            Input FASTA file path.

        Returns
        -------
        dict
            Dictionary of {sequence_id: sequence}.
        """

        seqs = {}
        name = None
        chunks = []

        with FastaParser.open_file(path) as handle:
            for line in handle:
                line = line.strip()

                if not line:
                    continue

                # Start a new FASTA record.
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(chunks).upper()

                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)

        # Save the final FASTA record.
        if name is not None:
            seqs[name] = "".join(chunks).upper()

        return seqs
