# RiboParser documentation

**RiboParser** is a modular command-line toolkit for comprehensive RNA-seq and Ribo-seq data analysis.

This documentation restores the detailed workflow originally described in the long README and reorganizes it into a multi-page website. The goal is to keep the GitHub repository homepage readable while preserving complete command examples, parameter descriptions, expected outputs, and result interpretation in `docs/`.

## Complete analysis workflow

1. Software installation
2. Reference file creation
3. Raw data download
4. Raw data cleaning
5. Data alignment
6. Sequencing quality analysis
7. Gene-level analysis
8. Codon-level analysis
9. smORF identification
10. Other utility modules

## Example dataset

The original tutorial uses public RNA-seq and Ribo-seq data from `GSE67387` as the demonstration dataset.

```text
Dataset: GSE67387
Reference:
Nedialkova DD, Leidel SA.
Optimization of Codon Translation Rates via tRNA Modifications Maintains Proteome Integrity.
Cell 2015 Jun 18;161(7):1606-18.
PMID: 26052047
```

## Quick test

```bash
riboparser -v
riboparser -c
riboparser -d
riboparser -m
```
