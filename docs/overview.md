# 1 Overview

RiboParser is designed for complete RNA-seq and ribosome profiling data analysis. This page only provides the conceptual overview; detailed commands are placed in the workflow pages to avoid duplication.

## What RiboParser does

The complete workflow consists of two parts: a **public pipeline** that can be handled by general bioinformatics tools, and the **RiboParser-specific analysis** described below.

Public pipeline (general tools):

- reference preparation for transcriptome-aware Ribo-seq analysis
- RNA-seq and Ribo-seq raw data cleaning
- contaminant classification against rRNA, tRNA, ncRNA, mRNA, and genome indexes
- splice-aware genome alignment
- transcriptome quantification

RiboParser-specific analysis:

- Ribo-seq quality control
- P-site offset inference
- RNA/Ribo read density construction
- merged density matrix generation
- periodicity, metaplot, coverage, and correlation analysis
- gene-level quantification and read-density retrieval
- codon-level pausing, occupancy, decoding time, selection time, variation, and odds-ratio analysis
- smORF scanning, clustering, Ribo-seq evidence evaluation, and quantification
- SeRP signal and peak analysis
- helper utilities for FASTA, FASTQ, bedGraph, Bowtie logs, RSEM tables, and merged Ribo-seq outputs


## Citation

```text
Ren, S., Li, Y. & Zhou, Z.
RiboParser/RiboShiny: An integrated platform for comprehensive analysis and visualization of ribo-seq data.
Journal of Genetics and Genomics (2025).
doi:10.1016/j.jgg.2025.04.010.
```
