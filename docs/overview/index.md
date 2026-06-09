# 1.1 Overview

RiboParser provides an integrated workflow for RNA-seq and ribosome profiling data analysis. The documentation is organized by practical workflow rather than by a single long README.

## Major functions

- software installation and dependency checking
- reference genome and transcriptome preparation
- Bowtie indexes for genome, mRNA, rRNA, tRNA, and ncRNA
- STAR genome index construction
- RSEM transcriptome index construction
- RNA-seq raw data cleaning, classification, alignment, and quantification
- Ribo-seq raw data cleaning, classification, alignment, and quantification
- Ribo-seq quality control
- P-site offset prediction
- RPF density generation and merging
- periodicity, metagene, coverage, and correlation analysis
- gene-level quantification
- codon-level pausing, occupancy, decoding time, selection time, CoV, and meta-codon analysis
- smORF scanning, filtering, evidence evaluation, and integration
- utility functions including shuffling, retrieval, and frame-shift detection

## Recommended strategy

Use the `Workflow` section as the main tutorial. Use `Other toolkits` as a command reference.
