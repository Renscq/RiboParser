# RiboParser

<p align="center">
  <strong>RiboParser</strong> is a modular command-line toolkit for comprehensive RNA-seq and Ribo-seq data analysis.
</p>

<p align="center">
  <a href="https://renscq.github.io/RiboParser/">Documentation</a> |
  <a href="https://github.com/Renscq/RiboParser/issues">Issues</a> |
  <a href="#citation">Citation</a> |
  <a href="#license">License</a>
</p>

---

## Overview

RiboParser provides an integrated workflow for RNA-seq and ribosome profiling (Ribo-seq) analysis, including reference preparation, raw data preprocessing, alignment, quantification, quality control, RPF density construction, gene-level analysis, codon-level analysis, SeRP analysis, smORF analysis, and downstream visualization with RiboShiny.

The complete workflow includes:

1. Software installation
2. Reference file preparation
3. Raw data download
4. Raw data cleaning
5. Alignment and quantification
6. RNA-seq and Ribo-seq quality control
7. Gene-level analysis
8. Codon-level analysis
9. smORF and SeRP-related analysis

## Features

- Reference preparation for RNA-seq and Ribo-seq workflows
- RPF quality control and P-site offset detection
- RPF density generation and sample-level merging
- Periodicity, metaplot, coverage, and correlation analysis
- Gene-level quantification and comparative analysis
- Codon pausing score, codon occupancy, codon decoding time, and codon selection time analysis
- smORF scanning, filtering, evidence evaluation, and integration
- SeRP-related signal and peak analysis
- Helper scripts for FASTA, FASTQ, Bowtie logs, RSEM outputs, bedGraph, and merged Ribo-seq outputs
- Interactive downstream visualization with RiboShiny

## Installation

### conda / micromamba

```bash
conda install riboparser -c rensc
# or
micromamba install riboparser -c rensc
```

### pip

```bash
pip install riboparser
```

### GitHub source

```bash
git clone https://github.com/Renscq/RiboParser.git
cd RiboParser
pip install build
python -m build
pip install .
```

## Quick test

```bash
riboparser -v
riboparser -c
riboparser -d
riboparser -m
```

## Documentation

The full documentation is available at:

<https://renscq.github.io/RiboParser/>

## Main command groups

| Group | Scope |
|---|---|
| Ribo-quality | Reference construction, QC, offset, density, merge, periodicity, metaplot, coverage, correlation, and quantification |
| Ribo-pausing | Codon pausing, occupancy, coefficient of variation, decoding time, selection time, odds ratio, and meta-codon analysis |
| Ribo-utils | RPF retrieval, BAM-to-bigWig conversion, gene plotting, shuffling, and frame-shift-related utilities |
| RNA | RNA-seq density and offset-related processing |
| SeRP | SeRP overlap, peak, and property analysis |
| smORF | smORF scanning, filtering, evidence evaluation, and integration |
| Helper scripts | FASTA, FASTQ, Bowtie, bedGraph, RSEM, RiboCode, RiboTISH, and merged Ribo-seq result utilities |

## Citation

Ren, S., Li, Y. & Zhou, Z.  
**RiboParser/RiboShiny: An integrated platform for comprehensive analysis and visualization of ribo-seq data.**  
*Journal of Genetics and Genomics* (2025).  
DOI: **10.1016/j.jgg.2025.04.010**

## License

GPL-3.0-or-later.
