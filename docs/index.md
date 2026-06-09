# RiboParser

**RiboParser** is a modular toolkit for comprehensive RNA-seq and Ribo-seq data analysis.

It provides a reproducible workflow for:

- reference preparation
- raw data download
- raw data cleaning
- alignment and quantification
- RNA-seq and Ribo-seq quality control
- RPF density construction
- gene-level analysis
- codon-level analysis
- SeRP analysis
- smORF analysis
- downstream visualization with RiboShiny

## Documentation structure

The documentation is organized according to the analysis workflow:

1. Overview
2. Installation
3. Quick start
4. Workflow
5. Other toolkits
6. Performance
7. License
8. Acknowledgements

## Recommended workflow

```text
Reference preparation
→ Raw data download
→ Raw data cleaning
→ Alignment and quantification
→ Quality control
→ Gene-level analysis
→ Codon-level analysis
→ smORF / SeRP analysis
→ RiboShiny visualization
```

## Quick test

```bash
riboparser -v
riboparser -c
riboparser -d
riboparser -m
```
