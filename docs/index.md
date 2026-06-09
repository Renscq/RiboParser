# RiboParser documentation

**RiboParser** is a modular command-line toolkit for comprehensive RNA-seq and Ribo-seq data analysis.

This documentation keeps the repository README short and moves the complete tutorial, parameter explanations, command templates, output interpretation, and helper scripts into a structured website.

## Complete workflow

```text
Reference preparation
→ Raw data download
→ Raw data cleaning
→ Alignment and quantification
→ Quality control
→ Gene-level analysis
→ Codon-level analysis
→ smORF analysis
→ SeRP analysis
→ Visualization and downstream analysis
```

## Main command groups

| Group | Commands |
|---|---|
| Ribo-quality | `rpf_Reference`, `rpf_Check`, `rpf_Digest`, `rpf_Offset`, `rpf_Density`, `rpf_Merge`, `rpf_Periodicity`, `rpf_Metaplot`, `rpf_Coverage`, `rpf_Corr`, `rpf_Quant`, `rpf_Percent` |
| RNA | `rna_Offset`, `rna_Density` |
| Ribo-pausing | `rpf_Pausing`, `rpf_Occupancy`, `rpf_CoV`, `rpf_Cumulative_CoV`, `rpf_CDT`, `rpf_CST`, `rpf_Odd_Ratio`, `rpf_Meta_Codon` |
| Ribo-utils | `rpf_Shuffle`, `rpf_Shift`, `rpf_Retrieve`, `rpf_Bam2bw`, `rpf_Geneplot` |
| smORF | `smorf_scanner`, `smorf_filter`, `smorf_evidence`, `smorf_integrate` |
| SeRP | `serp_overlap`, `serp_peak`, `serp_properties` |
| Helper scripts | FASTA, FASTQ, bedGraph, Bowtie, merge_ribo, RiboCode, RiboTISH, RSEM, Unix helpers |

## Example dataset

The original workflow uses `GSE67387` as an example dataset.

```text
Nedialkova DD, Leidel SA.
Optimization of Codon Translation Rates via tRNA Modifications Maintains Proteome Integrity.
Cell 2015 Jun 18;161(7):1606-18.
PMID: 26052047
```
