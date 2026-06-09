# 4 Workflow overview

The workflow is organized as a complete tutorial.

## Main workflow

| Step | Page | Purpose |
|---|---|---|
| 4.1 | Reference preparation | Build genome, transcriptome, Bowtie, STAR, RSEM, and normalized RiboParser references |
| 4.2 | Raw data download | Download public or user-provided RNA-seq/Ribo-seq data |
| 4.3 | Raw data cleaning | Trim adapters and filter reads |
| 4.4 | Alignment and quantification | Classify reads, align to genome, and quantify expression |
| 4.5 | Quality control | Check library quality and generate density matrices |
| 4.6 | Gene-level analysis | Quantification, coverage, correlation, and density retrieval |
| 4.7 | Codon-level analysis | Pausing, occupancy, decoding time, selection time, CoV, odds ratio, and meta-codon |
| 4.8 | smORF analysis | Scan, filter, evaluate, and integrate smORF candidates |
| 4.9 | SeRP analysis | SeRP overlap, peak, and property analysis |

## Example dataset

The example workflow uses `GSE67387`, but the same structure can be applied to any RNA-seq/Ribo-seq project.
