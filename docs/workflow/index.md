# 4 Workflow overview

The workflow is organized as a complete tutorial. Each step contains detailed commands, parameter explanations, and expected outputs.

## Main workflow

The first four steps belong to the **public pipeline** and can be handled by general bioinformatics tools (e.g. STAR, Bowtie, RSEM). The remaining steps are the **RiboParser-specific analysis**.

```text
# 4.1-4.4: Public pipeline (general tools)
Reference preparation
→ Raw data download
→ Raw data cleaning
→ Alignment and quantification

# 4.5-4.9: RiboParser-specific analysis
→ Quality control
→ Gene-level analysis
→ Codon-level analysis
→ smORF analysis
→ SeRP analysis
```

| Step | Page | Purpose |
|---|---|---|
| 4.1 | [Reference preparation](reference-preparation.md) | Build genome, transcriptome, Bowtie, STAR, RSEM, and normalized RiboParser references |
| 4.2 | [Raw data download](raw-data-download.md) | Download public or user-provided RNA-seq/Ribo-seq data |
| 4.3 | [Raw data cleaning](raw-data-cleaning.md) | Trim adapters and filter reads |
| 4.4 | [Alignment and quantification](alignment-and-quantification.md) | Classify reads, align to genome, and quantify expression |
| 4.5 | [Quality control](quality-control/index.md) | Check library quality and generate density matrices |
| 4.6 | [Gene-level analysis](gene-level/index.md) | Quantification, coverage, correlation, and density retrieval |
| 4.7 | [Codon-level analysis](codon-level/index.md) | Pausing, occupancy, decoding time, selection time, CoV, meta-codon, and odds ratio |
| 4.8 | [smORF analysis](smorf/index.md) | Scan, cluster, evaluate, and quantify smORF candidates |
| 4.9 | [SeRP analysis](serp/index.md) | Peak, overlap, and property analysis |

## Notes

- Steps 4.1-4.4 are shared with any standard RNA-seq/Ribo-seq project and can be handled by other general tools.
- Steps 4.5-4.9 are RiboParser-specific and should be run in order, because later modules consume the reference files, merged density matrices, and normalized annotations produced by earlier steps.
