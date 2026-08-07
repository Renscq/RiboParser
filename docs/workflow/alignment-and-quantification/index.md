# 4.4 Alignment and quantification

This module classifies cleaned reads by reference type, aligns the remaining reads to the genome, and quantifies transcript and gene expression.

The module is applied separately to RNA-seq and Ribo-seq data. The two pipelines use identical commands and differ only in the input directory, the output directory, and the sample prefix used when merging results.

## Workflow

Each data type is processed through the same five steps:

```text
Step 1: Classify reads with Bowtie
  ↓
Step 2: Merge Bowtie mapping statistics
  ↓
Step 3: Align mRNA reads with STAR
  ↓
Step 4: Quantify expression with RSEM
  ↓
Step 5: Merge RSEM results
```

## Read classification order

```text
rRNA → tRNA → ncRNA → mRNA → genome
```

Reads are classified against sequential references so that contaminant reads (rRNA, tRNA, ncRNA) are removed first, transcript-derived reads are then kept for quantification, and the remaining reads are aligned to the genome.

## Sections

| Section | Input | Output |
|---|---|---|
| [4.4.1 RNA-seq](rna-seq.md) | cleaned RNA-seq reads in `./sce/3.rna-seq/1.cleandata/` | RNA-seq mapping statistics, BAM files, and expression tables |
| [4.4.2 Ribo-seq](ribo-seq.md) | cleaned Ribo-seq reads in `./sce/4.ribo-seq/1.cleandata/` | Ribo-seq mapping statistics, BAM files, and expression tables |

## Output interpretation

| Output | Meaning |
|---|---|
| classified BAM files | reads assigned to rRNA, tRNA, ncRNA, mRNA, or genome |
| `Aligned.sortedByCoord.out.bam` | genome-aligned BAM file |
| `Aligned.toTranscriptome.out.bam` | transcriptome-aligned BAM file for RSEM and downstream RiboParser analyses |
| `*.genes.results` | RSEM gene-level quantification results |
| `*.isoforms.results` | RSEM isoform-level quantification results |
| `gene.*.txt`, `isoforms.*.txt` | merged quantification tables |

## Notes

- Build the STAR index and RSEM index for the reference genome before running this module.
- Run Step 1 and Step 2 before aligning, because the mapping statistics tell you whether the sequencing data meets expectations.
- The transcriptome-aligned BAM files (`Aligned.toTranscriptome.out.bam`) are the main input for the quality-control module (Section 4.5).
