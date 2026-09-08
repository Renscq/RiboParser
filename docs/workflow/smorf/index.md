# 4.8 smORF analysis

## Function

The smORF analysis module provides a complete workflow for scanning candidate small open reading frames (smORFs), filtering and clustering candidates with sequence and Kozak-context criteria, evaluating Ribo-seq translation evidence, and quantifying P-site density for reliable smORFs.

The workflow is designed for transcript-centric smORF discovery. Candidate ORFs are first generated from a genome FASTA and genePred annotation, then filtered and evaluated using P-site density tracks from Ribo-seq data.

## Workflow

```text
smorf_scanner → smorf_cluster → smorf_evidence → smorf_quant
```

| Step | Command | Main purpose |
|---|---|---|
| 1 | `smorf_scanner` | Scan transcript-centric candidate ORFs from genome FASTA and genePred annotation. |
| 2 | `smorf_cluster` | Filter scanned ORFs by sequence and Kozak-context criteria, then cluster them into non-redundant ORF families. |
| 3 | `smorf_evidence` | Evaluate family-aware Ribo-seq translation evidence and classify reliable smORFs. |
| 4 | `smorf_quant` | Quantify raw P-site density counts for reliable smORFs from per-sample density tracks. |


## Suggested directory structure

The smORF modules are organized under `sce/5.smorf/` (see the project layout in [3 New project](../../new-project.md)). Create the smORF sub-directories (`01.scanner` to `04.quant`, matching steps 1-4 above) as follows:

```bash
mkdir -p 01.scanner 02.cluster 03.evidence 04.quant
```


## Input dependency

The smORF workflow combines two kinds of external inputs. `smorf_scanner` reconstructs transcript sequences from a genome FASTA and a genePred annotation, while `smorf_evidence` and `smorf_quant` consume strand-specific or unstranded P-site density tracks (bedGraph/WIG) that can be generated from genome-aligned Ribo-seq BAM files with `rpf_Bam2bw` ([rpf_Bam2bw](../others/Bam2bw.md)). Inside the smORF branch the modules form a strict chain — `smorf_scanner` → `smorf_cluster` → `smorf_evidence` → `smorf_quant` — with each step taking the output of the previous module. The per-module input and output are summarized below:

| Step | Main input | Main output |
|---|---|---|
| Scanner (4.8.1) | genome FASTA (`-g`) and genePred annotation (`-a`) | ORF annotation `*.genePred`, candidate ORF table `*.message.txt`, and nucleotide/peptide FASTA `*.nt.fa` / `*.pep.fa` |
| Cluster (4.8.2) | scanner message table (`-i`) and genePred annotation (`-a`) | clustered family tables `*.family.message.txt` / `*.family.members.txt`, removed-ORF table, and cluster summary |
| Evidence (4.8.3) | clustered family tables, scanner message table and ORF genePred, and a density list pointing to P-site tracks (`-l`) | family evidence table, reliable smORF table and genePred, and evidence summary |
| Quant (4.8.4) | reliable smORF genePred (`-g`) and density list (`-l`) | ORF-by-sample raw P-site count matrix `<prefix>.density_quant.txt` |

## Main evidence levels

| Evidence | Meaning |
|---|---|
| ORF sequence | Start codon, stop codon, ORF length, strand, and category. |
| Kozak context | Start-codon context scored by annotated, built-in, PWM, or sequence-derived Kozak models. |
| Ribo-seq signal | Total RPF signal, covered nucleotides/codons, and coverage ratio. |
| Periodicity | Frame-specific signal distribution, especially frame-0 enrichment. |
| Start site | Resolved start site, leading-window support, and noncanonical extension. |
| Coverage shape | Uniform, disperse, or skewed RPF distribution across the ORF. |
| Multi-sample support | Reproducibility and support level across multiple Ribo-seq samples. |
| Quantification | Raw P-site density count matrix for reliable smORFs across samples. |


## Notes

- `smorf_scanner` produces many candidate ORFs. Scanner output alone should not be treated as evidence of translation.
- `smorf_cluster` filters and clusters candidates according to sequence features and Kozak context, but it does not use Ribo-seq evidence.
- `smorf_evidence` requires a density list pointing to the P-site density tracks of every sample.
- `smorf_quant` is most useful when multiple Ribo-seq samples or replicates are available.
- Strand-specific P-site density files can be generated with `rpf_Bam2bw` and then supplied to `smorf_evidence` as bedGraph or WIG tracks.
