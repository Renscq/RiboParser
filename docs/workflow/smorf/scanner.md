# 4.8.1 smORF scanner

## Purpose

`smorf_scanner` scans candidate ORFs from genome and genePred annotation. It can scan sense, antisense, or both strands and classify overlapping ORFs.

## Input files

| Input | Description |
|---|---|
| genome FASTA | genome sequence |
| genePred annotation | transcript annotation in genePred format |
| start codon list | allowed initiation codons |
| length range | minimum and maximum ORF amino-acid length |

## Parameters

| Parameter | Meaning |
|---|---|
| `-g / --genome` | input genome FASTA |
| `-a / --annotation` | input genePred annotation |
| `-o / --out-prefix` | output prefix |
| `--orf-prefix` | prefix added to ORF IDs |
| `--start-codons` | comma-separated start codons, such as `ATG,CTG,GTG,TTG` |
| `--min-aa` | minimum ORF length in amino acids |
| `--max-aa` | maximum ORF length in amino acids |
| `--scan-strand` | `sense`, `antisense`, or `both` |
| `--kozak-up` | upstream nucleotides included for Kozak context |
| `--kozak-down` | downstream nucleotides after start codon |
| `-t / --threads` | number of worker processes |
| `--mark-overlap` | mark nested or overlapping ORFs |
| `--remove-discarded` | remove same-frame internal ORFs or discarded records |
| `--include-stop` | keep stop codon symbol in peptide sequence |

## Example

```bash
smorf_scanner \
  --genome ../genome/GCF_mine_genomic.fna \
  --annotation ../norm/mine.genepred \
  --out-prefix mine \
  --start-codons ATG \
  --min-aa 8 \
  --max-aa 10000 \
  --scan-strand both \
  --kozak-up 6 \
  --kozak-down 6 \
  --mark-overlap \
  --threads 20
```

## Output files

| Output | Description |
|---|---|
| `mine.message.txt` | complete scanned ORF information |
| `mine.sequence.fa` | candidate ORF nucleotide sequences |
| `mine.peptide.fa` | candidate ORF peptide sequences |
| `log file` | running log |

## Result interpretation

Scanner output contains many candidates. Do not treat scanner-only ORFs as translated ORFs; use filtering and Ribo-seq evidence evaluation.
