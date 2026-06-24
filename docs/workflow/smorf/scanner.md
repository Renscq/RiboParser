# 4.8.1 smORF scanner

## Function

`smorf_scanner` scans transcript-centric candidate ORFs from a genome FASTA file and a genePred annotation file.

It reconstructs spliced transcript sequences, searches candidate ORFs using user-defined start codons and ORF length limits, classifies ORF categories, optionally marks overlapping/nested ORFs, and writes ORF annotation and sequence files.

## Input files

| Input | Required | Description |
|---|---|---|
| Genome FASTA | Yes | Genome sequence used to reconstruct transcript sequences. |
| genePred annotation | Yes | Transcript annotation in genePred format. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `-g`, `--genome` | Yes | Input genome FASTA file. |
| `-a`, `--annotation` | Yes | Input transcript annotation in genePred format. |
| `-o`, `--out-prefix` | No | Prefix of output files. Default: `ORF`. |
| `-p`, `--orf-prefix` | No | Prefix used to generate stable ORF IDs. Default: `ORF`. |
| `-s`, `--start-codons` | No | Comma-separated start codons used for ORF scanning. Default: `ATG`. |
| `-m`, `--min-aa` | No | Minimum ORF peptide length in amino acids. Default: `8`. |
| `-M`, `--max-aa` | No | Maximum ORF peptide length in amino acids. Default: `10000`. |
| `-x`, `--scan-strand` | No | Strand mode for scanning: `sense`, `antisense`, or `both`. Default: `sense`. |
| `-u`, `--kozak-up` | No | Number of upstream nucleotides extracted for Kozak context. Default: `6`. |
| `-d`, `--kozak-down` | No | Number of downstream nucleotides after the start codon extracted for Kozak context. Default: `6`. |
| `-t`, `--threads` | No | Number of worker processes used for ORF scanning. Default: `1`. |
| `-O`, `--mark-overlap` | No | Mark nested or overlapping ORFs after scanning. Default: `False`. |
| `-R`, `--remove-discarded` | No | Remove same-frame internal ORFs marked as discarded. Default: `False`. |
| `-I`, `--include-stop` | No | Keep the stop-codon symbol in peptide sequences. Default: `False`. |

## Output files

Assuming `-o smorf.raw`, the command writes:

| Output | Description |
|---|---|
| `smorf.raw.genePred` | ORF annotation in a genePredExt-like format. |
| `smorf.raw.message.txt` | Full ORF information table. This is the main input for `smorf_filter`. |
| `smorf.raw.nt.fa` | Candidate ORF nucleotide sequences in FASTA format. |
| `smorf.raw.pep.fa` | Candidate ORF peptide sequences in FASTA format. |

## Example

```bash
smorf_scanner \
  -g genome.fa \
  -a annotation.genePred \
  -o smorf.raw \
  -p smORF \
  -s ATG,CTG,GTG,TTG \
  -m 8 \
  -M 1000 \
  -x both \
  -u 6 \
  -d 6 \
  -O \
  -t 8
```

## Notes

- `--scan-strand sense` scans ORFs on the annotated transcript strand only.
- `--scan-strand antisense` scans antisense transcript sequences.
- `--scan-strand both` is useful when antisense ORFs are also of interest, but it may substantially increase the number of candidates.
- `--mark-overlap` is recommended when downstream filtering should distinguish primary, nested, or overlapping candidates.
- Scanner output is a candidate set. Use `smorf_filter` and `smorf_evidence` before prioritizing translated smORFs.
