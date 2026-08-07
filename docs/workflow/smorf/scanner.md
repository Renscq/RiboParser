# 4.8.1 smORF scanner

## Function

`smorf_scanner` scans transcript-centric candidate small open reading frames (smORFs) from a genome FASTA file and a genePred annotation file.

It reconstructs spliced transcript sequences from the genome, searches candidate ORFs using user-defined start codons and ORF length limits, classifies ORF categories by their position relative to the annotated CDS, optionally marks nested and overlapping ORFs, and writes four output files: an ORF annotation table, a full ORF information table, and the nucleotide and peptide FASTA files.

The scanner is designed for large genomes: each reading frame is translated exactly once and candidate ORFs reuse slices of the cached frame peptide, so nested ORFs are not translated independently. Multiprocessing (`-t`) splits the transcript set across workers.

## Workflow

```text
genome FASTA + genePred annotation → smorf_scanner → smorf.raw.genePred / message.txt / nt.fa / pep.fa → smorf_cluster
```

The command runs in two steps:

1. **Checking the input arguments** – validate the genome and annotation files and the ORF scanning parameters.
2. **Scanning ORFs from transcript sequences** – reconstruct transcript sequences, scan the reading frames, classify categories, and write all output files. With `-t 1` a single-process streaming path is used; with `-t > 1` transcripts are processed in parallel and the output is written in balanced partitions.

## Input files

| Input | Required | Description |
|---|---|---|
| Genome FASTA | Yes | Genome sequence used to reconstruct spliced transcript sequences. |
| genePred annotation | Yes | Transcript annotation in genePred format. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `-g`, `--genome` | Yes | Input genome sequence file in FASTA or FASTA.GZ format. |
| `-a`, `--annotation` | Yes | Input transcript annotation file in genePred format. |
| `-o`, `--out-prefix` | No | Prefix of all output files. Default: `ORF`. |
| `-p`, `--orf-prefix` | No | Prefix used to generate stable ORF identifiers. Default: `ORF`. |
| `-s`, `--start-codons` | No | Comma-separated start codons used for ORF scanning. RNA `U` is converted to `T`. Default: `ATG`. |
| `-m`, `--min-aa` | No | Minimum ORF peptide length in amino acids. Default: `8`. |
| `-M`, `--max-aa` | No | Maximum ORF peptide length in amino acids. Default: `10000`. |
| `-x`, `--scan-strand` | No | Strand mode for scanning: `sense`, `antisense`, or `both`. Default: `sense`. |
| `-u`, `--kozak-up` | No | Number of upstream nucleotides extracted for the Kozak context. Default: `6`. |
| `-d`, `--kozak-down` | No | Number of downstream nucleotides after the start codon extracted for the Kozak context. Default: `6`. |
| `-I`, `--include-stop` | No | Keep the stop-codon symbol in peptide sequences. Default: `False`. |
| `-O`, `--mark-overlap` | No | Mark nested or overlapping ORFs after scanning. Default: `False`. |
| `-R`, `--remove-discarded` | No | Remove same-frame internal ORFs that the classifier labels as discarded. Default: `False`. |
| `-t`, `--thread` | No | Number of worker processes. Default: `1`. |

## Output files

Assuming `-o gmx4`, the command writes:

| Output | Description |
|---|---|
| `gmx4.genePred` | ORF annotation in a genePredExt-like format. |
| `gmx4.message.txt` | Full ORF information table. This is the main input for `smorf_cluster`. |
| `gmx4.nt.fa` | Candidate ORF nucleotide sequences in FASTA format. |
| `gmx4.pep.fa` | Candidate ORF peptide sequences in FASTA format. |

### ORF information table (`gmx4.message.txt`)

One row per candidate ORF, tab-separated:

```text
orf_id  gene_id  transcript_id  chrom  strand  source_strand  category  priority  overlap_type
frame  tx_orf_start  tx_orf_end  genomic_start  genomic_end  start_codon  stop_codon
nt_length  aa_length  kozak_seq  completeness  exon_count  exon_starts  exon_ends
kozak_start_index  ambiguous_codon_count
```

| Column | Description |
|---|---|
| `orf_id` | Generated stable ORF identifier (for example `ORF00000001`). |
| `gene_id` / `transcript_id` | Source gene and transcript of the ORF. |
| `chrom` / `strand` | Genomic location and strand of the ORF. |
| `source_strand` | Scanning orientation relative to the transcript: `sense` or `antisense`. |
| `category` | Positional ORF category assigned by the classifier (see below). |
| `priority` | Overlap-filtering priority: `primary`, `secondary`, or `discarded`. |
| `overlap_type` | Relationship to overlapping ORFs (only filled with `-O`). |
| `frame` | Reading frame in the scanned sequence. |
| `tx_orf_start` / `tx_orf_end` | ORF boundaries in transcript coordinates. |
| `genomic_start` / `genomic_end` | Minimum/maximum genomic coordinate covered by the ORF. |
| `start_codon` / `stop_codon` | Detected start codon; stop codon, or `NA` for a partial ORF. |
| `nt_length` / `aa_length` | ORF nucleotide length (including the terminal stop codon when present) and peptide length. |
| `kozak_seq` | Fixed-width, `N`-padded start-codon context. |
| `completeness` | `complete` or `3prime_partial`. |
| `exon_count` / `exon_starts` / `exon_ends` | ORF exon structure in ascending genomic order. |
| `kozak_start_index` | Zero-based start-codon index in `kozak_seq`. |
| `ambiguous_codon_count` | Number of ORF codons containing non-ACGT characters. |

### Category classification

ORFs are classified by their position relative to the annotated CDS of the transcript:

| Category | Meaning |
|---|---|
| `annotated_ORF` | Matches the annotated CDS (highest priority). |
| `uORF` | Fully upstream of the CDS (5' UTR). |
| `dORF` | Fully downstream of the CDS (3' UTR). |
| `emORF` | Contains the CDS start or end (extended / merged ORF). |
| `same_frame_iORF` | Internal to the CDS, same reading frame (marked `discarded` by default). |
| `iORF` | Internal to the CDS, different reading frame. |
| `overlap_uORF` | Overlaps the CDS start. |
| `overlap_dORF` | Overlaps the CDS end. |
| `lncORF` | Found on a non-coding transcript. |
| `antisense_ORF` | Found on the antisense strand. |
| `other_ORF` | Any other configuration. |

### Sequence FASTA files

The FASTA headers carry the ORF metadata, for example:

```text
>ORF00000001 gene=GlmaCp001 transcript=GlmaCp001 type=annotated_ORF strand=- length=1062
ATGACTGCAATTTTAGAGAGACGCGAGAGCGAAAGCCTATGGGGTCGCTTCTGTAACTGG
```

- `gmx4.nt.fa` – nucleotide sequences in coding orientation; the header includes `length=`.
- `gmx4.pep.fa` – translated peptide sequences; the header includes `aa_length=`. The terminal stop symbol is kept only with `-I`.

## Examples

### Use the common scanning setup on a large genome

```bash
cd ./sce/5.smorf/01.scanner

smorf_scanner \
  --genome ~/gmx/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna \
  --annotation ~/gmx/norm/gmx4.genepred \
  --out-prefix gmx4 \
  --start-codons ATG,CTG,GTG,TTG,ACG,ATA,ATT,ATC \
  --min-aa 8 \
  --max-aa 10000 \
  --scan-strand sense \
  --kozak-up 6 \
  --kozak-down 6 \
  --mark-overlap \
  --remove-discarded \
  --thread 20 \
  &> gmx4.scan.log
```

This example scans the soybean genome (93,168 transcripts). With 20 workers it finished in about a minute and wrote 13,347,430 candidate ORFs into `gmx4.genePred`, `gmx4.message.txt`, `gmx4.nt.fa`, and `gmx4.pep.fa` (approximately 1.5–2.8 GB each).

## Notes

- `--scan-strand sense` scans ORFs on the annotated transcript strand only; `--scan-strand antisense` scans the reverse-complemented transcript sequences; `--scan-strand both` covers both but substantially increases the number of candidates.
- `-O/--mark-overlap` is recommended when downstream filtering should distinguish primary, nested, or overlapping candidates; without it `overlap_type` is left at `none`.
- `-R/--remove-discarded` removes complete same-frame internal ORFs that the classifier labels as `discarded`; keep it enabled for cleaner candidate sets.
- The scanner only generates candidate ORFs: scanner output alone should not be treated as evidence of translation.
- Scanner output is a candidate set. Use `smorf_cluster` (4.8.2) and `smorf_evidence` (4.8.3) before prioritizing translated smORFs.
