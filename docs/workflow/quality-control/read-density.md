# 4.5.4 Read density

## Purpose

Generate transcript-level RNA-seq or Ribo-seq density files. Density files are the central input for merging, periodicity, metaplot, coverage, correlation, gene-level, and codon-level analyses.

## `rna_Density`

### Function

Generate RNA-seq transcript-level read density using a constant RNA offset table.

### Input files

| Input | Description |
|---|---|
| BAM | filtered BAM from `rpf_Check` |
| offset table | `sample_offset.txt` from `rna_Offset` |
| transcript annotation | `gene.norm.txt` |
| transcript FASTA | `gene.norm.rna.fa` |

### Parameters

| Parameter | Meaning |
|---|---|
| `-b / --bam` | input BAM |
| `-p / --psite` | RNA offset table |
| `-t / --transcript` | normalized transcript annotation |
| `-s / --sequence` | transcript FASTA |
| `-o / --output` | output prefix |
| `-m / --min` | minimum read length |
| `-M / --max` | maximum read length |
| `-l` | retain longest transcript per gene |
| `--thread` | number of worker threads |
| `--silence` | suppress verbose output |

### Example

```bash
for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rna_Density \
    -b $bam \
    -m 27 \
    -M 50 \
    -l \
    --thread 10 \
    -p ../03.offset/$prefix_name"_offset.txt" \
    -s ../../../1.reference/norm/gene.norm.rna.fa \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

### Output

| Output | Description |
|---|---|
| `sample_rna.txt` | RNA-seq transcript density |
| `sample.log` | running log |

## `rpf_Density`

### Function

Generate Ribo-seq P-site density after offset correction and optional periodicity filtering.

### Input files

| Input | Description |
|---|---|
| BAM | filtered Ribo-seq BAM |
| P-site offset table | usually `sample_SSCBM_offset.txt` |
| transcript annotation | `gene.norm.txt` |
| transcript FASTA | `gene.norm.rna.fa` |

### Parameters

| Parameter | Meaning |
|---|---|
| `-b / --bam` | input BAM |
| `-p / --psite` | P-site offset table |
| `-t / --transcript` | normalized transcript annotation |
| `-s / --sequence` | transcript FASTA |
| `-o / --output` | output prefix |
| `-m / --min` | minimum read length |
| `-M / --max` | maximum read length |
| `--period` | minimum 3-nt periodicity threshold for retained reads/transcripts |
| `-l` | retain longest transcript per gene |
| `--thread` | number of worker threads |
| `--silence` | suppress verbose output |

### Example

```bash
for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rpf_Density \
    -b $bam \
    -m 27 \
    -M 33 \
    --period 40 \
    -l \
    --thread 12 \
    -p ../03.offset/$prefix_name"_SSCBM_offset.txt" \
    -s ../../../1.reference/norm/gene.norm.rna.fa \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

### Output

| Output | Description |
|---|---|
| `sample_rpf.txt` | Ribo-seq P-site density |
| `sample.log` | running log |

## Interpretation

RNA-seq density represents transcript coverage, whereas Ribo-seq density represents ribosome-protected fragment positions after offset correction. Ribo-seq density should be frame-aware and should only be used for codon-level analysis after periodicity validation.
