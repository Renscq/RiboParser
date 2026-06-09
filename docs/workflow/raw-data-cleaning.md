# 4.3 Raw data cleaning

## Purpose

Remove adapters, discard very short reads, and prepare clean FASTQ files for read classification and alignment.

## RNA-seq cleaning

```bash
mkdir -p ./sce/3.rna-seq/1.cleandata/
cd ./sce/3.rna-seq/1.cleandata/

for fq in ../../2.rawdata/rna-seq/*fastq.gz
do
  cutadapt --match-read-wildcards \
    -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC \
    -m 25 \
    -O 6 \
    -j 12 \
    -o $(basename $fq fastq.gz)clean.fastq.gz \
    $fq \
    &>> $(basename $fq .fastq.gz)".cutadapt.log"
done
```

## Ribo-seq cleaning

```bash
mkdir -p ./sce/4.ribo-seq/1.cleandata/
cd ./sce/4.ribo-seq/1.cleandata/

for fq in ../../2.rawdata/ribo-seq/*fastq.gz
do
  cutadapt --match-read-wildcards \
    -a AAAAAAAA \
    -m 25 \
    -O 6 \
    -j 10 \
    -o $(basename $fq fastq.gz)clean.fastq.gz \
    $fq \
    &>> $(basename $fq .fastq.gz)".cutadapt.log"
done
```

## Parameter explanation

| Parameter | Meaning |
|---|---|
| `--match-read-wildcards` | allow wildcard matching in reads |
| `-a` | adapter sequence |
| `-m` | minimum retained read length |
| `-O` | minimum adapter overlap length |
| `-j` | threads |
| `-o` | output FASTQ |

## Outputs

| File | Description |
|---|---|
| `*.clean.fastq.gz` | cleaned FASTQ |
| `*.cutadapt.log` | adapter trimming statistics |

## QC interpretation

Check:

- percentage of reads with adapters
- retained read number
- post-trimming read length distribution
- Ribo-seq dominant RPF length range
- RNA-seq read length consistency
