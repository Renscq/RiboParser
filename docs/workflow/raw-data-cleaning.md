# 4.3 Raw data cleaning

Raw data cleaning removes adapters, filters short reads, and prepares clean FASTQ files.

## RNA-seq example

```bash
cutadapt --match-read-wildcards \
  -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC \
  -m 25 -O 6 -j 12 \
  -o sample.clean.fastq.gz sample.fastq.gz
```

## Ribo-seq example

```bash
cutadapt --match-read-wildcards \
  -a AAAAAAAA \
  -m 25 -O 6 -j 10 \
  -o sample.clean.fastq.gz sample.fastq.gz
```

## Recommended checks

- read length distribution
- adapter trimming rate
- rRNA/tRNA contamination
- mapping rate
- final usable read count
