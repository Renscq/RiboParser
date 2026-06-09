# 4.3 Raw data cleaning

Raw data cleaning removes adapters and short reads. The `GSE67387` data are already cleaned, so the commands below are provided as general examples.

## RNA-seq data cleaning

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
    &>> $fq".log"
done
```

## Ribo-seq data cleaning

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
    &>> $fq".log"
done
```

## Recommended checks after cleaning

- adapter removal rate
- read length distribution
- retained read count
- Ribo-seq read length enrichment
- rRNA/tRNA/ncRNA contamination after classification
