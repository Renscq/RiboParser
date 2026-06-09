# 4.1 Reference preparation

Reference preparation creates all reference files required for RNA-seq and Ribo-seq analysis.

## Create directories

```bash
mkdir -p ./sce/1.reference/
cd ./sce/1.reference/

mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index
```

## Download reference files from NCBI

The tutorial uses the yeast reference genome as an example.

```bash
# genome sequence
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# GTF or GFF3
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz

# cDNA sequence
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz

# feature table
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz

# decompression
gunzip *.gz
```

## Generate cDNA sequence from GFF3

```bash
gffread \
  -g GCF_000146045.2_R64_genomic.fna \
  GCF_000146045.2_R64_genomic.gff \
  -F \
  -w cdna.fa
```

## Create genome index using Bowtie

```bash
bowtie-build \
  ../GCF_000146045.2_R64_genomic.fna \
  ./genome/genome \
  --threads 12 \
  &>> ./genome/genome_build.log
```

## Create mRNA index using Bowtie

A helper script is used to extract FASTA sequences by ID.

```bash
retrieve_seq -h
```

```bash
# filter the mRNA sequence
grep -i 'gbkey=mRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./mrna/mrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./mrna/mrna.ids \
  -o ./mrna/mrna.fa \
  &>> ./mrna/mrna_build.log

# build the mRNA index
bowtie-build ./mrna/mrna.fa ./mrna/mrna --threads 12 &>> ./mrna/mrna_build.log
```

## Create rRNA index using Bowtie

```bash
# filter the rRNA sequence
grep -i 'gbkey=rRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./rrna/rrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./rrna/rrna.ids \
  -o ./rrna/rrna.fa \
  &>> ./rrna/rrna_build.log

# build the rRNA index
bowtie-build ./rrna/rrna.fa ./rrna/rrna --threads 12 &>> ./rrna/rrna_build.log
```

## Create tRNA index using Bowtie

```bash
# filter the tRNA sequence
grep -i 'gbkey=tRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./trna/trna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./trna/trna.ids \
  -o ./trna/trna.fa \
  &>> ./trna/trna_build.log

# build the tRNA index
bowtie-build ./trna/trna.fa ./trna/trna --threads 12 &>> ./trna/trna_build.log
```

## Create ncRNA index using Bowtie

```bash
# filter ncRNA sequences
grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' ./cdna.fa \
  | cut -d ' ' -f 1 \
  | cut -c 2- \
  > ./ncrna/ncrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./ncrna/ncrna.ids \
  -o ./ncrna/ncrna.fa \
  &>> ./ncrna/ncrna_build.log

# build the ncRNA index
bowtie-build ./ncrna/ncrna.fa ./ncrna/ncrna --threads 12 &>> ./ncrna/ncrna_build.log
```

## Standardize GTF or GFF3 files with `rpf_Reference`

```bash
rpf_Reference -h
```

```text
usage: rpf_Reference [-h] -g GENOME -t GTF -o OUTPUT [-u UTR] [-c] [-l] [-w]

Required arguments:
  -g GENOME   input genome sequence
  -t GTF      input GTF/GFF3 annotation file
  -o OUTPUT   output prefix

Options:
  -u UTR      add pseudo UTR to leaderless transcripts
  -c          only retain protein-coding transcripts
  -l          only retain the longest protein-coding transcript
  -w          output whole message
```

Run:

```bash
rpf_Reference \
  -g ../GCF_000146045.2_R64_genomic.fna \
  -t ../GCF_000146045.2_R64_genomic.gff \
  -u 30 \
  -o ./norm/gene \
  &>> ./norm/norm_build.log
```

Typical outputs:

```text
gene.norm.gtf
gene.norm.txt
gene.norm.rna.fa
gene.norm.cds.fa
```

## Create genome index using STAR

```bash
STAR \
  --genomeSAindexNbases 11 \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir ./star-index \
  --genomeFastaFiles GCF_000146045.2_R64_genomic.fna \
  --sjdbGTFfile ./norm/gene.norm.gtf
```

## Create transcriptome index using RSEM

```bash
rsem-prepare-reference \
  -p 12 \
  --gtf ../norm/gene.norm.gtf \
  ../GCF_000146045.2_R64_genomic.fna \
  ./rsem-index/rsem
```
