# 4.1 Reference preparation

## Purpose

Reference preparation builds all files needed for read classification, alignment, quantification, RPF density construction, and codon-level Ribo-seq analysis.

## Inputs

| Input | Description |
|---|---|
| genome FASTA | genomic sequence |
| GTF/GFF3 annotation | gene, transcript, exon, CDS, and UTR annotation |
| cDNA/transcript FASTA | transcript sequence for rRNA/tRNA/ncRNA/mRNA extraction |
| feature table | optional NCBI feature table |
| external tools | Bowtie, STAR, RSEM, gffread |

## Create directories

```bash
mkdir -p ./sce/1.reference/
cd ./sce/1.reference/

mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index star-index
```

## Download reference files

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz

gunzip *.gz
```

## Generate cDNA from genome and GFF3

```bash
gffread \
  -g GCF_000146045.2_R64_genomic.fna \
  GCF_000146045.2_R64_genomic.gff \
  -F \
  -w cdna.fa
```

## Build Bowtie genome index

```bash
bowtie-build \
  ../GCF_000146045.2_R64_genomic.fna \
  ./genome/genome \
  --threads 12 \
  &>> ./genome/genome_build.log
```

## Build mRNA index

```bash
grep -i 'gbkey=mRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./mrna/mrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./mrna/mrna.ids \
  -o ./mrna/mrna.fa \
  &>> ./mrna/mrna_build.log

bowtie-build ./mrna/mrna.fa ./mrna/mrna --threads 12 &>> ./mrna/mrna_build.log
```

## Build rRNA index

```bash
grep -i 'gbkey=rRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./rrna/rrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./rrna/rrna.ids \
  -o ./rrna/rrna.fa \
  &>> ./rrna/rrna_build.log

bowtie-build ./rrna/rrna.fa ./rrna/rrna --threads 12 &>> ./rrna/rrna_build.log
```

## Build tRNA index

```bash
grep -i 'gbkey=tRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./trna/trna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./trna/trna.ids \
  -o ./trna/trna.fa \
  &>> ./trna/trna_build.log

bowtie-build ./trna/trna.fa ./trna/trna --threads 12 &>> ./trna/trna_build.log
```

## Build ncRNA index

```bash
grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' ./cdna.fa \
  | cut -d ' ' -f 1 \
  | cut -c 2- \
  > ./ncrna/ncrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./ncrna/ncrna.ids \
  -o ./ncrna/ncrna.fa \
  &>> ./ncrna/ncrna_build.log

bowtie-build ./ncrna/ncrna.fa ./ncrna/ncrna --threads 12 &>> ./ncrna/ncrna_build.log
```

## `rpf_Reference` full help

### Function

Normalize annotation into RiboParser-compatible transcript tables and sequences.

### Command

```bash
rpf_Reference \
  -g ../GCF_000146045.2_R64_genomic.fna \
  -t ../GCF_000146045.2_R64_genomic.gff \
  -u 30 \
  -o ./norm/gene \
  &>> ./norm/norm_build.log
```

### Required parameters

| Parameter | Meaning |
|---|---|
| `-g`, `--genome` | genome FASTA file |
| `-t`, `--gtf` | GTF or GFF3 annotation file |
| `-o`, `--output` | output prefix |

### Optional parameters

| Parameter | Meaning |
|---|---|
| `-u`, `--utr` | add pseudo UTR length for leaderless transcripts |
| `-c` | retain only protein-coding transcripts |
| `-l` | retain only the longest protein-coding transcript per gene |
| `-w` | output full message table |

### Outputs

| File | Description |
|---|---|
| `gene.norm.gtf` | normalized GTF |
| `gene.norm.txt` | RiboParser transcript annotation table |
| `gene.norm.rna.fa` | normalized transcript FASTA |
| `gene.norm.cds.fa` | normalized CDS FASTA |

## Build STAR index

```bash
STAR \
  --genomeSAindexNbases 11 \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir ./star-index \
  --genomeFastaFiles GCF_000146045.2_R64_genomic.fna \
  --sjdbGTFfile ./norm/gene.norm.gtf
```

## Build RSEM index

```bash
rsem-prepare-reference \
  -p 12 \
  --gtf ../norm/gene.norm.gtf \
  ../GCF_000146045.2_R64_genomic.fna \
  ./rsem-index/rsem
```

## Quality checklist

- Bowtie indexes exist for genome, mRNA, rRNA, tRNA, and ncRNA.
- STAR index directory is complete.
- RSEM index files are complete.
- `gene.norm.txt`, `gene.norm.rna.fa`, and `gene.norm.cds.fa` are generated.
- Chromosome names are consistent across genome, GTF/GFF3, BAM, and normalized reference.
