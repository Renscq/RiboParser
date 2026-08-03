# 4.1 Reference preparation

## Purpose

Reference preparation builds all files needed for read classification, alignment, quantification, RPF density construction, and codon-level Ribo-seq analysis.

The workflow is divided into five steps:

```text
Step 1: Prepare directories and reference files
Step 2: Build Bowtie indexes (genome, mRNA, rRNA, tRNA, ncRNA)
Step 3: Build the normalized RiboParser reference (rpf_Reference)
Step 4: Build the STAR index
Step 5: Build the RSEM index
```

All commands are run inside the `1.reference` directory created in Step 1.

## Inputs

| Input | Description | Required |
|---|---|---|
| genome FASTA | genomic sequence | yes |
| GTF/GFF3 annotation | gene, transcript, exon, CDS, and UTR annotation | yes |
| feature table | optional NCBI feature table | no |
| external tools | Bowtie, STAR, RSEM, gffread | yes |

## Step 1: Prepare directories and reference files

### 1.1 Create directories

```bash
mkdir -p ./sce/1.reference/
cd ./sce/1.reference/

mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index star-index
```

### 1.2 Download reference files

Download the genome sequence, annotations, and optional feature table from NCBI. Use the URLs corresponding to your species.

```bash
# genome sequence
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
# GTF or GFF3 annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
# optional feature table
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz

gunzip *.gz
```

## Step 2: Build Bowtie indexes

### 2.1 Generate cDNA from the genome and GFF3

The transcript FASTA used for mRNA/rRNA/tRNA/ncRNA extraction is generated from the genome and GFF3 with `gffread`, not downloaded directly.

```bash
gffread \
  -g ./GCF_000146045.2_R64_genomic.fna \
  ./GCF_000146045.2_R64_genomic.gff \
  -F \
  -w ./cdna.fa
```

### 2.2 Build the genome index

```bash
bowtie-build \
  ./GCF_000146045.2_R64_genomic.fna \
  ./genome/genome \
  --threads 12 \
  &>> ./genome/genome_build.log
```

### 2.3 Build the mRNA index

```bash
grep -i 'gbkey=mRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./mrna/mrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./mrna/mrna.ids \
  -o ./mrna/mrna.fa \
  &>> ./mrna/mrna_build.log

bowtie-build ./mrna/mrna.fa ./mrna/mrna --threads 12 &>> ./mrna/mrna_build.log
```

### 2.4 Build the rRNA index

```bash
grep -i 'gbkey=rRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./rrna/rrna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./rrna/rrna.ids \
  -o ./rrna/rrna.fa \
  &>> ./rrna/rrna_build.log

bowtie-build ./rrna/rrna.fa ./rrna/rrna --threads 12 &>> ./rrna/rrna_build.log
```

### 2.5 Build the tRNA index

```bash
grep -i 'gbkey=tRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./trna/trna.ids

retrieve_seq \
  -i ./cdna.fa \
  -n ./trna/trna.ids \
  -o ./trna/trna.fa \
  &>> ./trna/trna_build.log

bowtie-build ./trna/trna.fa ./trna/trna --threads 12 &>> ./trna/trna_build.log
```

### 2.6 Build the ncRNA index

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

## Step 3: Build the normalized RiboParser reference

`rpf_Reference` normalizes the GTF/GFF3 annotation into RiboParser-compatible transcript tables and sequences, which are required by most downstream modules (e.g. quality control, gene-level analysis, and codon-level analysis).

### 3.1 Command

```bash
rpf_Reference \
  -g ./GCF_000146045.2_R64_genomic.fna \
  -t ./GCF_000146045.2_R64_genomic.gff \
  -u 30 \
  -o ./norm/gene \
  &>> ./norm/norm_build.log
```

### 3.2 Required parameters

| Parameter | Meaning |
|---|---|
| `-g`, `--genome` | genome FASTA file |
| `-t`, `--gtf` | GTF or GFF3 annotation file |
| `-o`, `--output` | output prefix |

### 3.3 Optional parameters

| Parameter | Meaning |
|---|---|
| `-u`, `--utr` | add pseudo UTR length for leaderless transcripts |
| `-c` | retain only protein-coding transcripts |
| `-l` | retain only the longest protein-coding transcript per gene |
| `-w` | output full message table |

### 3.4 Outputs

| File | Description |
|---|---|
| `gene.norm.gtf` | normalized GTF |
| `gene.norm.txt` | RiboParser transcript annotation table |
| `gene.norm.rna.fa` | normalized transcript FASTA |
| `gene.norm.cds.fa` | normalized CDS FASTA |

## Step 4: Build the STAR index

The STAR index is built from the genome sequence and the normalized GTF, and is used for RNA-seq alignment.

```bash
STAR \
  --genomeSAindexNbases 11 \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir ./star-index \
  --genomeFastaFiles ./GCF_000146045.2_R64_genomic.fna \
  --sjdbGTFfile ./norm/gene.norm.gtf
```

## Step 5: Build the RSEM index

The RSEM index is used for transcript-level quantification.

```bash
rsem-prepare-reference \
  -p 12 \
  --gtf ./norm/gene.norm.gtf \
  ./GCF_000146045.2_R64_genomic.fna \
  ./rsem-index/rsem
```

## Quality checklist

- Bowtie indexes exist for genome, mRNA, rRNA, tRNA, and ncRNA.
- STAR index directory is complete.
- RSEM index files are complete.
- `gene.norm.txt`, `gene.norm.rna.fa`, and `gene.norm.cds.fa` are generated.
- Chromosome names are consistent across genome, GTF/GFF3, BAM, and normalized reference.
