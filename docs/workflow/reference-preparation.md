# 4.1 Reference preparation

## Create directories

```bash
mkdir -p ./sce/1.reference/
cd ./sce/1.reference/
mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index star-index
```

## Download reference files from NCBI

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz
gunzip *.gz
```

## Generate cDNA sequence from GFF3

```bash
gffread -g GCF_000146045.2_R64_genomic.fna GCF_000146045.2_R64_genomic.gff -F -w cdna.fa
```

## Build Bowtie indexes

```bash
bowtie-build ../GCF_000146045.2_R64_genomic.fna ./genome/genome --threads 12 &>> ./genome/genome_build.log
```

```bash
grep -i 'gbkey=mRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./mrna/mrna.ids
retrieve_seq -i ./cdna.fa -n ./mrna/mrna.ids -o ./mrna/mrna.fa &>> ./mrna/mrna_build.log
bowtie-build ./mrna/mrna.fa ./mrna/mrna --threads 12 &>> ./mrna/mrna_build.log
```

```bash
grep -i 'gbkey=rRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./rrna/rrna.ids
retrieve_seq -i ./cdna.fa -n ./rrna/rrna.ids -o ./rrna/rrna.fa &>> ./rrna/rrna_build.log
bowtie-build ./rrna/rrna.fa ./rrna/rrna --threads 12 &>> ./rrna/rrna_build.log
```

```bash
grep -i 'gbkey=tRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./trna/trna.ids
retrieve_seq -i ./cdna.fa -n ./trna/trna.ids -o ./trna/trna.fa &>> ./trna/trna_build.log
bowtie-build ./trna/trna.fa ./trna/trna --threads 12 &>> ./trna/trna_build.log
```

```bash
grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./ncrna/ncrna.ids
retrieve_seq -i ./cdna.fa -n ./ncrna/ncrna.ids -o ./ncrna/ncrna.fa &>> ./ncrna/ncrna_build.log
bowtie-build ./ncrna/ncrna.fa ./ncrna/ncrna --threads 12 &>> ./ncrna/ncrna_build.log
```

## Normalize annotation

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
  -l          only retain the longest protein-coding transcript for each gene
  -w          output whole message
```

```bash
rpf_Reference -g ../GCF_000146045.2_R64_genomic.fna -t ../GCF_000146045.2_R64_genomic.gff -u 30 -o ./norm/gene &>> ./norm/norm_build.log
```

## STAR and RSEM indexes

```bash
STAR --genomeSAindexNbases 11 --runThreadN 12 --runMode genomeGenerate --genomeDir ./star-index --genomeFastaFiles GCF_000146045.2_R64_genomic.fna --sjdbGTFfile ./norm/gene.norm.gtf
```

```bash
rsem-prepare-reference -p 12 --gtf ../norm/gene.norm.gtf ../GCF_000146045.2_R64_genomic.fna ./rsem-index/rsem
```
