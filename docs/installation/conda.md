# 2.2 conda

## Purpose

Install RiboParser and external bioinformatics dependencies in a reproducible environment.

## Create environment

```bash
conda create -n ribo python=3.12
conda activate ribo
```

or:

```bash
micromamba create -n ribo python=3.12
micromamba activate ribo
```

## Install external dependencies

```bash
conda install bowtie samtools cutadapt star bedtools subread rsem gffread sra-tools \
  ucsc-genepredtogtf ucsc-gtftogenepred ucsc-gff3togenepred ucsc-bedgraphtobigwig ucsc-bedsort \
  -c bioconda

conda install pigz -c conda-forge
```

## Install RiboParser

```bash
conda install riboparser -c rensc
```

or:

```bash
micromamba install riboparser -c rensc
```

## Main external tools

| Tool | Use |
|---|---|
| Bowtie | classification against rRNA/tRNA/ncRNA/mRNA/genome indexes |
| SAMtools | BAM sorting, indexing, and filtering |
| Cutadapt | adapter trimming and read filtering |
| STAR | splice-aware genome alignment |
| RSEM | transcript/gene expression quantification |
| gffread | cDNA extraction from genome and annotation |
| SRA Toolkit | public dataset download |
| UCSC utilities | bedGraph/bigWig and annotation conversion |
