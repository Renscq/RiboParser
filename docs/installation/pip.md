# 2.2 pip

## Purpose

Install the released Python package into an existing Python environment.

## Requirements

- Linux/POSIX environment
- Python >= 3.12

## Create environment

```bash
conda create -n ribo python=3.12
conda activate ribo
```

## Install

```bash
pip install riboparser
```

## Install external dependencies

pip only installs the RiboParser package itself. The external tools required by the full workflow (e.g. Bowtie, STAR, SAMtools, RSEM) must be installed separately, for example via conda:

```bash
conda install bowtie samtools cutadapt star bedtools subread rsem gffread sra-tools \
  ucsc-genepredtogtf ucsc-gtftogenepred ucsc-gff3togenepred ucsc-bedgraphtobigwig ucsc-bedsort \
  -c bioconda
```

## Update

```bash
pip install -U riboparser
```
