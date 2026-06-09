# 2.2 conda

First install `miniconda` or `micromamba`.

## Create environment

```bash
conda create -n ribo
conda activate ribo
```

## Install RiboParser directly

```bash
# conda
conda install riboparser -c rensc

# micromamba
micromamba install riboparser -c rensc
```

## Install software dependencies step by step

```bash
conda install bowtie samtools cutadapt star bedtools subread rsem gffread sra-tools \
  ucsc-genepredtogtf ucsc-gtftogenepred ucsc-gff3togenepred ucsc-bedgraphtobigwig ucsc-bedsort \
  -c bioconda

conda install pigz -c conda-forge
```

Then install RiboParser:

```bash
pip install riboparser
```
