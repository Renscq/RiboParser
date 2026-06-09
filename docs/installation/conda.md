# 2.2 conda

Create an environment:

```bash
conda create -n ribo python=3.12
conda activate ribo
```

Install RiboParser:

```bash
conda install riboparser -c rensc
```

Or use micromamba:

```bash
micromamba create -n ribo python=3.12
micromamba activate ribo
micromamba install riboparser -c rensc
```

## Common external dependencies

```bash
conda install bowtie samtools cutadapt star bedtools subread rsem gffread sra-tools \
  ucsc-genepredtogtf ucsc-gtftogenepred ucsc-gff3togenepred ucsc-bedgraphtobigwig ucsc-bedsort \
  -c bioconda

conda install pigz -c conda-forge
```
