# 4.1 Reference preparation

Reference preparation creates the genome, transcriptome, ncRNA, rRNA, tRNA, STAR, RSEM, and normalized RiboParser annotation files.

## Recommended directory structure

```text
1.reference/
├── genome/
├── mrna/
├── ncrna/
├── norm/
├── rrna/
├── rsem-index/
├── star-index/
└── trna/
```

## Main command

```bash
rpf_Reference \
  -g genome.fa \
  -t annotation.gtf \
  -u 30 \
  -o ./norm/gene
```

## Typical outputs

```text
gene.norm.gtf
gene.norm.txt
gene.norm.rna.fa
gene.norm.cds.fa
```

## External indexes

Common reference indexes include:

```text
Bowtie genome index
Bowtie mRNA index
Bowtie rRNA index
Bowtie tRNA index
Bowtie ncRNA index
STAR genome index
RSEM transcriptome index
```
