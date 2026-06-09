# 4.4 Alignment and quantification

The typical alignment order is:

```text
rRNA → tRNA → ncRNA → mRNA → genome
```

## Bowtie classification

```bash
bowtie -p 12 -v 1 \
  --un sample.norrna.fq \
  --al sample.rrna.fq \
  -x rrna sample.clean.fastq.gz \
  -S sample.rrna.sam
```

## STAR alignment

```bash
STAR \
  --runThreadN 12 \
  --readFilesCommand zcat \
  --genomeDir ./1.reference/star-index \
  --readFilesIn sample.clean.fastq.gz \
  --outSAMtype BAM SortedByCoordinate
```

## RNA-seq quantification

```bash
rsem-calculate-expression \
  -p 10 \
  --no-bam-output \
  --alignments \
  sample.Aligned.toTranscriptome.out.bam \
  ./1.reference/rsem-index/rsem \
  sample
```

## Merge RSEM results

```bash
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt
```
