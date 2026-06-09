# 4.5.1 Check

```bash
rpf_Check -h
```

```text
usage: rpf_Check [-h] -t TRANSCRIPT -b BAM -o OUTPUT [--thread THREAD] [-g {0,1}] [-a {star,hisat2,bowtie2}] [-r] [-l] [-s]
Required arguments:
  -t TRANSCRIPT    transcript annotation in TXT format
  -b BAM           mapping file in BAM format
  -o OUTPUT        output prefix
Options:
  --thread THREAD  number of threads
  -g {0,1}         0 = all reads, 1 = uniquely mapped reads
  -a ALIGNER       star, hisat2, or bowtie2
  -r               count reads aligned to negative strand
  -l               only retain longest transcripts
  -s               calculate read/gene saturation
```

```bash
for bam in ../../3.star/*Aligned.toTranscriptome.out.bam
do
  prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)
  rpf_Check -b $bam -s --thread 10 -t ../../../1.reference/norm/gene.norm.txt -o $prefix_name &>> $prefix_name".log"
done
```

## Merge quality-check results

```bash
merge_length -h
```

```text
usage: merge_length [-h] -l LIST [LIST ...] -o OUTPUT
Required arguments:
  -l LIST     length distribution files, for example '*length_distribution.txt'
  -o OUTPUT   output prefix
```

```bash
merge_saturation -h
```

```text
usage: merge_saturation [-h] -l LIST [LIST ...] -o OUTPUT
Required arguments:
  -l LIST     saturation files, for example '*gene_saturation.txt'
  -o OUTPUT   output prefix
```

```bash
merge_length -l *length_distribution.txt -o RIBO
merge_saturation -l *gene_saturation.txt -o RIBO
```
