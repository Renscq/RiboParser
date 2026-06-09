# 4.3 Raw data cleaning

```bash
for fq in ../../2.rawdata/rna-seq/*fastq.gz
do
  cutadapt --match-read-wildcards -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC -m 25 -O 6 -j 12 -o $(basename $fq fastq.gz)clean.fastq.gz $fq &>> $fq".log"
done
```

```bash
for fq in ../../2.rawdata/ribo-seq/*fastq.gz
do
  cutadapt --match-read-wildcards -a AAAAAAAA -m 25 -O 6 -j 10 -o $(basename $fq fastq.gz)clean.fastq.gz $fq &>> $fq".log"
done
```
