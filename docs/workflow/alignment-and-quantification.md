# 4.4 Alignment and quantification

## Bowtie classification

```bash
rrna='../../sce/1.reference/rrna/rrna'
trna='../../sce/1.reference/trna/trna'
ncrna='../../sce/1.reference/ncrna/ncrna'
mrna='../../sce/1.reference/mrna/mrna'
chrom='../../sce/1.reference/genome/genome'
threads=12
mismatch=1

for fq in ../1.cleandata/*fastq.gz
do
  fqname=$(basename $fq .fastq.gz)
  bowtie -p $threads -v $mismatch --un="$fqname".norrna.fq --al="$fqname".rrna.fq -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log
  bowtie -p $threads -v $mismatch --un="$fqname".notrna.fq --al="$fqname".trna.fq -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log
  bowtie -p $threads -v $mismatch --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log
  bowtie -p $threads -v $mismatch --un="$fqname".nomrna.fq --al="$fqname".mrna.fq -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log
  bowtie -p $threads -v $mismatch --un="$fqname".nogenome.fq --al="$fqname".genome.fq -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log
  pigz *fq
  for sam in *.sam; do samtools view -h -F 4 $sam | samtools sort -@ $threads -o $(basename $sam sam)bam; rm $sam; done
done
```

## Merge Bowtie logs

```bash
merge_bwt_log -h
```

```text
usage: merge_bwt_log [-h] -l LIST [LIST ...] -o OUTPUT [-n NAME]
Required arguments:
  -l LIST    Bowtie mapping log files
  -o OUTPUT  output prefix
Options:
  -n NAME    comma-separated database names; example: rRNA,tRNA,ncRNA,mRNA,Genome
```

```bash
merge_bwt_log -n rRNA,tRNA,ncRNA,mRNA,Genome -l *log -o RNA_seq &>> merge_bowtie.log
```

## STAR alignment and RSEM quantification

```bash
STAR --runThreadN 12 --readFilesCommand zcat --genomeDir ../../1.reference/star-index/ --readFilesIn sample.noncrna.fq.gz --outFileNamePrefix sample --outSAMtype BAM Unsorted --outFilterType BySJout --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --outSAMattributes All --alignEndsType Local --outFilterMultimapNmax 3 --outFilterMismatchNmax 1 --alignIntronMax 10000 --outFilterMatchNmin 20
```

```bash
rsem-calculate-expression -p 10 --no-bam-output --alignments -q sample.Aligned.toTranscriptome.out.bam ../../1.reference/rsem-index/rsem sample
```

## Merge RSEM results

```bash
merge_rsem -h
```

```text
usage: merge_rsem [-h] -l LIST [LIST ...] -o OUTPUT [-c {expected_count,TPM,FPKM}]
```

```bash
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt
```
