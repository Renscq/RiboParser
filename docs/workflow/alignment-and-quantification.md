# 4.4 Alignment and quantification

## Purpose

Classify reads by reference type, align remaining reads to the genome, and quantify transcript/gene abundance.

## Read classification order

```text
rRNA → tRNA → ncRNA → mRNA → genome
```

This order helps separate contaminant reads, transcript-derived reads, and genome-mapped reads.

## Bowtie classification template

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

  bowtie -p $threads -v $mismatch --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
    -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

  bowtie -p $threads -v $mismatch --un="$fqname".notrna.fq --al="$fqname".trna.fq \
    -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

  bowtie -p $threads -v $mismatch --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
    -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

  bowtie -p $threads -v $mismatch --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
    -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

  bowtie -p $threads -v $mismatch --un="$fqname".nogenome.fq --al="$fqname".genome.fq \
    -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log

  pigz *fq

  for sam in *.sam
  do
    samtools view -h -F 4 $sam | samtools sort -@ $threads -o $(basename $sam sam)bam
    rm $sam
  done
done
```

## `merge_bwt_log` full help

### Function

Merge multiple Bowtie mapping log files into a single summary table.

### Parameters

| Parameter | Meaning |
|---|---|
| `-l`, `--list` | Bowtie log files |
| `-o`, `--output` | output prefix |
| `-n`, `--name` | comma-separated database labels, such as `rRNA,tRNA,ncRNA,mRNA,Genome` |

### Example

```bash
merge_bwt_log \
  -n rRNA,tRNA,ncRNA,mRNA,Genome \
  -l *log \
  -o RNA_seq \
  &>> merge_bowtie.log
```

### Output

| File | Description |
|---|---|
| `RNA_seq.txt` | merged mapping statistics |
| `merge_bowtie.log` | running log |

## STAR alignment

```bash
genome='../../1.reference/star-index/'
threads=12

for fastq in ../2.bowtie/*.noncrna.fq.gz
do
  output=$(basename $fastq .noncrna.fq.gz)

  STAR --runThreadN $threads \
    --readFilesCommand zcat \
    --genomeDir $genome \
    --readFilesIn $fastq \
    --outFileNamePrefix $output \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --quantMode TranscriptomeSAM GeneCounts \
    --outReadsUnmapped Fastx \
    --outSAMattributes All \
    --alignEndsType Local \
    --outFilterMultimapNmax 3 \
    --outFilterMismatchNmax 1 \
    --alignIntronMax 10000 \
    --outFilterMatchNmin 20

  pigz *mate1

  samtools sort -@ $threads \
    $output"Aligned.out.bam" \
    -o $output"Aligned.sortedByCoord.out.bam"

  samtools index -@ $threads $output"Aligned.sortedByCoord.out.bam"
  rm $output"Aligned.out.bam"
done
```

## RSEM quantification

```bash
for bam in ../3.star/*Aligned.toTranscriptome.out.bam
do
  rsem-calculate-expression -p 10 \
    --no-bam-output \
    --alignments \
    -q $bam \
    ../../1.reference/rsem-index/rsem \
    $(basename $bam Aligned.toTranscriptome.out.bam)
done
```

## `merge_rsem` full help

### Function

Merge selected columns from multiple RSEM result files.

### Parameters

| Parameter | Meaning |
|---|---|
| `-l`, `--list` | RSEM result files, such as `*.genes.results` or `*.isoforms.results` |
| `-o`, `--output` | output file |
| `-c`, `--column` | one of `expected_count`, `TPM`, or `FPKM` |

### Examples

```bash
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt

merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt
```

## Output interpretation

| Output | Meaning |
|---|---|
| classified BAM files | reads assigned to rRNA/tRNA/ncRNA/mRNA/genome |
| `Aligned.sortedByCoord.out.bam` | genome-aligned BAM |
| `Aligned.toTranscriptome.out.bam` | transcriptome-aligned BAM for RSEM/RiboParser |
| `*.genes.results` | RSEM gene-level quantification |
| `*.isoforms.results` | RSEM isoform-level quantification |
