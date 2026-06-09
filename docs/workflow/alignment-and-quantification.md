# 4.4 Alignment and quantification

This step classifies reads against different reference databases and then aligns remaining reads to the genome with STAR.

## RNA-seq read classification with Bowtie

This step is usually more important for rRNA-depletion RNA-seq libraries. For oligo(dT)-based RNA-seq, most reads are expected to be mRNA.

```bash
mkdir -p ./sce/3.rna-seq/2.bowtie/
cd ./sce/3.rna-seq/2.bowtie/

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

## Ribo-seq read classification with Bowtie

The Ribo-seq workflow uses the same classification strategy:

```bash
mkdir -p ./sce/4.ribo-seq/2.bowtie/
cd ./sce/4.ribo-seq/2.bowtie/

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

## Merge Bowtie mapping logs

```bash
merge_bwt_log -h
```

```text
usage: merge_bwt_log [-h] -l LIST [LIST ...] -o OUTPUT [-n NAME]

Required arguments:
  -l LIST    list of Bowtie mapping log files
  -o OUTPUT  output prefix

Options:
  -n NAME    names of each database; default: rRNA,tRNA,ncRNA,mRNA,Genome
```

```bash
merge_bwt_log \
  -n rRNA,tRNA,ncRNA,mRNA,Genome \
  -l *log \
  -o RNA_seq \
  &>> merge_bowtie.log
```

## RNA-seq STAR alignment

```bash
mkdir -p ./sce/3.rna-seq/3.star/
cd ./sce/3.rna-seq/3.star/

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

## Ribo-seq STAR alignment

```bash
mkdir -p ./sce/4.ribo-seq/3.star/
cd ./sce/4.ribo-seq/3.star/

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

## RNA-seq expression quantification with RSEM

```bash
mkdir -p ./sce/3.rna-seq/4.quantification/
cd ./sce/3.rna-seq/4.quantification/

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

## Ribo-seq expression quantification with RSEM

```bash
mkdir -p ./sce/4.ribo-seq/4.quantification/
cd ./sce/4.ribo-seq/4.quantification/

for bam in ../3.star/*Aligned.toTranscriptome.out.bam
do
  rsem-calculate-expression -p 12 \
    --no-bam-output \
    --alignments \
    -q $bam \
    ../../1.reference/rsem-index/rsem \
    $(basename $bam Aligned.toTranscriptome.out.bam)
done
```

## Merge RSEM results

```bash
merge_rsem -h
```

```text
usage: merge_rsem [-h] -l LIST [LIST ...] -o OUTPUT [-c {expected_count,TPM,FPKM}]

Required arguments:
  -l LIST     list of RSEM result files
  -o OUTPUT   output file name

Options:
  -c COLUMN   expected_count, TPM, or FPKM
```

```bash
# merge gene expression
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt &>> merge_rsem.log
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt &>> merge_rsem.log
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt &>> merge_rsem.log

# merge isoform expression
merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt &>> merge_rsem.log
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt &>> merge_rsem.log
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt &>> merge_rsem.log
```
