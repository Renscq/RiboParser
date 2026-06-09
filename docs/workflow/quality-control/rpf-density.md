# 4.5.4 RPF density

This step converts a BAM file into a transcript-level read density table.

## RNA-seq density

```bash
cd ./3.rna-seq/5.riboparser/04.density/

for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rna_Density \
    -b $bam \
    -m 27 \
    -M 33 \
    -l \
    --thread 10 \
    -p ../03.offset/$prefix_name"_offset.txt" \
    -s ../../../1.reference/norm/gene.norm.rna.fa \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

## Ribo-seq command help

```bash
rpf_Density -h
```

```text
usage: rpf_Density [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -p PSITE -o OUTPUT
                   [-l] [-m MIN] [-M MAX] [--period PERIODICITY]
                   [--silence] [--thread THREAD]

Required arguments:
  -t TRANSCRIPT        transcript annotation in TXT format
  -s SEQUENCE          transcript sequence in FASTA format
  -b BAM               mapping file in BAM format
  -p PSITE             P-site offset table
  -o OUTPUT            output prefix

Options:
  -l                   only retain transcript with longest CDS of each gene
  -m MIN               minimum read length to keep
  -M MAX               maximum read length to keep
  --period PERIODICITY minimum 3-nt periodicity to keep
  --thread THREAD      number of threads
```

## Ribo-seq density

```bash
cd ./4.ribo-seq/5.riboparser/04.density/

for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rpf_Density \
    -b $bam \
    -m 27 \
    -M 33 \
    --period 40 \
    -l \
    --thread 12 \
    -p ../03.offset/$prefix_name"_SSCBM_offset.txt" \
    -s ../../../1.reference/norm/gene.norm.rna.fa \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

## Output files

```text
SRR1944912.log
SRR1944912_rpf.txt
```
