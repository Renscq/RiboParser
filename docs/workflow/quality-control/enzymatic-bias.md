# 4.5.2 Enzymatic bias

Ribonuclease digestion and ligation steps can introduce sequence bias in Ribo-seq libraries.

## Biological rationale

`Ribonuclease digestion bias`: ribonucleases may preferentially cleave specific sequence or structural contexts, generating uneven fragment representation.

`Ligation bias`: ligation efficiency can vary by sequence context, fragment length, or RNA secondary structure.

These biases should be evaluated because they can affect codon-level and positional RPF density interpretation.

## Command help

```bash
rpf_Digest -h
```

```text
usage: rpf_Digest [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -o OUTPUT
                  [-l] [--scale] [-m MIN] [-M MAX]

Required arguments:
  -t TRANSCRIPT    transcript annotation in TXT format
  -s SEQUENCE      transcript sequence in FASTA format
  -b BAM           mapping file in BAM format
  -o OUTPUT        output prefix

Options:
  -l               only retain transcript with longest CDS of each gene
  --scale          scale motif matrix
  -m MIN           minimum read length to keep
  -M MAX           maximum read length to keep
```

## RNA-seq example

```bash
cd ./3.rna-seq/5.riboparser/02.digestion/

for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rpf_Digest \
    -b $bam \
    -m 25 \
    -M 50 \
    --scale \
    -s ../../../1.reference/norm/gene.norm.rna.fa \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done

merge_digestion -l *pwm.txt -o RNA
```

## Ribo-seq example

```bash
cd ./4.ribo-seq/5.riboparser/02.digestion/

for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rpf_Digest \
    -b $bam \
    -m 27 \
    -M 33 \
    --scale \
    -s ../../../1.reference/norm/gene.norm.rna.fa \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done

merge_digestion -l *pwm.txt -o RIBO
```

## Output files

```text
SRR1944912_3end_counts.txt
SRR1944912_3end_pwm.txt
SRR1944912_3end_seqlogo2.pdf
SRR1944912_5end_counts.txt
SRR1944912_5end_pwm.txt
SRR1944912_5end_seqlogo2.pdf
SRR1944912_digestion_sites.txt
SRR1944912_scaled_digestion_sites_plot.pdf
SRR1944912.log
```
