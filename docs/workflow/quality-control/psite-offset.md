# 4.5.3 P-site offset

Codon-resolution Ribo-seq analysis requires accurate assignment of the ribosomal A-, P-, and E-site codons for each RPF.

The offset is the distance from the 5' end of the RPF to the first nucleotide of the P-site codon. RiboParser supports two common methods:

- `RSBM`: ribosome structure-based model
- `SSCBM`: start/stop codon-based model

## RNA-seq offset table

Offset prediction is unnecessary for RNA-seq. A constant offset can be assigned to all read lengths.

```bash
cd ./3.rna-seq/5.riboparser/03.offset/

for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rna_Offset \
    -m 27 \
    -M 50 \
    -e 12 \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

## Ribo-seq command help

```bash
rpf_Offset -h
```

```text
usage: rpf_Offset [-h] -t TRANSCRIPT -b BAM -o OUTPUT
                  [--mode {SSCBM,RSBM}] [-a {both,tis,tts}]
                  [-l] [-m MIN] [-M MAX] [-p EXP_PEAK] [-s SHIFT]
                  [--silence] [-d]

Required arguments:
  -t TRANSCRIPT    transcript annotation in TXT format
  -b BAM           mapping file in BAM format
  -o OUTPUT        output prefix

Options:
  --mode           SSCBM or RSBM
  -a               align reads to both, TIS, or TTS
  -l               only retain transcript with longest CDS of each gene
  -m MIN           minimum read length
  -M MAX           maximum read length
  -p EXP_PEAK      expected RPF length fitted to ribosome structure
  -s SHIFT         P-site shift for different RPF lengths
  -d               output offset details
```

## Ribo-seq offset prediction

```bash
cd ./4.ribo-seq/5.riboparser/03.offset/

for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rpf_Offset \
    -b $bam \
    -m 27 \
    -M 33 \
    -p 30 \
    -d \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

## Merge offset results

```bash
merge_offset_detail -l *end.txt -o RIBO
merge_offset -l *SSCBM_offset.txt -o RIBO_SSCBM
merge_offset -l *RSBM_offset.txt -o RIBO_RSBM
```

## Output files

```text
SRR1944912_RSBM_offset.pdf
SRR1944912_RSBM_offset.png
SRR1944912_RSBM_offset.txt
SRR1944912_SSCBM_offset.pdf
SRR1944912_SSCBM_offset.png
SRR1944912_SSCBM_offset_scale.pdf
SRR1944912_SSCBM_offset_scale.png
SRR1944912_SSCBM_offset.txt
SRR1944912_tis_3end.txt
SRR1944912_tis_5end.txt
SRR1944912_tts_3end.txt
SRR1944912_tts_5end.txt
SRR1944912.log
```
