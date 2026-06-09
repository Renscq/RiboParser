# 4.5.2 Enzymatic bias

## Purpose

`rpf_Digest` evaluates 5' and 3' end sequence preference caused by nuclease digestion, ligation, or library construction bias.

## Input files

| Input | Description |
|---|---|
| BAM | filtered BAM from `rpf_Check` |
| transcript annotation | `gene.norm.txt` |
| transcript sequence | `gene.norm.rna.fa` |
| read length range | Ribo-seq usually 27–33 nt; RNA-seq can use a broader range |

## Parameters

| Parameter | Meaning |
|---|---|
| `-b / --bam` | input BAM file |
| `-t / --transcript` | normalized transcript annotation |
| `-s / --sequence` | transcript FASTA sequence |
| `-o / --output` | output prefix |
| `-m / --min` | minimum read length |
| `-M / --max` | maximum read length |
| `-l` | use the longest transcript per gene |
| `--scale` | scale motif matrix and generate scaled digestion plot |

## Example

```bash
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
```

## Output files

| Output | Description |
|---|---|
| `sample_5end_counts.txt` | 5' end nucleotide counts |
| `sample_5end_pwm.txt` | 5' end PWM |
| `sample_5end_seqlogo2.pdf` | 5' end sequence logo |
| `sample_3end_counts.txt` | 3' end nucleotide counts |
| `sample_3end_pwm.txt` | 3' end PWM |
| `sample_3end_seqlogo2.pdf` | 3' end sequence logo |
| `sample_digestion_sites.txt` | digestion site distribution |
| `sample_scaled_digestion_sites_plot.pdf` | scaled digestion profile |
| `sample.log` | running log |

## Result interpretation

Strong end-sequence preference indicates digestion or ligation bias. This is especially important for codon-level interpretation because local sequence bias can mimic ribosome pausing.

## Merge related results

### `merge_digestion`

| Parameter | Meaning |
|---|---|
| `-l` | digestion PWM files, for example `*pwm.txt` |
| `-o` | output prefix |

```bash
merge_digestion -l *pwm.txt -o RIBO
```
