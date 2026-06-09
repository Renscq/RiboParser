# 4.5.6 Periodicity

Tri-nucleotide periodicity is a critical quality metric for Ribo-seq data. High-quality Ribo-seq should show strong frame preference.

## Interpretation

- High periodicity supports codon-resolution analysis.
- Low periodicity introduces frame ambiguity.
- Frame ambiguity can misassign A/P/E-site codons and generate false-positive pause sites.

## Command help

```bash
rpf_Periodicity -h
```

```text
usage: rpf_Periodicity [-h] -r RPF -o OUTPUT [-t TRANSCRIPT] [-m MIN]
                       [--tis TIS] [--tts TTS]

Required arguments:
  -r RPF         input RPF density file
  -o OUTPUT      output prefix
  -t TRANSCRIPT  transcript annotation in TXT format

Options:
  -m MIN         retain transcript with more than minimum RPFs
  --tis TIS      discard codons after TIS
  --tts TTS      discard codons before TTS
```

## RNA-seq periodicity check

```bash
cd ./3.rna-seq/5.riboparser/06.periodicity/

rpf_Periodicity \
  -r ../05.merge/RNA_merged.txt \
  -m 30 \
  --tis 0 \
  --tts 0 \
  -o RNA \
  &>> RNA.log
```

## Ribo-seq periodicity check

```bash
cd ./4.ribo-seq/5.riboparser/06.periodicity/

rpf_Periodicity \
  -r ../05.merge/RIBO_merged.txt \
  -m 30 \
  --tis 0 \
  --tts 0 \
  -o RIBO \
  &>> RIBO.log
```

## Output files

```text
RIBO_count_periodicity_plot.pdf
RIBO_count_periodicity_plot.png
RIBO_ratio_periodicity_plot.pdf
RIBO_ratio_periodicity_plot.png
RIBO_periodicity.txt
RIBO.log
```
