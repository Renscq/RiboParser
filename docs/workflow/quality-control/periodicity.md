# 4.5.6 Periodicity

## Purpose

`rpf_Periodicity` calculates frame-specific read distribution and 3-nt periodicity. It is a key Ribo-seq QC metric.

## Input files

| Input | Description |
|---|---|
| density file | `RIBO_merged.txt` or `RNA_merged.txt` |
| transcript annotation | optional `gene.norm.txt` |
| minimum read count | threshold for included transcripts |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-o / --output` | output prefix |
| `-t / --transcript` | normalized transcript annotation |
| `-m / --min` | minimum read count |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |

## Example

```bash
rpf_Periodicity \
  -r ../05.merge/RIBO_merged.txt \
  -m 30 \
  --tis 0 \
  --tts 0 \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_periodicity.txt` | frame-specific periodicity table |
| `RIBO_count_periodicity_plot.pdf/png` | periodicity count plot |
| `RIBO_ratio_periodicity_plot.pdf/png` | periodicity ratio plot |
| `RIBO.log` | running log |

## Result interpretation

High-quality Ribo-seq data should show strong in-frame enrichment. Low periodicity indicates that codon-level analyses should be interpreted cautiously.

## Merge related results

### `merge_period`

| Parameter | Meaning |
|---|---|
| `-l` | input `*periodicity.txt` files |
| `-o` | output prefix |

```bash
merge_period -l *periodicity.txt -o RIBO
```
