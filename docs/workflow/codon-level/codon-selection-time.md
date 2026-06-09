# 4.7.4 Codon selection time

## Purpose

`rpf_CST` estimates codon selection time from RNA-seq and Ribo-seq density profiles.

## Input files

| Input | Description |
|---|---|
| Ribo-seq density | `RIBO_merged.txt` |
| RNA-seq density | `RNA_merged.txt` |
| gene list/annotation | `gene.norm.txt` |
| iteration setting | number of CST iterations if needed |

## Parameters

| Parameter | Meaning |
|---|---|
| `--rpf` | input Ribo-seq density file |
| `--rna` | input RNA-seq density file |
| `-l / --list` | gene list or annotation |
| `-o / --output` | output prefix |
| `-s` | ribosomal site |
| `-f` | reading frame |
| `-m` | minimum RPF count |
| `-t` | iteration number |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `--scale` | scale method |
| `--stop` | remove stop codon |

## Example

```bash
for sites in E P A
do
  rpf_CST \
    -l ../../../1.reference/norm/gene.norm.txt \
    --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
    --rpf ../05.merge/RIBO_merged.txt \
    --stop \
    -m 50 \
    -f 0 \
    -s $sites \
    --tis 10 \
    --tts 5 \
    -o "$sites"_site \
    &>> "$sites"_site.log
done
```

## Output files

| Output | Description |
|---|---|
| `A_site_codon_selection_time.txt` | codon selection time table |
| `A_site_iterative_codon_selection_time.txt` | iterative CST table |
| `A_site_cst_corr.txt` | CST correlation table |
| `A_site_cst_corrplot.pdf/png` | CST correlation plot |
| `A_site_cst_heatplot.pdf/png` | CST heatmap |

## Result interpretation

CST is useful for evaluating codon-level selection or elongation-related differences across experimental groups.

## Merge related results

### `merge_cst`

```bash
merge_cst -l *cst.txt -o RIBO
```
