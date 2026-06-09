# 4.7.3 Codon decoding time

## Purpose

`rpf_CDT` estimates codon decoding time by combining Ribo-seq and RNA-seq density.

## Input files

| Input | Description |
|---|---|
| Ribo-seq density | `RIBO_merged.txt` |
| RNA-seq density | `RNA_merged.txt` |
| gene list/annotation | `gene.norm.txt` |
| site/frame settings | E/P/A site and reading frame |

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
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `--scale` | scale method |
| `--stop` | remove stop codon |

## Example

```bash
for sites in E P A
do
  rpf_CDT \
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
| `A_site_cdt.txt` | codon decoding time table |
| `A_site_cdt_corr.txt` | CDT correlation table |
| `A_site_cdt_corrplot.pdf/png` | CDT correlation plot |
| `A_site_cdt_heatplot.pdf/png` | CDT heatmap |

## Result interpretation

CDT attempts to normalize translational density by transcript abundance. It is useful for comparing codon-level elongation patterns across conditions.

## Merge related results

### `merge_cdt`

```bash
merge_cdt -l *cdt.txt -o RIBO
```
