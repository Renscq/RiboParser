# 4.7.7 Codon odds ratio

## Purpose

`rpf_Odd_Ratio` calculates enrichment or depletion of codon-associated RPF signal as an odds-ratio-style metric.

## Input files

| Input | Description |
|---|---|
| Ribo-seq density | `RIBO_merged.txt` |
| gene list/annotation | `gene.norm.txt` or selected gene list |
| site setting | E, P, or A site |
| frame setting | usually frame 0 for high-quality Ribo-seq |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-l / --list` | gene list or annotation |
| `-o / --output` | output prefix |
| `-s` | ribosomal site |
| `-f` | reading frame |
| `-m` | minimum read count |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `--stop` | remove stop codon |

## Example

```bash
for sites in E P A
do
  rpf_Odd_Ratio \
    -l ../../../1.reference/norm/gene.norm.txt \
    -r ../05.merge/RIBO_merged.txt \
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
| `A_site_odd_ratio.txt` | codon odds ratio table |
| `A_site_odd_ratio_plot.pdf/png` | odds ratio plot when generated |
| `A_site.log` | running log |

## Result interpretation

Odds ratio values help summarize codon enrichment relative to background. Interpret together with pausing and occupancy results.

## Merge related results

### `merge_odd_ratio`

```bash
merge_odd_ratio -l *odd_ratio.txt -o RIBO
```
