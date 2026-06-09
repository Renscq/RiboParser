# 4.7.2 Codon occupancy

## Purpose

`rpf_Occupancy` calculates codon occupancy at E/P/A sites and compares codon-specific RPF density across samples.

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
| `-r / --rpf` | input Ribo-seq density file |
| `-l / --list` | gene list or annotation |
| `-o / --output` | output prefix |
| `-s` | ribosomal site: `E`, `P`, or `A` |
| `-f` | reading frame |
| `-m` | minimum RPF count |
| `-n` | normalize to RPM |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `--scale` | scale method |
| `--stop` | remove stop codon |
| `--all` | output all RPF density |

## Example

```bash
for sites in E P A
do
  rpf_Occupancy \
    -l ../../../1.reference/norm/gene.norm.txt \
    -r ../05.merge/RIBO_merged.txt \
    -m 30 \
    -s "$sites" \
    -f 0 \
    --stop \
    --scale minmax \
    -o "$sites"_site \
    &>> "$sites"_site.log
done
```

## Output files

| Output | Description |
|---|---|
| `A_site_codon_density.txt` | codon-level density table |
| `A_site_codon_occupancy.txt` | codon occupancy table |
| `A_site_occupancy_corr.txt` | occupancy correlation table |
| `A_site_occupancy_corrplot.pdf/png` | correlation plot |
| `A_site_occupancy_heatplot.pdf/png` | absolute occupancy heatmap |
| `A_site_occupancy_relative_heatplot.pdf/png` | relative occupancy heatmap |
| `A_site_occupancy_relative_lineplot.pdf/png` | relative occupancy line plot |

## Result interpretation

Codon occupancy reflects codon-specific ribosome density. It should be interpreted after confirming frame periodicity and sufficient read depth.

## Merge related results

### `merge_occupancy`

```bash
merge_occupancy -l *occupancy.txt -o RIBO
```
