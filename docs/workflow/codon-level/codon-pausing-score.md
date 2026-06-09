# 4.7.1 Codon pausing score

## Purpose

`rpf_Pausing` calculates relative codon pausing score by comparing codon-level density against gene-level background.

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
| `-l / --list` | gene list or annotation; default can be whole transcript set |
| `-o / --output` | output prefix |
| `-s` | ribosomal site: `E`, `P`, or `A` |
| `-f` | reading frame: `0`, `1`, `2`, or `all` |
| `-b` | background codon number |
| `-m` | minimum RPF count |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `-n` | normalize to RPM |
| `--scale` | scale method: `zscore` or `minmax` |
| `--stop` | remove stop codon |
| `--fig` | figure output type |
| `--all` | output all gene-level pausing scores |

## Example

```bash
for sites in E P A
do
  rpf_Pausing \
    -l ../../../1.reference/norm/gene.norm.txt \
    -r ../05.merge/RIBO_merged.txt \
    -b 0 \
    --stop \
    -m 30 \
    -s $sites \
    -f 0 \
    --scale minmax \
    -o "$sites"_site \
    &>> "$sites"_site.log
done
```

## Output files

| Output | Description |
|---|---|
| `A_site_cds_codon_pausing_score.txt` | CDS codon-level pausing scores |
| `A_site_cds_pausing_score.txt` | gene-level pausing table |
| `A_site_sum_codon_pausing_score.txt` | summed codon pausing scores |
| `A_site_total_pausing_heatplot.pdf/png` | all codon heatmap |
| `A_site_valid_pausing_heatplot.pdf/png` | filtered codon heatmap |

## Result interpretation

High pausing scores indicate locally enriched RPF density. Interpret candidate pauses together with coverage, digestion bias, and replicate consistency.

## Merge related results

### `merge_pausing`

```bash
merge_pausing -l *pausing_score.txt -o RIBO
```
