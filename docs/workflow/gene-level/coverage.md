# 4.6.2 Gene coverage

## Purpose

`rpf_Coverage` and `rpf_Percent` evaluate gene-body coverage, regional read distribution, and frame-aware coverage.

## Input files

| Input | Description |
|---|---|
| density file | merged RNA/Ribo density |
| transcript annotation | `gene.norm.txt` |
| bin settings | number of bins for 5'UTR, CDS, and 3'UTR |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-t / --transcript` | normalized transcript annotation |
| `-o / --output` | output prefix |
| `-f` | reading frame |
| `-m` | minimum read count |
| `-b` | bin settings such as `10,150,10` |
| `-n` | normalize to RPM |
| `--outlier` | remove outlier profiles |
| `--heat` | draw heatmap |
| `--bar` | draw barplot |

## Example

```bash
rpf_Coverage \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  --outlier \
  -b 10,150,10 \
  -n \
  --heat \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_*_coverage.txt` | binned coverage matrix |
| `RIBO_*_heat_plot.pdf/png` | coverage heatmap |
| `RIBO_*_coverage_bar_plot.pdf/png` | coverage barplot |
| `RIBO_*_coverage_line_plot.pdf/png` | average coverage profile |

## Result interpretation

Coverage profiles reveal whether reads are concentrated in CDS, 5' ends, 3' ends, or uniformly distributed.
