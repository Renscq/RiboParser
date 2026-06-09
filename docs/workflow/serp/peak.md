# 4.9.2 SeRP peak

## Purpose

`serp_peak` identifies enriched SeRP signal peaks from signal profiles.

## Input files

| Input | Description |
|---|---|
| signal file | bedGraph, WIG, or density-like signal file |
| annotation | optional gene/transcript annotation |
| thresholds | peak calling thresholds |

## Parameters

| Parameter | Meaning |
|---|---|
| `-i / --input` | input signal file |
| `-o / --output` | output prefix |
| `-a / --annotation` | annotation file |
| `--min-signal` | minimum peak signal |
| `--min-width` | minimum peak width |
| `--max-width` | maximum peak width |
| `--merge-distance` | distance for merging adjacent peaks |
| `--strand` | strand-aware peak calling |

## Example

```bash
serp_peak \
  -i serp_signal.bedgraph \
  -a gene.norm.txt \
  --min-signal 5 \
  --min-width 3 \
  -o serp_peak
```

## Output files

| Output | Description |
|---|---|
| `serp_peak.txt` | peak table |
| `serp_peak.bed` | peak regions |
| `serp_peak.log` | running log |

## Result interpretation

Peak calling thresholds should be chosen based on library depth, background signal, and biological replicate consistency.
