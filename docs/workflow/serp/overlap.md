# 4.9.1 SeRP overlap

## Purpose

`serp_overlap` evaluates overlap among SeRP-related regions, peaks, or signal blocks.

## Input files

| Input | Description |
|---|---|
| region files | BED-like or SeRP peak files |
| annotation | optional gene or transcript annotation |
| output prefix | analysis prefix |

## Parameters

| Parameter | Meaning |
|---|---|
| `-i / --input` | input SeRP region or peak file |
| `-a / --annotation` | optional annotation file |
| `-o / --output` | output prefix |
| `--min-overlap` | minimum overlap length or fraction |
| `--strand` | strand-aware overlap mode |
| `--merge` | merge overlapping regions before comparison |

## Example

```bash
serp_overlap \
  -i serp_regions.bed \
  -a gene.norm.txt \
  -o serp_overlap
```

## Output files

| Output | Description |
|---|---|
| `serp_overlap.txt` | overlap summary table |
| `serp_overlap.bed` | overlapping regions if generated |
| `serp_overlap.log` | running log |

## Result interpretation

Use overlap analysis to compare SeRP peaks across samples or to intersect peaks with annotated features.
