# 4.6.4 Gene density retrieval

## Purpose

`rpf_Retrieve` extracts formatted gene-level density for downstream visualization, custom plotting, and external analysis.

## Input files

| Input | Description |
|---|---|
| density file | `RNA_merged.txt` or `RIBO_merged.txt` |
| gene/transcript annotation | `gene.norm.txt` or selected gene list |
| minimum read count | optional filtering threshold |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-o / --output` | output prefix |
| `-l / --list` | gene list or transcript annotation |
| `-m / --min` | minimum read count |
| `-n` | normalize read count to RPM |
| `-f` | melt three-column data of each sample to one column |
| `-s` | split gene density into individual TXT files |

## Example

```bash
rpf_Retrieve \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  -f \
  -n \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_retrieve.txt` | formatted retrieved density table |
| `RIBO.log` | running log |
| `split gene files` | optional per-gene density files when `-s` is used |

## Result interpretation

Use retrieved density for custom visualization, IGV-like plotting, gene-level inspection, or external statistical modeling.
