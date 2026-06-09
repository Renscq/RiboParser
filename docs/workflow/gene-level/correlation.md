# 4.6.3 Gene correlation

## Purpose

`rpf_Corr` evaluates sample reproducibility at gene and RPF density levels.

## Input files

| Input | Description |
|---|---|
| density file | `RNA_merged.txt` or `RIBO_merged.txt` |
| sample names | stored in merged density file |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-o / --output` | output prefix |

## Example

```bash
rpf_Corr \
  -r ../05.merge/RIBO_merged.txt \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_gene_correlation_plot.pdf/png` | gene-level correlation heatmap |
| `RIBO_rpf_correlation_plot.pdf/png` | RPF-level correlation heatmap |
| `RIBO_gene_corr_*.txt` | gene-level correlation tables |
| `RIBO_rpf_corr_*.txt` | RPF-level correlation tables |

## Result interpretation

Use this page to determine whether replicates cluster together and whether outlier samples should be removed.
