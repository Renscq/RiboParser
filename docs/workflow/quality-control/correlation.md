# 4.5.9 Correlation

## Purpose

`rpf_Corr` calculates sample correlation based on gene-level and/or RPF-level density.

## Input files

| Input | Description |
|---|---|
| density file | `RNA_merged.txt` or `RIBO_merged.txt` |
| output prefix | analysis prefix |

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
| `RIBO_gene_corr_f0.txt` | gene-level frame 0 correlation |
| `RIBO_gene_corr_f1.txt` | gene-level frame 1 correlation |
| `RIBO_gene_corr_f2.txt` | gene-level frame 2 correlation |
| `RIBO_gene_corr_frame.txt` | gene-level frame-combined correlation |
| `RIBO_gene_correlation_plot.pdf/png` | gene-level correlation heatmap |
| `RIBO_rpf_corr_f0.txt` | RPF-level frame 0 correlation |
| `RIBO_rpf_corr_f1.txt` | RPF-level frame 1 correlation |
| `RIBO_rpf_corr_f2.txt` | RPF-level frame 2 correlation |
| `RIBO_rpf_corr_frame.txt` | RPF-level frame-combined correlation |
| `RIBO_rpf_correlation_plot.pdf/png` | RPF-level correlation heatmap |

## Result interpretation

Biological replicates should show high correlation. Low correlation may reflect batch effects, poor alignment, poor periodicity, or biological heterogeneity.
