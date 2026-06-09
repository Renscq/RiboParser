# 4.6.1 Gene quantification

## Purpose

`rpf_Quant` quantifies Ribo-seq density in CDS, UTR, or selected transcript regions. It is typically used after density merging.

## Input files

| Input | Description |
|---|---|
| density file | `RIBO_merged.txt` |
| output prefix | analysis prefix |
| TIS/TTS trimming | codon trimming parameters to remove initiation/termination artifacts |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-o / --output` | output prefix |
| `-f` | reading frame: `0`, `1`, `2`, or `all` |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `--utr5` | quantify 5'UTR |
| `--utr3` | quantify 3'UTR |

## Example

```bash
rpf_Quant \
  -r ../05.merge/RIBO_merged.txt \
  --tis 15 \
  --tts 5 \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_cds_rpf_quant.txt` | raw CDS RPF count |
| `RIBO_cds_rpm_quant.txt` | CDS RPM table |
| `RIBO_cds_rpkm_quant.txt` | CDS RPKM table |
| `RIBO_cds_tpm_quant.txt` | CDS TPM table |
| `RIBO_cds_rpm_bar_plot.pdf/png` | RPM barplot |
| `RIBO_cds_rpm_cdf_plot.pdf/png` | RPM cumulative distribution |
| `RIBO_cds_rpm_heatmap.pdf/png` | RPM heatmap |
| `RIBO_cds_rpm_pca_plot.pdf/png` | PCA plot |
| `RIBO_cds_rpm_pca.txt` | PCA coordinates |
| `RIBO_total.txt` | total read summary |

## Result interpretation

Exclude codons near start and stop codons to avoid initiation/termination artifacts. For standard gene-level translation quantification, frame 0 and CDS-only density are usually preferred.

## Merge related results

### `merge_quant`

| Parameter | Meaning |
|---|---|
| `-l` | quantification result files |
| `-o` | output prefix |

```bash
merge_quant -l *quant.txt -o RIBO
```
