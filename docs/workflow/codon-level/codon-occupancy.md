# 4.7.2 Codon occupancy

## Function

`rpf_Occupancy` calculates codon occupancy from a merged RPF density table. The command normalizes codon-level density by the average density of each gene, summarizes codon occupancy across genes, and generates correlation and heatmap plots.

## Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| Gene list | Optional transcript or gene list. If not provided, all genes in the density table are used. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `-r` | Yes | Input RPF density table in TXT format. |
| `-l` | No | Gene name or transcript ID list in TXT format. If omitted, all genes are used. |
| `-o` | Yes | Output prefix. |
| `-s` | No | Ribosomal site used for occupancy calculation. Choices: `E`, `P`, `A`. Default: `P`. |
| `-f` | No | Reading frame used for occupancy calculation. Choices: `0`, `1`, `2`, `all`. Default: `all`. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `30`. |
| `-n` | No | Normalize RPF counts to RPM before calculating occupancy. Disabled by default. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `15`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `5`. |
| `--scale` | No | Scaling method for relative occupancy. Choices: `zscore`, `minmax`. Default: `minmax`. |
| `--stop` | No | Remove stop codons from codon-level summary tables and heatmaps. Disabled by default. |
| `--all` | No | Output all filtered RPF density and per-gene codon density tables. Disabled by default. |

## Output files

| Output | Description |
|---|---|
| `<prefix>_codon_occupancy.txt` | Codon-level occupancy summary table containing total density, absolute occupancy, and relative occupancy. |
| `<prefix>_occupancy_corr.txt` | Correlation matrix calculated from absolute codon occupancy. |
| `<prefix>_occupancy_corrplot.pdf` / `.png` | Occupancy correlation heatmap. |
| `<prefix>_occupancy_heatplot.pdf` / `.png` | Absolute codon occupancy heatmap. |
| `<prefix>_occupancy_relative_heatplot.pdf` / `.png` | Relative codon occupancy heatmap. |
| `<prefix>_occupancy_relative_lineplot.pdf` / `.png` | Relative codon occupancy line plot. |
| `<prefix>_rpf_density.txt` | Filtered RPF density table. Generated only with `--all`. |
| `<prefix>_codon_density.txt` | Per-gene and per-codon density table. Generated only with `--all`. |

## Example

```bash
for site in E P A
do
    rpf_Occupancy \
        -l ../../../1.reference/norm/gene.norm.txt \
        -r ../05.merge/RIBO_merged.txt \
        -m 30 \
        -s ${site} \
        -f 0 \
        --stop \
        --scale minmax \
        -o ${site}_site \
        &> ${site}_site.log
done
```

## Notes

- Occupancy values reflect codon-specific ribosome density after gene-level normalization.
- Use `--stop` when stop codons should be excluded from summary plots.
- Use `--all` when you need intermediate density tables for troubleshooting or downstream analysis.

## Merge related results

Use `merge_occupancy` to merge multiple codon occupancy tables:

```bash
merge_occupancy -l *_codon_occupancy.txt -o RIBO
```
