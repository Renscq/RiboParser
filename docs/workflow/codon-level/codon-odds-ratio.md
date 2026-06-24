# 4.7.7 Codon odds ratio

## Function

`rpf_Odd_Ratio` compares codon-associated RPF enrichment between control and treatment sample groups. It builds two-dimensional count tables for codon-level signal, calculates odds-ratio-style statistics, corrects p-values, and summarizes significant codon enrichment or depletion.

## Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| Gene list | Optional transcript or gene list. If not provided, all genes are used. |
| Control and treatment sample names | Comma-separated sample names. These names must match columns in the merged RPF density table. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `-r` | Yes | Input RPF density table in TXT format. |
| `-l` | No | Gene name or transcript ID list in TXT format. If omitted, all genes are used. |
| `-o` | Yes | Output prefix. |
| `--thread` | No | Number of threads used for p-value calculation. Default: `1`. |
| `-s` | No | Ribosomal site used for odds-ratio calculation. Choices: `E`, `P`, `A`. Default: `P`. |
| `-f` | No | Reading frame used for odds-ratio calculation. Choices: `0`, `1`, `2`, `all`. Default: `all`. |
| `-c` | Yes | Comma-separated control sample names. |
| `-t` | Yes | Comma-separated treatment sample names. |
| `-n` | No | Normalize RPF counts to RPM before odds-ratio calculation. Disabled by default. |
| `-z` | No | Fill empty codon values with `0.01 * minimal value`. Disabled by default. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `50`. |
| `--fdr` | No | Significance column used for filtering. Choices: `bhfdr`, `pvalue`. Default: `bhfdr`. |
| `-v` | No | Significance threshold for the selected p-value or FDR column. Default: `0.05`. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `0`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `0`. |
| `--stop` | No | Remove stop codons from the codon-level summary table. Disabled by default. |
| `--scale` | No | Scaling method for relative odds-ratio summary values. Choices: `zscore`, `minmax`. Default: `minmax`. |

## Output files

| Output | Description |
|---|---|
| `<prefix>_codon_odd_ratio.txt` | Significant codon-level odds-ratio result table after filtering by the selected p-value or FDR threshold. |
| `<prefix>_sum_codon_odd_ratio.txt` | Codon-level summary table containing control/treatment counts, proportions, and relative values. |
| `<prefix>_odd_lineplot.pdf` / `.png` | Line plot summarizing significant differential codon counts and proportions. |
| `<prefix>_odd_scatter.pdf` / `.png` | Scatter plot comparing control and treatment codon counts/proportions. |

## Example

```bash
rpf_Odd_Ratio \
    -l ../../../1.reference/norm/gene.norm.txt \
    -r ../05.merge/RIBO_merged.txt \
    -c WT_ribo_YPD1,WT_ribo_YPD2,WT_ribo_YPD3 \
    -t ncs2d_ribo_YPD1,ncs2d_ribo_YPD2,ncs2d_ribo_YPD3 \
    --stop \
    -m 50 \
    -f 0 \
    -s A \
    --tis 10 \
    --tts 5 \
    --fdr bhfdr \
    -v 0.05 \
    -o A_site \
    &> A_site_odd_ratio.log
```

## Notes

- The sample names provided to `-c` and `-t` must exactly match sample columns in the merged density table.
- Use `--fdr bhfdr` for Benjamini-Hochberg FDR filtering or `--fdr pvalue` for raw p-value filtering.
- The command generates line and scatter plots in the current implementation. The internal barplot function exists but is not called by the command-line workflow.
- Interpret odds-ratio results together with pausing, occupancy, and read-depth information.

## Merge related results

Use `merge_odd_ratio` to merge multiple codon odds-ratio summary tables:

```bash
merge_odd_ratio -l *_sum_codon_odd_ratio.txt -o RIBO
```
