# 4.7.4 Codon selection time

## Function

`rpf_CST` estimates codon selection time using matched Ribo-seq and RNA-seq density tables. It iteratively calculates codon-level selection-time values and outputs both detailed iterative results and final codon-level summaries.

## Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| RNA density table | Merged RNA-seq density table, usually `RNA_merged.txt`. |
| Gene list | Optional transcript or gene list. If not provided, all genes shared by the two density tables are used. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `--rpf` | Yes | Input RPF density table in TXT format. |
| `--rna` | Yes | Input RNA-seq density table in TXT format. |
| `-l` | No | Gene name or transcript ID list in TXT format. If omitted, all genes are used. |
| `-o` | Yes | Output prefix. |
| `-s` | No | Ribosomal site used for CST calculation. Choices: `E`, `P`, `A`. Default: `P`. |
| `-f` | No | Reading frame used for CST calculation. Choices: `0`, `1`, `2`, `all`. Default: `all`. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `30`. |
| `-t` | No | Number of iterations used for CST calculation. Default: `10`. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `0`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `0`. |
| `--scale` | No | Scaling method for relative CST values. Choices: `zscore`, `minmax`. Default: `minmax`. |
| `--stop` | No | Remove stop codons from codon-level summary tables and heatmaps. Disabled by default. |

## Output files

| Output | Description |
|---|---|
| `<prefix>_iterative_codon_selection_time.txt` | Detailed CST values for each iteration and sample. |
| `<prefix>_codon_selection_time.txt` | Final codon-level CST table containing absolute and relative CST values. |
| `<prefix>_cst_corr.txt` | Correlation matrix calculated from absolute CST values. |
| `<prefix>_cst_corrplot.pdf` / `.png` | CST correlation heatmap. |
| `<prefix>_cst_heatplot.pdf` / `.png` | CST heatmap. |

## Example

```bash
for site in E P A
do
    rpf_CST \
        -l ../../../1.reference/norm/gene.norm.txt \
        --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
        --rpf ../05.merge/RIBO_merged.txt \
        --stop \
        -m 50 \
        -f 0 \
        -s ${site} \
        --tis 10 \
        --tts 5 \
        -o ${site}_site \
        &> ${site}_site.log
done
```

## Notes

- CST requires matched Ribo-seq and RNA-seq density tables. Gene names must be consistent between the two files.
- Increase `-t` when more iterative refinement is needed, but this also increases runtime.
- CST and CDT are related but not identical. CDT reports a direct decoding-time estimate, while CST performs iterative codon selection-time estimation.

## Merge related results

Use `merge_cst` to merge multiple CST result tables:

```bash
merge_cst -l *_codon_selection_time.txt -o RIBO
```
