# 4.7.1 Codon pausing score

## Function

`rpf_Pausing` calculates relative codon pausing scores from a merged RPF density table. The command compares codon-level RPF density against either the gene-level average density or the local upstream/downstream background, then summarizes pausing scores by gene and codon.

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
| `-s` | No | Ribosomal site used for pausing calculation. Choices: `E`, `P`, `A`. Default: `P`. |
| `-f` | No | Reading frame used for pausing calculation. Choices: `0`, `1`, `2`, `all`. Default: `all`. |
| `-b` | No | Number of upstream and downstream codons used as local background. Default: `2`. Set `0` to use the gene-level mean CDS density as background. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `50`. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `10`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `5`. |
| `-n` | No | Normalize RPF counts to RPM before calculating pausing scores. Disabled by default. |
| `--scale` | No | Scaling method for codon-level summary values. Choices: `zscore`, `minmax`. Default: `minmax`. |
| `--stop` | No | Remove stop codons from codon-level summary tables and heatmaps. Disabled by default. |
| `--ind` | No | Summarize valid pausing scores using codons with non-zero density in each individual sample. Disabled by default. |
| `--fig` | No | Draw per-gene pausing plots. Choices: `none`, `png`, `pdf`. Default: `none`. This can take a long time for large datasets. |
| `--all` | No | Output all codon-level pausing scores for each gene. Disabled by default. |

## Output files

| Output | Description |
|---|---|
| `<prefix>_cds_pausing_score.txt` | Gene-level CDS pausing score table. |
| `<prefix>_cds_codon_pausing_score.txt` | Per-gene and per-codon pausing score summary. |
| `<prefix>_sum_codon_pausing_score.txt` | Codon-level summary table containing codon counts, RPF counts, absolute pausing scores, and relative pausing scores. |
| `<prefix>_all_pausing_score.txt` | Full codon-level pausing table. Generated only with `--all`. |
| `<prefix>_total_pausing_heatplot.pdf` / `.png` | Heatmap of total codon pausing scores. |
| `<prefix>_valid_pausing_heatplot.pdf` / `.png` | Heatmap of valid codon pausing scores. |
| `<gene>_pausing_plot.pdf` / `.png` | Per-gene pausing plot. Generated only when `--fig pdf` or `--fig png` is used. |

## Example

```bash
for site in E P A
do
    rpf_Pausing \
        -l ../../../1.reference/norm/gene.norm.txt \
        -r ../05.merge/RIBO_merged.txt \
        -b 0 \
        --stop \
        -m 30 \
        -s ${site} \
        -f 0 \
        --scale minmax \
        -o ${site}_site \
        &> ${site}_site.log
done
```

## Notes

- High pausing scores indicate local RPF enrichment relative to the selected background.
- `-b 0` uses the average CDS density of each gene as the background. Positive `-b` values use a local window around each codon.
- Interpret candidate pausing codons together with RPF coverage, digestion bias, periodicity, and replicate consistency.

## Merge related results

Use `merge_pausing` to merge multiple pausing result tables:

```bash
merge_pausing -l *_sum_codon_pausing_score.txt -o RIBO
```
