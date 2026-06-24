# 4.7.3 Codon decoding time

## Function

`rpf_CDT` estimates codon decoding time by combining matched Ribo-seq and RNA-seq density tables. It normalizes codon-level RPF signal by RNA-seq signal and reports codon-level decoding-time summaries and plots.

## Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| RNA density table | Merged RNA-seq density table, usually `RNA_merged.txt`. |
| Gene annotation table | Normalized gene annotation table, for example `gene.norm.txt`. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `--rpf` | Yes | Input RPF density table in TXT format. |
| `--rna` | Yes | Input RNA-seq density table in TXT format. |
| `-l` | Yes | Gene annotation or gene list file in TXT format. |
| `-o` | Yes | Output prefix. |
| `-s` | No | Ribosomal site used for CDT calculation. Choices: `E`, `P`, `A`. Default: `P`. |
| `-f` | No | Reading frame used for CDT calculation. Choices: `0`, `1`, `2`, `all`. Default: `all`. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `30`. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `15`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `5`. |
| `--scale` | No | Scaling method for relative CDT values. Choices: `zscore`, `minmax`. Default: `minmax`. |
| `--stop` | No | Remove stop codons from codon-level summary tables and heatmaps. Disabled by default. |

## Output files

| Output | Description |
|---|---|
| `<prefix>_cdt.txt` | Codon decoding time summary table. |
| `<prefix>_cdt_corr.txt` | Correlation matrix calculated from normalized CDT values. |
| `<prefix>_cdt_corrplot.pdf` / `.png` | CDT correlation heatmap. |
| `<prefix>_cdt_heatplot.pdf` / `.png` | CDT heatmap. |

## Example

```bash
for site in E P A
do
    rpf_CDT \
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

- CDT requires matched Ribo-seq and RNA-seq density tables. Gene names must be consistent between the two files.
- The command uses the gene annotation file to define transcript structure and transcript length information.
- CDT is useful for comparing codon-level elongation-related patterns across samples or conditions.

## Merge related results

Use `merge_cdt` to merge multiple CDT result tables:

```bash
merge_cdt -l *_cdt.txt -o RIBO
```
