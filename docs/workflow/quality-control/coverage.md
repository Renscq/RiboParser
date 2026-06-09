# 4.5.8 Coverage

## Purpose

Coverage analysis evaluates how read density is distributed across normalized gene bodies or transcript regions.

## `rpf_Coverage`

### Input files

| Input | Description |
|---|---|
| density file | `RNA_merged.txt` or `RIBO_merged.txt` |
| transcript annotation | `gene.norm.txt` |
| output prefix | analysis prefix |

### Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-t / --transcript` | normalized transcript annotation |
| `-o / --output` | output prefix |
| `-f` | reading frame: `0`, `1`, `2`, or `all` |
| `-m` | minimum read count |
| `-b` | bins for 5'UTR, CDS, and 3'UTR, for example `10,150,10` |
| `-n` | normalize to RPM |
| `--thread` | number of threads |
| `--outlier` | remove outlier genes/transcripts |
| `--set` | `intersect` or `union` strategy for sample filtering |
| `--heat` | output heatmap |
| `--bar` | output barplot |

### Example

```bash
rpf_Coverage \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  --outlier \
  -b 10,150,10 \
  -n \
  --heat \
  -o RIBO \
  &>> RIBO.log
```

### Outputs

| Output | Description |
|---|---|
| `*_coverage.txt` | binned coverage matrix |
| `*_heat_plot.pdf/png` | gene-body coverage heatmap |
| `*_coverage_bar_plot.pdf/png` | coverage barplot |
| `*_coverage_line_plot.pdf/png` | average coverage line plot |

## `rpf_Percent`

### Function

Calculate the percentage of reads assigned to different transcript regions or frames.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-t / --transcript` | normalized transcript annotation |
| `-o / --output` | output prefix |
| `-f` | reading frame |
| `-m` | minimum read count |
| `-n` | normalize to RPM |

### Example

```bash
rpf_Percent \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -n \
  -m 50 \
  -f 0 \
  -o RIBO \
  &>> RIBO.log
```

## Merge coverage results

| Command | Meaning |
|---|---|
| `merge_coverage -l *coverage.txt -o RIBO` | merge coverage tables |

## Interpretation

Uniform coverage across CDS supports reliable quantification. Strong 5' or 3' bias may indicate RNA degradation, nuclease bias, mapping artifacts, or library-specific structure.
