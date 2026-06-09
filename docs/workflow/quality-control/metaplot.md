# 4.5.7 Metaplot

## Purpose

`rpf_Metaplot` generates metagene profiles around TIS and TTS regions to evaluate initiation/termination patterns and positional density distribution.

## Input files

| Input | Description |
|---|---|
| density file | merged RNA/Ribo density |
| transcript annotation | `gene.norm.txt` |
| minimum read count | threshold for included transcripts |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-t / --transcript` | normalized transcript annotation |
| `-o / --output` | output prefix |
| `-m / --min` | minimum read count |
| `--utr5` | codon/bin number in 5' UTR window |
| `--cds` | codon/bin number in CDS window |
| `--utr3` | codon/bin number in 3' UTR window |
| `-n` | normalize read count to RPM |
| `--mode` | plot mode: `line` or `bar` |

## Example

```bash
rpf_Metaplot \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  --mode bar \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_tis_tts_metaplot.txt` | merged TIS/TTS metaplot table |
| `RIBO_*_meta_bar_plot.pdf/png` | sample-level metaplot figures |
| `RIBO_*_tis_tts_metaplot.txt` | sample-level metaplot tables |
| `RIBO.log` | running log |

## Result interpretation

Ribo-seq metagene profiles should show interpretable density around TIS/TTS. Strong artifacts at ends may indicate nuclease bias, offset issues, or low-quality libraries.

## Merge related results

### `merge_metagene`

| Parameter | Meaning |
|---|---|
| `-l` | metagene result files |
| `-o` | output prefix |

```bash
merge_metagene -l *metaplot.txt -o RIBO
```
