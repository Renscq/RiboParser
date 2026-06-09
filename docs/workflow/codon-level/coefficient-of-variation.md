# 4.7.5 Coefficient of Variation

## Purpose

Coefficient of Variation analysis models gene-level RPF variability as a function of coverage.

The model used in the original README is:

```text
log2(CV) = 1/2 * log2(beta / mu + alpha)
```

where:

| Symbol | Meaning |
|---|---|
| `CV` | coefficient of variation in the ribosome profile of a gene |
| `mu` | mean coverage, measured as RPF reads per codon |
| `alpha`, `beta` | fitting parameters |

## `rpf_CoV`

### Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-g / --group` | sample group design file |
| `-l / --list` | gene list or annotation |
| `-o / --output` | output prefix |
| `-f` | reading frame |
| `-m` | minimum RPF count |
| `-n` | normalize to RPM |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `--fig` | draw fitted figure |

### Example

```bash
rpf_CoV \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -f 0 \
  -m 30 \
  --tis 10 \
  --tts 5 \
  --fig \
  -g design.txt \
  -o RIBO \
  &>> RIBO.log
```

### Group design example

```text
Name                         Group
WT_ribo_YPD1                 WT_ribo_YPD
WT_ribo_YPD2                 WT_ribo_YPD
WT_ribo_YPD3                 WT_ribo_YPD
ncs2d_ribo_YPD1              ncs2d_ribo_YPD
ncs2d_ribo_YPD2              ncs2d_ribo_YPD
ncs2d_ribo_YPD3              ncs2d_ribo_YPD
```

## `rpf_Cumulative_CoV`

### Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-o / --output` | output prefix |
| `-l / --list` | gene list |
| `-m` | minimum RPF count |
| `-n` | normalize to RPM |
| `-t / --trim` | trim transcript to specific length |
| `-s` | split gene density into individual files |
| `-z` | set start site to zero |

### Example

```bash
rpf_Cumulative_CoV \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  -n \
  -z \
  -t 300 \
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `gene_CoV.txt` | gene-level CoV table |
| `gene_compared_CoV.txt` | group-comparison CoV table |
| `*_CoV_fitplot.pdf/png` | model fitting plot |
| `RIBO.log` | running log |

## Interpretation

Low-coverage genes often show inflated variability. CoV modeling helps distinguish real translational heterogeneity from coverage-dependent noise.
