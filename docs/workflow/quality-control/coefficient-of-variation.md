# 4.5.15 Coefficient of Variation

Codon-level changes may reflect translation elongation, but gene-level validation can be affected by variable ribosome profiling coverage.

Low-coverage genes may appear to have stronger pausing due to coverage-dependent noise. RiboParser implements coefficient-of-variation-based analysis to model this dependence.

The model used in the README is:

```text
log2(CV) = 1/2 * log2(beta / mu + alpha)
```

where:

- `CV`: coefficient of variation in the ribosome profile of a gene
- `mu`: mean coverage, measured as RPF reads per codon
- `alpha`, `beta`: fitting parameters

When `alpha = 0` and `beta = 1`, the equation corresponds to a Poisson distribution. When `alpha > 0` and `beta = 1`, it corresponds to a negative binomial-like behavior.

## `rpf_CoV` help

```bash
rpf_CoV -h
```

```text
usage: rpf_CoV [-h] -r RPF [-g GROUP] [-l LIST] -o OUTPUT
               [-f {0,1,2,all}] [-m MIN] [-n]
               [--tis TIS] [--tts TTS] [--fig]

Required arguments:
  -r RPF       input RPF density file
  -g GROUP     sample group file
  -l LIST      gene list; default: whole
  -o OUTPUT    output prefix

Options:
  -f           reading frame
  -m MIN       minimum RPF count
  -n           normalize to RPM
  --tis        discard codons after TIS
  --tts        discard codons before TES
  --fig        draw fitted figure
```

## Group design file

```text
Name                         Group
WT_ribo_YPD1                 WT_ribo_YPD
WT_ribo_YPD2                 WT_ribo_YPD
WT_ribo_YPD3                 WT_ribo_YPD
ncs2d_ribo_YPD1              ncs2d_ribo_YPD
ncs2d_ribo_YPD2              ncs2d_ribo_YPD
ncs2d_ribo_YPD3              ncs2d_ribo_YPD
elp6d_ribo_YPD1              elp6d_ribo_YPD
elp6d_ribo_YPD2              elp6d_ribo_YPD
elp6d_ribo_YPD3              elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD1        ncs2d_elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD2        ncs2d_elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD3        ncs2d_elp6d_ribo_YPD
```

## Example

```bash
cd ./4.ribo-seq/5.riboparser/15.coefficient_of_variation/

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

## Cumulative CoV

```bash
rpf_Cumulative_CoV -h
```

```text
usage: rpf_Cumulative_CoV [-h] -r RPF [-o OUTPUT] [-l LIST]
                          [-m MIN] [-n] [-t TRIM] [-s] [-z]

Required arguments:
  -r RPF       input RPF density file
  -o OUTPUT    output prefix

Options:
  -l LIST      gene list
  -m MIN       minimum RPF count
  -n           normalize to RPM
  -t TRIM      trim transcript with specified length
  -s           split gene RPF to each TXT file
  -z           set start site to zero
```

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

```text
gene_CoV.txt
gene_compared_CoV.txt
gene_WT_ribo_YPD_vs_ncs2d_ribo_YPD_CoV_fitplot.pdf
gene_WT_ribo_YPD_vs_ncs2d_ribo_YPD_CoV_fitplot.png
gene.log
```
