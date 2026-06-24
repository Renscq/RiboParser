# 4.7.5 Coefficient of Variation

## Function

The coefficient-of-variation modules evaluate variability in RPF density across CDS regions. RiboParser provides two related commands:

- `rpf_CoV` calculates gene-level CoV values from CDS RPF density and can compare CoV distributions between sample groups.
- `rpf_Cumulative_CoV` calculates cumulative CoV along CDS positions and can output a meta table for positional visualization.

The model described for coverage-dependent CoV is:

```text
log2(CV) = 1/2 * log2(beta / mu + alpha)
```

where `CV` is the coefficient of variation, `mu` is the mean coverage, and `alpha` / `beta` are fitting parameters.

## 4.7.5.1 rpf_CoV

### Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| Group design file | Optional two-column sample group file. Required only for group comparison. |
| Gene list | Optional transcript or gene list. If not provided, all genes are used. |

### Parameters

| Parameter | Required | Description |
|---|---|---|
| `-r` | Yes | Input RPF density table in TXT format. |
| `-g` | No | Sample group design file. If provided, group-level CoV comparison is performed. |
| `-l` | No | Gene name or transcript ID list in TXT format. If omitted, all genes are used. |
| `-o` | Yes | Output prefix. The main output is `<prefix>_CoV.txt`. |
| `-f` | No | Reading frame used for CoV calculation. Choices: `0`, `1`, `2`, `all`. Default: `all`. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `5`. |
| `-n` | No | Normalize RPF counts to RPM before CoV calculation. Disabled by default. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `15`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `5`. |
| `--fig` | No | Draw CoV fitting plots. Disabled by default. |

### Group design file

The group design file should contain at least two columns: `Name` and `Group`.

```text
Name	Group
WT_ribo_YPD1	WT_ribo_YPD
WT_ribo_YPD2	WT_ribo_YPD
WT_ribo_YPD3	WT_ribo_YPD
ncs2d_ribo_YPD1	ncs2d_ribo_YPD
ncs2d_ribo_YPD2	ncs2d_ribo_YPD
ncs2d_ribo_YPD3	ncs2d_ribo_YPD
```

### Output files

| Output | Description |
|---|---|
| `<prefix>_CoV.txt` | Gene-level RPF sum, mean, and CoV table. |
| `<prefix>_compared_CoV.txt` | Group comparison statistics. Generated only when `-g` is provided. |
| `<prefix>_<group1>_vs_<group2>_CoV_fitplot.pdf` / `.png` | CoV fitting plot for each group comparison. |

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
    &> RIBO_CoV.log
```

## 4.7.5.2 rpf_Cumulative_CoV

### Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| Gene list | Optional transcript or gene list. If not provided, all genes are used. |

### Parameters

| Parameter | Required | Description |
|---|---|---|
| `-r` | Yes | Input RPF density table in TXT format. |
| `-o` | No | Output prefix. If omitted, the command uses the input file name-derived prefix. |
| `-l` | No | Gene name or transcript ID list in TXT format. If omitted, all genes are used. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `0`. |
| `-n` | No | Normalize RPF counts to RPM before cumulative CoV calculation. Disabled by default. |
| `-t` | No | Trim transcripts to this nucleotide length for the meta cumulative CoV table. Default: `50`. |
| `-s` | No | Split gene-level cumulative CoV results into one file per gene. Disabled by default. |
| `-z` | No | Reset each CDS start position to zero-based or CDS-start-relative coordinates before output. Disabled by default. |

### Output files

| Output | Description |
|---|---|
| `<prefix>_cumulative_CoV.txt` | Cumulative CoV table for all retained genes. |
| `<gene>_cumulative_CoV.txt` | Per-gene cumulative CoV table. Generated only with `-s`. |
| `<prefix>_meta<trim>_cumulative_CoV.txt` | Meta cumulative CoV table for the first `trim` nucleotides. |

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
    &> RIBO_cumulative_CoV.log
```

## Notes

- Low-coverage transcripts often show inflated CoV. Choose `-m` carefully according to the sequencing depth.
- The `rpf_CoV` command uses P-site density internally for CoV calculation.
- `rpf_Cumulative_CoV` imports all frames from the P-site density and then reshapes the density table into nucleotide-level cumulative CoV profiles.
