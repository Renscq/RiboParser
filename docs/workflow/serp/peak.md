# 4.9.1 SeRP peak

## `serp_peak`

`serp_peak` detects enriched SeRP/IP signal peaks from a merged RPF coverage table by comparing immunoprecipitation samples with control samples.

### Function

Use this command when you want to:

- compare control and IP RPF density profiles;
- normalize sample read counts to RPM;
- scan smoothed IP/control enrichment ratios along transcripts;
- identify enriched binding or collision peak regions;
- export peak tables, BED regions, peak-associated sequences, ratio tables, and optional figures.

### Input

The main input should be a RiboParser merged RPF coverage table. The command expects the first annotation columns and then frame-specific sample columns.

A typical table has the following structure:

```text
name  now_nt  from_tis  from_tts  region  codon  CK1_f0  CK1_f1  CK1_f2  IP1_f0  IP1_f1  IP1_f2  ...
```

The sample names supplied to `--ck` and `--ip` should match the sample prefixes before `_f0`, `_f1`, and `_f2`.

Optional inputs:

| Input | Description |
|---|---|
| Normalization table | Two-column text file containing sample name and total RPF count. If omitted, total counts are calculated from the input RPF table. |
| Gene annotation file | Optional annotation file used to add gene names to transcript IDs. |

### Parameters

| Parameter | Required | Meaning |
|---|---:|---|
| `-r` | yes | Input RPF coverage table. |
| `-n` | no | Total RPF count table used for normalization. If omitted, total counts are calculated from the input table. |
| `--scale` | no | Normalization scale. Default: `1000000` for RPM. |
| `-a` | no | Gene annotation file in TXT format. |
| `--ck` | yes | Control sample names, separated by commas, for example `CK1,CK2`. |
| `--ip` | yes | Immunoprecipitation sample names, separated by commas, for example `IP1,IP2`. |
| `-m` | no | Minimum gene-level RPF coverage required for a gene to be retained. Default: `50`. |
| `--corr` | no | Minimum replicate correlation required within each group. Default: `0.3`. |
| `-f` | no | Missing-value filling mode. Default: `30`, meaning the first 30 amino-acid positions are used as the background region. Use `-1` for the current-gene mean, or `0` for the global gene mean. |
| `-s` | no | Savitzky-Golay smoothing window size. The value should be odd. Default: `3`. |
| `-k` | no | Savitzky-Golay polynomial order. Default: `1`. |
| `-w` | no | Minimum binding peak width in amino acids. Default: `5`. |
| `-e` | no | Enrichment threshold for peak height. Default: `2.0`. |
| `-c` | no | Enrichment threshold used for collision/edge extension around a peak. Default: `1.5`. |
| `-g` | no | Maximum allowed internal gap width inside a peak. Default: `1`. |
| `-p` | no | Maximum allowed gap proportion inside a peak. Default: `0.2`. |
| `--back` | no | Background-region mode for filtering enrichment. Default: `0`, meaning no 5-prime background filter is applied. |
| `--bf` | no | Keep the background-fold filter enabled. In the current implementation this option defaults to `True`. |
| `--up` | no | Number of upstream codons retrieved around each peak. Default: `10`. |
| `--down` | no | Number of downstream codons retrieved around each peak. Default: `10`. |
| `-o` | no | Output prefix. Default: `results`. |
| `--all` | no | Output all peak regions, including overlapping candidates. By default, only the optimal non-overlapping peak is retained. |
| `--rpm` | no | Registered command-line option for peak RPM output. The current implementation does not write a separate peak RPM file because that output block is commented out. |
| `--ratio` | no | Output original enrichment ratios for each gene. |
| `--fig` | no | Draw demo figures for peak scanning results. This can take a long time for large datasets. |

### Output

The output prefix is controlled by `-o`.

| Output | Description |
|---|---|
| `<prefix>_peaks.log` | Peak-scanning log. |
| `<prefix>_peaks.txt` | Main peak result table with peak statistics and adjusted significance values. |
| `<prefix>_peaks.bed` | BED-like peak regions for all reported peaks. |
| `<prefix>_sig_peaks.bed` | BED-like peak regions with `P_Value < 0.05`. |
| `<prefix>_peaks_sequence.txt` | Upstream, peak, and downstream nucleotide/amino-acid sequences. |
| `<prefix>_peaks_ratio.txt` | Smoothed enrichment-ratio table and peak/collision annotations. |
| `<prefix>_enrich_ratio.txt` | Original per-gene ratio output. Written only when `--ratio` is used. |
| `<prefix>_peaks_scripts.m` | MATLAB script output for peak regions. |
| `<prefix>_peaks_none_scripts.m` | MATLAB script output for genes without peaks or filtered genes. |
| `<prefix>_figures/` | Optional figure directory. Written only when `--fig` is used. |

### Examples

Run peak detection with two control and two IP samples:

```bash
serp_peak \
  -r RIBO_merged.txt \
  --ck CK1,CK2 \
  --ip IP1,IP2 \
  -a gene_annotation.txt \
  -o SeRP_IP
```

Run peak detection with a custom normalization table:

```bash
serp_peak \
  -r RIBO_merged.txt \
  -n total_rpf_counts.txt \
  --ck CK1,CK2 \
  --ip IP1,IP2 \
  -m 50 \
  --corr 0.3 \
  -o SeRP_IP.norm
```

Use stricter enrichment and peak-width filters:

```bash
serp_peak \
  -r RIBO_merged.txt \
  --ck Mock1,Mock2 \
  --ip Flag1,Flag2 \
  -e 2.5 \
  -w 8 \
  -g 1 \
  -p 0.2 \
  -o SeRP_strict
```

Export original ratio values and demo figures:

```bash
serp_peak \
  -r RIBO_merged.txt \
  --ck CK1,CK2 \
  --ip IP1,IP2 \
  --ratio \
  --fig \
  -o SeRP_with_figures
```

### Notes

- `--ck` and `--ip` must use the same sample prefixes as the frame-specific columns in the RPF table.
- The smoothing window supplied to `-s` should be an odd integer when smoothing is enabled.
- `--fig` can be slow because it draws per-gene peak-scanning figures.
- The `--rpm` option is present in the command-line parser, but the separate peak-RPM output section is currently commented out in the implementation.
- `--bf` currently defaults to `True`, so passing the flag does not switch the background-fold filter from `False` to `True`; it is already enabled.
