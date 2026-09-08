# 4.9.1 SeRP peak

## `serp_peak`

`serp_peak` detects enriched SeRP/IP signal peaks from RPF density profiles by comparing immunoprecipitation samples with control samples.

### Function

Use this command when you want to:

- compare control and IP RPF density profiles;
- normalize sample read counts to RPM;
- scan smoothed IP/control enrichment ratios along transcripts;
- identify enriched binding or collision peak regions;
- export peak tables, BED regions, peak-associated sequences, ratio tables, and reusable enrichment profiles for later plotting.

### Input

The main input is a RiboParser merged RPF coverage table in TXT format, or a JSON/JSONL density file. The command expects the first annotation columns and then frame-specific sample columns.

A typical TXT table has the following structure:

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
| `-r` | yes | Input RPF density file in TXT, JSON, JSONL, or compressed JSON format. |
| `--ck` | yes | Control sample names, separated by commas, for example `CK1,CK2`. In consensus mode, the sample order defines biological replicate pairing with `--ip`. |
| `--ip` | yes | Immunoprecipitation sample names, separated by commas, for example `IP1,IP2`. In consensus mode, the sample order defines biological replicate pairing with `--ck`. |
| `-o` | yes | Output file prefix. |
| `-n` | no | Total RPF count table used for normalization. If omitted, total counts are calculated from the input table. |
| `-a` | no | Gene annotation file in TXT format. |
| `--method` | no | Peak-calling method. `legacy` uses the corrected historical algorithm; `consensus` uses matched biological replicates and peak-overlap support. Default: `legacy`. |
| `-m` | no | Minimum gene-level RPF count required in every sample for a gene to be retained. Default: `50`. |
| `--corr` | no | Minimum replicate correlation required within each group. Method defaults are `0.3` for legacy and `0.5` for consensus. |
| `--scale` | no | Normalization scale. Default: `1000000` for RPM. |
| `-f` | no | Legacy-only control zero-fill mode: `0` uses the global CDS background, `-1` uses the current-gene CDS mean, and positive values use a leading-codon window. Default: `30`. |
| `--back` | no | Leading CDS codons treated as background/non-callable region. Legacy mode supports `0` or `30`; consensus mode accepts any non-negative value. Default: `0`, meaning no 5-prime background filter is applied. |
| `--bf` | no | Legacy-only switch for maximum-background fold filtering. Defaults to `True`. |
| `-s` | no | Legacy-only Savitzky-Golay smoothing window size; use `0` to disable smoothing. The value should be a positive odd integer when smoothing is enabled. Default: `3`. |
| `-k` | no | Legacy-only Savitzky-Golay polynomial order. Default: `1`. |
| `-w` | no | Minimum binding peak width in amino acids. Default: `5`. |
| `-e` | no | Enrichment threshold for peak height. Default: `2.0`. |
| `-c` | no | Lower enrichment threshold used to bridge/extend peak edges. `--collision` is retained as a compatibility alias. Default: `1.5`. |
| `-g` | no | Maximum consecutive gap length retained within a candidate peak. Default: `1`. |
| `-p` | no | Maximum gap proportion retained within a candidate peak. Default: `0.2`. |
| `--all` | no | Legacy-only: retain all qualified peak-region permutations. By default, only the optimal non-overlapping peak set is retained. |
| `--consensus-window` | no | Centered rolling-sum window in codons used for each matched replicate pair. Default: `5`. |
| `--pseudocount` | no | RPM pseudocount added before local IP/control ratio calculation. Default: `0.1`. |
| `--min-support` | no | Minimum fraction of matched replicate pairs supporting each consensus peak. Default: `1.0`. |
| `--min-overlap` | no | Minimum replicate/consensus peak overlap in codons; `0` uses `ceil(width/2)`. Default: `0`. |
| `--stop-trim` | no | Number of terminal CDS codons excluded from consensus QC and peak calling. Default: `5`. |
| `--min-codon-rpf` | no | Consensus-only minimum mean raw RPF count per analyzed CDS codon required in every CK/IP sample; `0` disables this filter. Default: `0.0`. |
| `--max-edge-extension` | no | Consensus-only maximum lower-threshold edge extension per peak side in codons. Default: `10`. |
| `--up` | no | Number of upstream codons retrieved around each peak. Default: `10`. |
| `--down` | no | Number of downstream codons retrieved around each peak. Default: `10`. |
| `--ratio` | no | Also output method-specific replicate enrichment ratios: all-pairwise ratios for legacy mode or matched local ratios for consensus mode. |

### Output

The output prefix is controlled by `-o`.

| Output | Description |
|---|---|
| `<prefix>_peaks.log` | Peak-scanning log. |
| `<prefix>_peaks.txt` | Main peak result table with peak statistics and adjusted significance values. |
| `<prefix>_peaks.bed` | BED-like peak regions for all reported peaks. |
| `<prefix>_sig_peaks.bed` | BED-like peak regions with `P_Value < 0.05`. |
| `<prefix>_peaks_sequence.txt` | Upstream, peak, and downstream nucleotide/amino-acid sequences. |
| `<prefix>_peaks_ratio.txt` | Smoothed enrichment-ratio table with peak/collision annotations. It can be reused by `serp_plot` for figure generation. |
| `<prefix>_enrich_ratio.txt` | Method-specific replicate enrichment ratios. Written only when `--ratio` is used. |

### Examples

First, move to the working directory:

```bash
cd ./sce/6.serp/01.peak
```


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

Run the replicate-consensus caller with matched biological replicate pairs and export matched local ratios:

```bash
serp_peak \
  -r RIBO_merged.txt \
  --ck CK1,CK2 \
  --ip IP1,IP2 \
  --method consensus \
  --ratio \
  -o SeRP_consensus
```

### Notes

- `--ck` and `--ip` must use the same sample prefixes as the frame-specific columns in the RPF table.
- `--method consensus` requires equal numbers of `--ck` and `--ip` samples (at least two matched biological replicate pairs), and the sample order defines the pairing.
- The smoothing window supplied to `-s` should be `0` or a positive odd integer when smoothing is enabled.
- `--method legacy` only supports `--back 0` or `--back 30`; consensus mode accepts any non-negative `--back` value.
- `--bf` is a boolean switch that defaults to `True`, so the maximum-background fold filter is enabled by default; pass `--no-bf` to disable it.
- Figure generation is no longer part of `serp_peak`; use the independent `serp_plot` command with the reusable enrichment profiles written to `<prefix>_peaks_ratio.txt`.
