# 4.7.7 Meta-codon analysis

## Function

`rpf_Meta_Codon` extracts average RPF density around selected codons or codon motifs. It can process single codons such as `AAA` or multi-codon motifs whose length is a multiple of three.

## Input files

| Input | Description |
|---|---|
| RPF density table | Merged Ribo-seq density table, usually `RIBO_merged.txt`. |
| Codon list | Optional one-column TXT file containing codons or codon motifs. If omitted, all codons from the codon table are used. |
| Gene list | Optional transcript or gene list. If not provided, all genes are used. |

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `-r` | Yes | Input RPF density table in TXT format. |
| `-l` | No | Gene name or transcript ID list in TXT format. If omitted, all genes are used. |
| `-c` | No | One-column codon or codon-motif list in TXT format. RNA `U` is converted to DNA `T`. If omitted, all codons are used. |
| `-o` | Yes | Output prefix. |
| `-f` | No | Reading frame used for meta-codon analysis. Choices in the parser are `0`, `1`, `2`; the script default is `all`. |
| `-a` | No | Number of upstream and downstream codons to retrieve around each target codon or motif. Default: `10`. |
| `-m` | No | Retain transcripts with more than this minimum number of RPFs. Default: `50`. |
| `--tis` | No | Number of codons after the translation initiation site to discard. Default: `0`. |
| `--tts` | No | Number of codons before the translation termination site to discard. Default: `0`. |
| `-n` | No | Normalize RPF counts to RPM before extracting meta-codon profiles. Disabled by default. |
| `-u` | No | Remove overlapping repeated target codons within the same window to reduce background noise. Disabled by default. |
| `-s` | No | Scale window density by gene-level density. Disabled by default. |
| `--smooth` | No | Smooth density values using a Savitzky-Golay filter. Format: `window,order`, for example `3,1`. Disabled by default. |
| `--thread` | No | Number of threads. Default: `1`. The current implementation loops over codons sequentially. |
| `--fig` | No | Draw meta-codon line plots. Disabled by default. |
| `--unit` | No | X-axis position unit, either `codon` or `nucleotide` (nt). Default: `codon`. With `--frame all`, `nucleotide` expands the density to one value per nt so that tri-nucleotide periodicity is visible in both the density table and the figure. |

## Output files

| Output | Description |
|---|---|
| `<prefix>_<codon>_<site_num>_<valid_site_num>_meta_density.txt` | Meta-density table around the selected codon or codon motif. |
| `<prefix>_<codon>_<site_num>_<valid_site_num>_meta_sequence.txt` | Sequence context table for the selected codon or codon motif. |
| `<prefix>_<codon>.pdf` / `.png` | Meta-codon line plot. Generated only with `--fig`. |

## Codon list example

```text
AAA
AAG
ATG
GCTGAA
```

## Example

```bash
rpf_Meta_Codon \
    -r ../05.merge/RIBO_merged.txt \
    -m 50 \
    -f 0 \
    -c codon_list.txt \
    -a 15 \
    -u \
    -n \
    --fig \
    -o RIBO \
    &> RIBO_meta_codon.log
```

## Notes

- The command imports A-site density internally from the merged RPF table.
- Target motifs must have lengths divisible by three; otherwise they are skipped.
- The output file name records both the total number of matched sites and the number of valid sites after optional filtering.
- Use `-u` when repeated target codons within the same window may inflate local background signal.
