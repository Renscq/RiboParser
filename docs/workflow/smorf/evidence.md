# 4.8.3 smORF evidence

## Purpose

`smorf_evidence` evaluates translation evidence for smORFs using strand-aware P-site density from bedGraph/WIG files.

## Input files

| Input | Description |
|---|---|
| ORF table | filtered smORF table from `smorf_filter` |
| genePred annotation | optional ORF exon-block annotation |
| density files | plus/minus or unstranded P-site density files |
| density list | sample-strand-path table |

## Parameters

| Parameter | Meaning |
|---|---|
| `-i / --orf-table` | filtered ORF table |
| `-o / --output` | output evidence table |
| `--genepred` | genePred annotation for ORF exon blocks |
| `--density-list` | TSV file with sample, strand, path, and format |
| `--density-plus` | plus-strand P-site density file |
| `--density-minus` | minus-strand P-site density file |
| `--density` | unstranded P-site density file |
| `--density-format` | `auto`, `wig`, or `bedgraph` |
| `--coord-mode` | `0based-half-open` or `1based-closed` |
| `--min-rpf-sum` | minimum total RPF count |
| `--min-covered-codon` | minimum covered codon number |
| `--min-coverage-ratio` | minimum codon coverage ratio |
| `--strong-periodicity` | strong periodicity threshold |
| `--moderate-periodicity` | moderate periodicity threshold |
| `--strong-start-pause` | strong start-pause threshold |
| `--moderate-start-pause` | moderate start-pause threshold |
| `--strong-stop-pause` | strong stop-pause threshold |
| `--moderate-stop-pause` | moderate stop-pause threshold |
| `--strong-release` | strong release threshold |
| `--moderate-release` | moderate release threshold |
| `--uniform-coverage-ratio` | coverage ratio for uniform-support label |
| `--uniform-gini` | Gini threshold for uniform density |
| `--uniform-max-to-mean` | max/mean threshold for uniform density |
| `--skewed-max-to-mean` | max/mean threshold for skewed density |
| `--skewed-top-fraction` | top-position fraction threshold |
| `--disperse-coverage-ratio` | coverage threshold for dispersed density |

## Example

```bash
smorf_evidence \
  -i mine.reliable.passed.message.txt \
  --genepred mine.genePred \
  --density-list ribo.bedgraph.list \
  -o mine.smorf.riboseq_evidence.txt
```

File `ribo.bedgraph.list` contains these message:

```bash
sample	strand	path	format
sample1 +	/project/ribo/bedgraph/sample1_plus.rpf.bedgraph	bedgraph
sample1 -	/project/ribo/bedgraph/sample1_minus.rpf.bedgraph	bedgraph
sample2 +	/project/ribo/bedgraph/sample2_plus.rpf.bedgraph	bedgraph
sample2 -	/project/ribo/bedgraph/sample2_minus.rpf.bedgraph	bedgraph
...

```

These `bedgraph` file are contain the rpf density can be generated with `rpf_Bam2Bw` function (see the 5. Other toolkits).


## Output files

| Output | Description |
|---|---|
| `mine.smorf.riboseq_evidence.txt` | long-format ORF-by-sample evidence table |
| `log file` | running log |

## Result interpretation

High-confidence translated smORFs should show sufficient RPF density, clear translational initiation and termination signatures, robust codon coverage, strong frame specificity, and reproducibility across samples. Treat single-evidence candidates cautiously.
