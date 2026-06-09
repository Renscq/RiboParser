# 4.8.4 smORF integrate

## Purpose

`smorf_integrate` integrates long-format smORF evidence into ORF-level matrices and summaries.

## Input files

| Input | Description |
|---|---|
| evidence table | long-format output from `smorf_evidence` |
| capture labels | labels defining captured/supported ORFs |
| pass labels | labels defining reliable ORFs |
| sample threshold | minimum supporting sample number |

## Parameters

| Parameter | Meaning |
|---|---|
| `-i / --input` | long-format smORF evidence table |
| `--output-matrix` | ORF-by-sample matrix output |
| `--output-integrated` | integrated ORF-level table |
| `--capture-labels` | evidence labels used to define captured samples |
| `--pass-labels` | evidence labels used to define reliable samples |
| `--excellent-min-samples` | minimum sample count for excellent support |

## Example

```bash
smorf_integrate \
  -i mine.smorf.riboseq_evidence.txt \
  --output-matrix mine.smorf.frame_density_matrix.txt \
  --output-integrated mine.smorf.integrated_evidence.txt
```

## Output files

| Output | Description |
|---|---|
| `mine.smorf.frame_density_matrix.txt` | ORF-by-sample evidence matrix |
| `mine.smorf.integrated_evidence.txt` | integrated ORF-level evidence table |

## Result interpretation

Use the integrated evidence table to prioritize smORFs with multi-sample, frame-specific, and robust translation evidence.
