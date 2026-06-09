# 4.8.2 smORF filter

## Purpose

`smorf_filter` filters scanned ORFs by length, category, strand, start codon, and Kozak sequence/PWM evidence.

## Input files

| Input | Description |
|---|---|
| scanner table | `*.message.txt` from `smorf_scanner` |
| Kozak source | annotated ORFs, built-in species consensus, external PWM, or sequence |
| filter settings | ORF category, length, start codon, and strand options |

## Parameters

| Parameter | Meaning |
|---|---|
| `-i / --input` | input scanner message table |
| `-o / --out-prefix` | output prefix |
| `--keep-start-codons` | allowed start codons |
| `--min-aa` | minimum ORF length |
| `--max-aa` | maximum ORF length |
| `--keep-categories` | ORF categories to keep |
| `--remove-categories` | ORF categories to remove |
| `--keep-antisense` | retain antisense ORFs |
| `--keep-secondary` | retain ORFs using secondary start codons |
| `--keep-partial` | retain incomplete or partial ORFs |
| `--kozak-mode` | `none`, `annotated`, `builtin`, `pwm`, or `sequence` |
| `--builtin-kozak` | built-in Kozak model, such as `yeast`, `plant`, `rice`, etc. |
| `--kozak-pwm` | external Kozak PWM file |
| `--kozak-seq` | external Kozak sequence motif |
| `--annotated-categories` | categories used to derive annotated Kozak background |
| `--min-annotated-kozak` | minimum annotated ORF count for deriving Kozak PWM |
| `--fallback-builtin-kozak` | fallback built-in Kozak model |
| `--no-kozak-fallback` | disable fallback Kozak model |
| `--min-kozak-pwm-score` | minimum Kozak PWM score |
| `--export-kozak-pwm` | export generated Kozak PWM |
| `--list-builtin-kozak` | list built-in Kozak models |

## Example

```bash
smorf_filter \
  -i mine.message.txt \
  -o mine.reliable \
  --min-aa 8 \
  --max-aa 10000 \
  --kozak-mode annotated \
  --keep-categories uORF,dORF,lncORF,overlap_uORF,overlap_dORF
```

## Output files

| Output | Description |
|---|---|
| `mine.reliable.passed.message.txt` | filtered ORFs that passed all criteria |
| `mine.reliable.failed.message.txt` | filtered-out ORFs and reasons |
| `Kozak PWM output` | optional exported PWM |
| `log file` | running log |

## Result interpretation

Use conservative filtering for functional smORF discovery. Keep candidate categories aligned with the biological question.
