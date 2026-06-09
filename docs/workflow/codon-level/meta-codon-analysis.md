# 4.7.6 Meta-codon analysis

## Purpose

`rpf_Meta_Codon` visualizes RPF density around selected codons or codon motifs.

## Input files

| Input | Description |
|---|---|
| density file | `RIBO_merged.txt` |
| codon list | single codons or codon motifs |
| gene list | optional selected genes |

## Parameters

| Parameter | Meaning |
|---|---|
| `-r / --rpf` | input density file |
| `-l / --list` | gene list |
| `-c / --codon` | codon list file |
| `-o / --output` | output prefix |
| `-f` | reading frame |
| `-a / --around` | upstream/downstream codon window length |
| `-m` | minimum RPF count |
| `--tis` | discard codons after TIS |
| `--tts` | discard codons before TTS |
| `-n` | normalize to RPM |
| `-u` | remove cross-repetition codons across windows |
| `-s` | scale window density by gene density |
| `--smooth` | smoothing window size |
| `--thread` | number of threads |
| `--fig` | output figure |

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
  -o RIBO \
  &>> RIBO.log
```

## Output files

| Output | Description |
|---|---|
| `RIBO_AAA_*_meta_density.txt` | meta-density around codon AAA |
| `RIBO_AAA_*_meta_sequence.txt` | sequence context table |
| `RIBO_AAA.pdf/png` | meta-codon plot |
| `RIBO.log` | running log |

## Result interpretation

Meta-codon analysis is useful for evaluating local density changes around codons or motifs. Use normalized and scaled profiles for cross-sample comparison.
