# 4.5.16 Meta-codon analysis

Meta-codon analysis visualizes RPF density within a defined window around selected codons or codon motifs.

It supports frame-specific visualization and smoothing to reduce unstable spike signals.

## Command help

```bash
rpf_Meta_Codon -h
```

```text
usage: rpf_Meta_Codon [-h] [-l LIST] -r RPF [-c CODON] -o OUTPUT
                      [-f {0,1,2}] [-a AROUND] [-m MIN]
                      [--tis TIS] [--tts TTS] [-n] [-u] [-s]
                      [--smooth SMOOTH] [--thread THREAD] [--fig]

Required arguments:
  -l LIST      gene list; default: whole
  -r RPF       input RPF density file
  -c CODON     codon list
  -o OUTPUT    output prefix

Options:
  -f           reading frame
  -a AROUND    upstream/downstream codon window length
  -m MIN       minimum RPF count
  --tis        discard codons after TIS
  --tts        discard codons before TTS
  -n           normalize to RPM
  -u           remove cross-repetition codons in different windows
  -s           scale window density by gene density
  --smooth     smoothing window
  --thread     number of threads
  --fig        output figure
```

## Codon list

```text
AAA
AAC
AAG
AAT
AAGAAG
ATGATG
CCCGGG
...
```

## Example

```bash
cd ./4.ribo-seq/5.riboparser/16.meta_codon/

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

```text
RIBO_AAA_97591_8146_meta_density.txt
RIBO_AAA_97591_8146_meta_sequence.txt
RIBO_AAA.pdf
RIBO_AAA.png
RIBO.log
```
