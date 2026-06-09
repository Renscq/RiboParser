# 4.5.7 Metaplot

Metagene analysis investigates read density around start and stop codons.

For Ribo-seq, metagene profiles can reveal:

- initiation-associated peak/trough patterns
- termination-associated accumulation
- translation dynamics around TIS/TTS
- periodicity across metagene windows

## Command help

```bash
rpf_Metaplot -h
```

```text
usage: rpf_Metaplot [-h] -t TRANSCRIPT -r RPF -o OUTPUT
                    [-m MIN] [--utr5 UTR5] [--cds CDS] [--utr3 UTR3]
                    [-n] [--mode {line,bar}]

Required arguments:
  -t TRANSCRIPT  transcript annotation in TXT format
  -r RPF         RPF density file
  -o OUTPUT      output prefix

Options:
  -m MIN         delete transcript with less than minimum RPFs
  --utr5 UTR5    codon number in 5' UTR region
  --cds CDS      codon number in CDS region
  --utr3 UTR3    codon number in 3' UTR region
  -n             normalize RPF count to RPM
  --mode         line or bar
```

## RNA-seq metagene analysis

```bash
cd ./3.rna-seq/5.riboparser/07.metaplot/

rpf_Metaplot \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RNA_merged.txt \
  -m 50 \
  --mode bar \
  -o RNA \
  &>> RNA.log
```

## Ribo-seq metagene analysis

```bash
cd ./4.ribo-seq/5.riboparser/07.metaplot/

rpf_Metaplot \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  --mode bar \
  -o RIBO \
  &>> RIBO.log
```

## Output files

```text
RIBO_tis_tts_metaplot.txt
RIBO_SRR1944912_meta_bar_plot.pdf
RIBO_SRR1944912_meta_bar_plot.png
RIBO_SRR1944912_tis_tts_metaplot.txt
RIBO.log
```
