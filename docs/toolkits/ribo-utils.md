# 5.4 Ribo-utils commands

## Data shuffling

Some analysis processes require randomly assigned data as controls.

```bash
rpf_Shuffle -h
```

```text
usage: rpf_Shuffle [-h] -r RPF -o OUTPUT [-l LIST] [-s SEED] [-i]

Required arguments:
  -r RPF       input RPF density file
  -o OUTPUT    output prefix

Options:
  -l LIST      gene list
  -s SEED      random seed
  -i           shuffle RPFs for each sample
```

Ribo-seq example:

```bash
cd ./sce/4.ribo-seq/5.riboparser/17.shuffle/

rpf_Shuffle \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -s 0 \
  -i \
  -o RIBO \
  &>> RIBO.log
```

RNA-seq example:

```bash
cd ./sce/3.rna-seq/5.riboparser/10.shuffle/

rpf_Shuffle \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RNA_merged.txt \
  -s 0 \
  -i \
  -o RNA \
  &>> RNA.log
```

## Retrieve and format gene density

```bash
rpf_Retrieve -h
```

```text
usage: rpf_Retrieve [-h] -r RPF [-o OUTPUT] [-l LIST] [-m MIN] [-n] [-f] [-s]

Required arguments:
  -r RPF       input RPF density file
  -o OUTPUT    output prefix

Options:
  -l LIST      gene list
  -m MIN       retain transcript with more than minimum RPFs
  -n           normalize RPF count to RPM
  -f           melt three-column data of each sample to one column
  -s           split gene RPF to each TXT file
```

## Frame-shift detection

A frameshift can be detected by ribosome occupancy in different reading frames.

```bash
rpf_Shift -h
```

```text
usage: rpf_Shift.py [-h] -r RPF -o OUTPUT [-t TRANSCRIPT] [-p PERIOD]
                    [-m MIN] [--tis TIS] [--tts TTS]

Required arguments:
  -r RPF         input RPF density file
  -o OUTPUT      output prefix
  -t TRANSCRIPT  transcript annotation in TXT format

Options:
  -p PERIOD      minimum in-frame value for frame-shift screening
  -m MIN         minimum RPF count
  --tis          discard codons after TIS
  --tts          discard codons before TTS
```

Example:

```bash
cd ./sce/4.ribo-seq/5.riboparser/19.frame_shift/

rpf_Shift \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  --tis 5 \
  --tts 5 \
  -m 50 \
  -p 45 \
  -o RIBO \
  &>> RIBO.log
```

Output:

```text
RIBO_gene_frame_shift_count_plot.pdf
RIBO_gene_frame_shift_count_plot.png
RIBO_gene_frame_shift_count.txt
RIBO_gene_periodicity.txt
RIBO_SRR1944912_gene_frame_shift.txt
RIBO.log
```
