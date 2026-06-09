# 4.5.4 Read density

## RNA-seq: `rna_Density`

```bash
rna_Density -h
```

```text
usage: rna_Density [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -p PSITE -o OUTPUT [-l] [-m MIN] [-M MAX] [--thread THREAD] [--silence]
Required arguments:
  -t TRANSCRIPT        transcript annotation in TXT format
  -s SEQUENCE          transcript sequence in FASTA format
  -b BAM               mapping file in BAM format
  -p PSITE             RNA offset table
  -o OUTPUT            output prefix
Options:
  -l                   only retain transcript with longest CDS of each gene
  -m MIN               minimum read length to keep
  -M MAX               maximum read length to keep
  --thread THREAD      number of threads
  --silence            suppress verbose output
```

```bash
rna_Density -b sample.bam -m 27 -M 33 -l --thread 10 -p sample_offset.txt -s gene.norm.rna.fa -t gene.norm.txt -o sample
```

## Ribo-seq: `rpf_Density`

```bash
rpf_Density -h
```

```text
usage: rpf_Density [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -p PSITE -o OUTPUT [-l] [-m MIN] [-M MAX] [--period PERIODICITY] [--silence] [--thread THREAD]
```

```bash
rpf_Density -b sample.bam -m 27 -M 33 --period 40 -l --thread 12 -p sample_SSCBM_offset.txt -s gene.norm.rna.fa -t gene.norm.txt -o sample
```
