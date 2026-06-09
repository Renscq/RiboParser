# 4.5.2 Enzymatic bias

```bash
rpf_Digest -h
```

```text
usage: rpf_Digest [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -o OUTPUT [-l] [--scale] [-m MIN] [-M MAX]
Required arguments:
  -t TRANSCRIPT    transcript annotation in TXT format
  -s SEQUENCE      transcript sequence in FASTA format
  -b BAM           mapping file in BAM format
  -o OUTPUT        output prefix
Options:
  -l               only retain transcript with longest CDS of each gene
  --scale          scale the motif matrix
  -m MIN           minimum read length to keep
  -M MAX           maximum read length to keep
```

```bash
rpf_Digest -b sample.bam -m 27 -M 33 --scale -s gene.norm.rna.fa -t gene.norm.txt -o sample
merge_digestion -l *pwm.txt -o RIBO
```

```text
usage: merge_digestion [-h] -l LIST [LIST ...] -o OUTPUT
```
