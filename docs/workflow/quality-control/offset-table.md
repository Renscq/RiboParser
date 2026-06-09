# 4.5.3 Offset table

## RNA-seq: `rna_Offset`

```bash
rna_Offset -h
```

```text
usage: rna_Offset [-h] -o OUTPUT [-m MIN] [-M MAX] [-e OFFSET]
Required arguments:
  -o OUTPUT    output prefix
Options:
  -m MIN       minimum read length
  -M MAX       maximum read length
  -e OFFSET    constant offset assigned to selected read lengths
```

```bash
rna_Offset -m 27 -M 50 -e 12 -o sample
```

## Ribo-seq: `rpf_Offset`

```bash
rpf_Offset -h
```

```text
usage: rpf_Offset [-h] -t TRANSCRIPT -b BAM -o OUTPUT [--mode {SSCBM,RSBM}] [-a {both,tis,tts}] [-l] [-m MIN] [-M MAX] [-p EXP_PEAK] [-s SHIFT] [--silence] [-d]
```

```bash
rpf_Offset -b sample.bam -m 27 -M 33 -p 30 -d -t gene.norm.txt -o sample
```

## Merge offset results

```bash
merge_offset_detail -h
merge_offset -h
```

```text
usage: merge_offset_detail [-h] -l LIST [LIST ...] -o OUTPUT
usage: merge_offset [-h] -l LIST [LIST ...] -o OUTPUT
```

```bash
merge_offset_detail -l *end.txt -o RIBO
merge_offset -l *SSCBM_offset.txt -o RIBO_SSCBM
merge_offset -l *RSBM_offset.txt -o RIBO_RSBM
```
