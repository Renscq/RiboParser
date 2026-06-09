# 5.5 RNA commands

| Command | Function |
|---|---|
| `rna_Offset` | Generate RNA-seq offset table with constant offset |
| `rna_Density` | Convert RNA-seq BAM to density file |

## RNA offset

```bash
rna_Offset \
  -m 27 \
  -M 50 \
  -e 12 \
  -o sample
```

## RNA density

```bash
rna_Density \
  -b sample.bam \
  -m 27 \
  -M 33 \
  -l \
  --thread 10 \
  -p sample_offset.txt \
  -s gene.norm.rna.fa \
  -t gene.norm.txt \
  -o sample
```
