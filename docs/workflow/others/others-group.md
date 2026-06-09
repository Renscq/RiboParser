## Ribo-Utility groups


## `rpf_Shuffle`

### Function

Generate shuffled density as random controls.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r` | input density file |
| `-o` | output prefix |
| `-l` | gene list |
| `-s` | random seed |
| `-i` | shuffle RPFs independently for each sample |

### Example

```bash
rpf_Shuffle \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -s 0 \
  -i \
  -o RIBO
```

## `rpf_Bam2bw`

### Function

Convert BAM-derived signal to bedGraph/bigWig-compatible tracks.

### Parameters

| Parameter | Meaning |
|---|---|
| `-b` | input BAM |
| `-o` | output prefix |
| `-p` | optional P-site offset table |
| `-m` | minimum read length |
| `-M` | maximum read length |
| `--strand` | output strand-specific signal |
| `--normalize` | normalize signal |
| `--thread` | number of threads |

## `rpf_Geneplot`

### Function

Plot Ribo-seq or RNA-seq density around selected genes.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r` | input density file |
| `-g` | gene ID or gene list |
| `-o` | output prefix |
| `-t` | transcript annotation |
| `-n` | normalize to RPM |
| `-f` | reading frame |
