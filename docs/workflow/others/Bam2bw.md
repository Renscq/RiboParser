## `rpf_Bam2bw`

### Function

Convert genome alignment BAM-derived signal to bedGraph/bigWig-compatible tracks.

### Parameters

| Parameter | Meaning |
|---|---|
| `-b` | input genome alignment BAM file |
| `-p` | p-site offset file in TXT format. (default: None). |
| `-o` | output prefix |
| `-f` | {bedgraph,wig}  output the file format. (default: bedgraph).|
| `-n` | normalise the RPFs to RPM. (default: False). |
| `-m` | merge minus and plus strand to one file. (default: False). |
| `-t` | set the reads multiple aligned times. (default: 3). |
| `--second` | discard the reads aligned secondary loci. (default: False). |
| `--supply` | discard the supplementary reads (default: False ).|



### Example

```bash
rpf_Bam2bw \
 -b sample1_Aligned.sortedByCoord.out.bam \
 -p sample1_RSBM_offset.txt \
 -t 1 \
 --second \
 --supply \
 -f bedgraph \
 -n \
 -o sample1

```
