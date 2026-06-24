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


## `rpf_Retrieve`

### Function

Convert genome alignment BAM-derived signal to bedGraph/bigWig-compatible tracks.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r` | the name of input RPFs density file in TXT format. |
| `-o` | prefix of output file name (default: filename + '_retrieve.txt'. |
| `-l` | the list of input genes for transcript id. |
| `-m` | retain transcript with more than minimum RPFs (default: 0). |
| `-n` | normalise the RPFs to RPM. (default: False). |
| `-f` | melt three column data of each sample to one column (default: False). |
| `-s` | split gene rpf to each TXT file (default: False). |


### Example

```bash
# filter the gene with minimal rpf
rpf_Retrieve \
 -m 50 \
 -r sample1_rpf_merged.txt \
 -o sample1_rpf_50

# filter the gene 
rpf_Retrieve \
 -l gene_id \ # or gene id list in txt format
 -r sample1_rpf_merged.txt \
 -o sample1_rpf_gene_id

 # convert the wider to longer format
rpf_Retrieve \
 -m 0 \
 -f \
 -r sample1_rpf_merged.txt \
 -o sample1_rpf
```


## `rpf_Geneplot`

### Function

Plot Ribo-seq or RNA-seq density around selected genes.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r` | the name of input RPFs file in TXT format. |
| `-g` | gene ID |
| `-l` | the gene name list in TXT format. |
| `--log` | {2, 10} set the y-axis to log scaling. (default: Not).|


### Example

```bash
rpf_Geneplot \
 -b sample1_rpf_merged.txt \
 -g gene_id \
 --log 2

```
