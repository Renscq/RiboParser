## `rpf_Retrieve`

### Function

Retrieve the gene rpf density from merged rpf density file.

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
