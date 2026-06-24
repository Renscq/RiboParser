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
