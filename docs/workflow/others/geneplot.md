# rpf_Geneplot

`rpf_Geneplot` draws per-gene RPF or RNA-seq density profiles from a merged density table.

### Function

Use this command to visualize the nucleotide-level density profile of one gene or a list of genes across samples. The command extracts the selected gene records from a merged density table, converts frame-specific density columns into a continuous nucleotide-level profile, and writes one plot per gene.

### Input

The input file should be a merged RPF/RNA density table with annotation columns and frame-specific sample columns, for example:

```text
name  now_nt  from_tis  from_tts  region  codon  sample1_f0  sample1_f1  sample1_f2  ...
```

Gene IDs supplied with `-g` or `-l` should match the values in the `name` column.

### Parameters

| Parameter | Required | Meaning |
|---|---:|---|
| `-r` | yes | Input merged RPF/RNA density file in TXT/TSV format. |
| `-g` | no | Single gene/transcript ID to plot. Mutually exclusive with `-l`. |
| `-l` | no | Gene/transcript ID list file. Each line should contain one ID. Mutually exclusive with `-g`. |
| `--log` | no | Set y-axis to log scale. Allowed values: `2` or `10`. Default: no log scaling. |

At least one of `-g` or `-l` should be provided.

### Output

For each gene/transcript, the command writes two plot files in the current working directory:

```text
<gene_id>_density_plot.pdf
<gene_id>_density_plot.png
```

The x-axis is nucleotide position relative to the selected transcript region, and the y-axis shows normalized density used by the plotting function.

### Examples

Draw the density profile of one gene:

```bash
rpf_Geneplot \
  -r RIBO_merged.txt \
  -g AT1G01010
```

Draw one profile for each gene in a list:

```bash
rpf_Geneplot \
  -r RIBO_merged.txt \
  -l target_gene_ids.txt
```

Use log2 scaling on the y-axis:

```bash
rpf_Geneplot \
  -r RIBO_merged.txt \
  -g AT1G01010 \
  --log 2
```

Use log10 scaling on the y-axis:

```bash
rpf_Geneplot \
  -r RIBO_merged.txt \
  -l target_gene_ids.txt \
  --log 10
```

### Notes

- The command does not take an output prefix; plot file names are generated from gene/transcript IDs.
- The previous example using `-b` should be replaced by `-r`, which is the actual input-density argument.
- If a gene/transcript ID is not found in the input table, the command reports it and skips that record.
- Very long genes or many samples can generate large figures. Run the command on a small gene list first when checking a new dataset.
