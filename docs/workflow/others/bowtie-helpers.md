## Bowtie/RSEM helpers

Bowtie/RSEM helpers summarize alignment logs and merge quantification tables from external tools.

### Command summary

| Command | Function | Main output |
|---|---|---|
| `merge_bwt_log` | Merge Bowtie mapping logs and calculate read counts/ratios for sequential databases. | `<prefix>_mapping.txt`. |
| `merge_rsem` | Merge selected columns from RSEM result files. | Expression matrix. |

## Bowtie helper

### `merge_bwt_log`

#### Function

Summarize reads mapped to multiple sequential Bowtie databases, such as rRNA, tRNA, ncRNA, mRNA, and genome indexes. The command parses Bowtie log files, extracts the total number of reads, the number of reads reported as mapped to each database, and calculates the remaining reads as `Others`.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-l`, `--list` | Bowtie mapping log files. Shell wildcards such as `*log` can be used. |
| `-o` | Output prefix. The final file is `<prefix>_mapping.txt`. |
| `-n`, `--name` | Comma-separated database names in the same order as the mapping logs were generated. Default: `rRNA,tRNA,ncRNA,mRNA,Genome`. |

#### Example

```bash
merge_bwt_log \
  -l sample1.log sample2.log sample3.log \
  -n rRNA,tRNA,ncRNA,mRNA,Genome \
  -o RIBO
```

Output file:

```text
RIBO_mapping.txt
```

Output columns:

```text
Sample    Database    Count    Ratio
```

#### Typical use case

When reads are sequentially aligned to contaminant and target databases, the helper can be used to build a mapping summary table:

```bash
merge_bwt_log \
  -l ../2.bowtie/*/*.log \
  -n rRNA,tRNA,ncRNA,mRNA,Genome \
  -o bowtie
```

### Notes for `merge_bwt_log`

- The database names supplied with `-n` must match the order represented in each Bowtie log file.
- `Others` is calculated as `Total - sum(mapped reads in listed databases)`.
- The sample name is inferred from the log file basename before the first period.

## RSEM helper

### `merge_rsem`

#### Function

Merge one selected expression column from multiple RSEM result files into one matrix. Supported columns are `expected_count`, `TPM`, and `FPKM`.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-l`, `--list` | RSEM result files, for example `*.genes.results` or `*.isoforms.results`. |
| `-o` | Output merged matrix file. |
| `-c`, `--column` | Column to merge: `expected_count`, `TPM`, or `FPKM`. Default: `expected_count`. |

#### Examples

Merge expected counts:

```bash
merge_rsem \
  -l *.genes.results \
  -c expected_count \
  -o gene.expected_count.txt
```

Merge TPM values:

```bash
merge_rsem \
  -l *.genes.results \
  -c TPM \
  -o gene.TPM.txt
```

Merge FPKM values:

```bash
merge_rsem \
  -l *.genes.results \
  -c FPKM \
  -o gene.FPKM.txt
```

### Notes for `merge_rsem`

- Input files are read as tab-delimited tables with the first column as the row index.
- The output matrix keeps the RSEM target IDs as rows.
- Column names in the output are inferred from the input file names. Rename RSEM result files before merging if cleaner sample names are needed.
