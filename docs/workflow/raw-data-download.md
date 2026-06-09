# 4.2 Raw data download

Raw data can be downloaded from public databases such as NCBI SRA.

## Example

```bash
mkdir -p ./2.rawdata/ribo-seq/
cd ./2.rawdata/ribo-seq/

prefetch -o SRR1944912.sra SRR1944912
fastq-dump SRR1944912.sra
pigz SRR1944912.fastq
```

## Suggested structure

```text
2.rawdata/
├── rna-seq/
└── ribo-seq/
```

For large projects, use a sample table and SLURM array jobs.
