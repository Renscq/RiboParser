# 4.2 Raw data download

## Purpose

Download public RNA-seq and Ribo-seq datasets, or organize user-provided FASTQ files, into the `2.rawdata` directory for downstream cleaning.

The workflow is divided into three steps:

```text
Step 1: Prepare the raw data directory
Step 2: Download RNA-seq raw data
Step 3: Download Ribo-seq raw data
```

All commands are run inside the `2.rawdata` directory created in Step 1.

## Input

| Input | Description | Required |
|---|---|---|
| SRA accession | SRR/SRX/SRP/GSE identifiers | yes (public data) |
| output directory | raw data storage directory | yes |
| SRA Toolkit | `prefetch`, `fasterq-dump` | yes (public data) |
| compression tool | `pigz` or `gzip` | yes |

## Step 1: Prepare the raw data directory

```bash
mkdir -p ./sce/2.rawdata/rna-seq/
mkdir -p ./sce/2.rawdata/ribo-seq/
cd ./sce/2.rawdata/
```

## Step 2: Download RNA-seq raw data

The following example downloads the RNA-seq samples of the public dataset used in this tutorial.

```bash
cd ./rna-seq/

for sra in SRR1944925 SRR1944926 SRR1944927 SRR1944928 SRR1944929 SRR1944930 \
           SRR1944931 SRR1944932 SRR1944933 SRR1944934 SRR1944935
do
  prefetch -o ${sra}.sra ${sra}
  fasterq-dump ${sra}.sra
  pigz ${sra}.fastq
done
```

## Step 3: Download Ribo-seq raw data

```bash
cd ./ribo-seq/

for sra in SRR1944912 SRR1944913 SRR1944914 SRR1944915 SRR1944916 SRR1944917 \
           SRR1944918 SRR1944919 SRR1944920 SRR1944921 SRR1944922 SRR1944923
do
  prefetch -o ${sra}.sra ${sra}
  fasterq-dump ${sra}.sra
  pigz ${sra}.fastq
done
```

## Notes

- `fasterq-dump` is the recommended replacement for the deprecated `fastq-dump`.
- For paired-end data, `fasterq-dump` produces `*.fastq` files with `_1`/`_2` suffixes; adjust the compression loop accordingly.
- If you already have user-provided FASTQ files, copy them directly into `2.rawdata/rna-seq/` or `2.rawdata/ribo-seq/` and skip the download steps.
- The example accessions are taken from a public dataset; replace them with the accessions of your own project.

## Output

| File | Description |
|---|---|
| `*.sra` | downloaded SRA archive |
| `*.fastq.gz` | compressed FASTQ file |
