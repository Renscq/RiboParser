# 4.2 Raw data download

The tutorial uses `GSE67387` as the example dataset.

## RNA-seq raw data

Use `prefetch` from `sra-tools` to download SRA-format data and convert it to FASTQ.

```bash
mkdir -p ./sce/2.rawdata/rna-seq/
cd ./sce/2.rawdata/rna-seq/

prefetch -o SRR1944925.sra SRR1944925
prefetch -o SRR1944926.sra SRR1944926
prefetch -o SRR1944927.sra SRR1944927
prefetch -o SRR1944928.sra SRR1944928
prefetch -o SRR1944929.sra SRR1944929
prefetch -o SRR1944930.sra SRR1944930
prefetch -o SRR1944931.sra SRR1944931
prefetch -o SRR1944932.sra SRR1944932
prefetch -o SRR1944933.sra SRR1944933
prefetch -o SRR1944934.sra SRR1944934
prefetch -o SRR1944935.sra SRR1944935

# decompression
for sra in *.sra
do
  fastq-dump $sra
  pigz *fastq
done
```

## Ribo-seq raw data

```bash
mkdir -p ./sce/2.rawdata/ribo-seq/
cd ./sce/2.rawdata/ribo-seq/

prefetch -o SRR1944912.sra SRR1944912
prefetch -o SRR1944913.sra SRR1944913
prefetch -o SRR1944914.sra SRR1944914
prefetch -o SRR1944915.sra SRR1944915
prefetch -o SRR1944916.sra SRR1944916
prefetch -o SRR1944917.sra SRR1944917
prefetch -o SRR1944918.sra SRR1944918
prefetch -o SRR1944919.sra SRR1944919
prefetch -o SRR1944920.sra SRR1944920
prefetch -o SRR1944921.sra SRR1944921
prefetch -o SRR1944922.sra SRR1944922
prefetch -o SRR1944923.sra SRR1944923

# decompression
for sra in *.sra
do
  fastq-dump $sra
  pigz *fastq
done
```

## Suggested directory

```text
2.rawdata/
├── rna-seq/
└── ribo-seq/
```
