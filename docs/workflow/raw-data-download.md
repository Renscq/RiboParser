# 4.2 Raw data download

```bash
mkdir -p ./sce/2.rawdata/rna-seq/
cd ./sce/2.rawdata/rna-seq/
for sra in SRR1944925 SRR1944926 SRR1944927 SRR1944928 SRR1944929 SRR1944930 SRR1944931 SRR1944932 SRR1944933 SRR1944934 SRR1944935
do
  prefetch -o ${sra}.sra ${sra}
  fastq-dump ${sra}.sra
  pigz ${sra}.fastq
done
```

```bash
mkdir -p ./sce/2.rawdata/ribo-seq/
cd ./sce/2.rawdata/ribo-seq/
for sra in SRR1944912 SRR1944913 SRR1944914 SRR1944915 SRR1944916 SRR1944917 SRR1944918 SRR1944919 SRR1944920 SRR1944921 SRR1944922 SRR1944923
do
  prefetch -o ${sra}.sra ${sra}
  fastq-dump ${sra}.sra
  pigz ${sra}.fastq
done
```
