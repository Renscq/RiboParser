# 4.5.5 Merge density

Density files from multiple batches and samples can be integrated for unified downstream analysis.

## Create RNA-seq sample list

```bash
cd ./3.rna-seq/5.riboparser/05.merge/

merge_dst_list \
  -l ../04.density/*_rna.txt \
  -o RNA.file.list

cat RNA.file.list
```

Example:

```text
Name                     File                                                        Type
wt_rna_YPD1              /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944924_rna.txt RNA
wt_rna_YPD2              /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944925_rna.txt RNA
wt_rna_YPD3              /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944926_rna.txt RNA
ncs2d_rna_YPD1           /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944927_rna.txt RNA
```

## Merge RNA-seq density

```bash
rpf_Merge \
  -l RNA.file.list \
  -o RNA \
  &>> RNA.log
```

## Create Ribo-seq sample list

```bash
cd ./4.ribo-seq/5.riboparser/05.merge/

merge_dst_list \
  -l ../04.density/*_rpf.txt \
  -o RIBO.file.list

cat RIBO.file.list
```

Example:

```text
Name                         File                                                        Type
wt_ribo_YPD1                 /home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944912_rpf.txt Ribo
wt_ribo_YPD2                 /home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944913_rpf.txt Ribo
wt_ribo_YPD3                 /home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944914_rpf.txt Ribo
ncs2d_ribo_YPD1              /home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944915_rpf.txt Ribo
```

## Command help

```bash
merge_dst_list -h
```

```text
usage: merge_dst_list [-h] -l LIST [LIST ...] [-o OUTPUT]

Required arguments:
  -l LIST    list of density files
  -o OUTPUT  output sample list file
```

```bash
rpf_Merge -h
```

```text
usage: rpf_Merge [-h] -l LIST -o OUTPUT

Required arguments:
  -l LIST    sample list in TXT format
  -o OUTPUT  output prefix
```

## Merge Ribo-seq density

```bash
rpf_Merge \
  -l RIBO.file.list \
  -o RIBO \
  &>> RIBO.log
```

## Output files

```text
RIBO.log
RIBO.file.list
RIBO_merged.txt
```
