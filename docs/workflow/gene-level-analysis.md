# 4.6 Gene-level analysis

Gene-level analysis includes quantification, coverage, correlation, and formatted extraction of gene density.

## Main commands

```bash
rpf_Quant -h
rpf_Coverage -h
rpf_Corr -h
rpf_Retrieve -h
```

## Ribo-seq gene-level quantification

```bash
cd ./4.ribo-seq/5.riboparser/10.quantification/

rpf_Quant \
  -r ../05.merge/RIBO_merged.txt \
  --tis 15 \
  --tts 5 \
  -o RIBO \
  &>> RIBO.log
```

## Retrieve and format Ribo-seq gene density

```bash
cd ./sce/4.ribo-seq/5.riboparser/18.retrieve/

rpf_Retrieve \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  -f \
  -n \
  -o RIBO \
  &>> RIBO.log
```

## Retrieve and format RNA-seq gene density

```bash
cd ./sce/3.rna-seq/5.riboparser/11.retrieve/

rpf_Retrieve \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RNA_merged.txt \
  -m 50 \
  -f \
  -n \
  -o RNA \
  &>> RNA.log
```

## Output files

```text
RIBO_retrieve.txt
RNA_retrieve.txt
RIBO_cds_rpf_quant.txt
RIBO_cds_rpm_quant.txt
RIBO_cds_rpkm_quant.txt
RIBO_cds_tpm_quant.txt
```
