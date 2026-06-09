# 4.5.9 Correlation

Correlation analysis evaluates reproducibility across samples.

## Ribo-seq interpretation

RiboParser provides a hierarchical framework:

- `Gene-level reproducibility`: total RPFs across gene bodies; suitable for global translation consistency
- `ORF-level reproducibility`: nucleotide-resolution RPF quantification within ORFs; useful for localized translational variation

## Command help

```bash
rpf_Corr -h
```

```text
usage: rpf_Corr [-h] -r RPF -o OUTPUT

Required arguments:
  -r RPF      input RPF density file
  -o OUTPUT   output prefix
```

## RNA-seq correlation

```bash
cd ./3.rna-seq/5.riboparser/09.correlation/

rpf_Corr \
  -r ../05.merge/RNA_merged.txt \
  -o RNA \
  &>> RNA.log
```

## Ribo-seq correlation

```bash
cd ./4.ribo-seq/5.riboparser/09.correlation/

rpf_Corr \
  -r ../05.merge/RIBO_merged.txt \
  -o RIBO \
  &>> RIBO.log
```

## Output files

```text
RIBO_gene_corr_f0.txt
RIBO_gene_corr_f1.txt
RIBO_gene_corr_f2.txt
RIBO_gene_corr_frame.txt
RIBO_gene_correlation_plot.pdf
RIBO_gene_correlation_plot.png
RIBO_rpf_corr_f0.txt
RIBO_rpf_corr_f1.txt
RIBO_rpf_corr_f2.txt
RIBO_rpf_corr_frame.txt
RIBO_rpf_correlation_plot.pdf
RIBO_rpf_correlation_plot.png
RIBO.log
```
