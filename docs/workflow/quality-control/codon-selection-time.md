# 4.5.14 Codon selection time

Codon usage is unequal across genomes. Codon selection time can be used to quantify codon-level translational selection patterns.

## Command help

```bash
rpf_CST -h
```

```text
usage: rpf_CST [-h] --rpf RPF --rna RNA [-l LIST] -o OUTPUT
               [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN]
               [-t TIMES] [--tis TIS] [--tts TTS]
               [--scale {zscore,minmax}] [--stop]

Required arguments:
  --rpf RPF     input Ribo-seq density file
  --rna RNA     input RNA-seq density file
  -l LIST       gene list; default: whole
  -o OUTPUT     output prefix

Options:
  -s            E/P/A site
  -f            reading frame
  -m MIN        minimum RPF count
  -t TIMES      iteration number
  --tis         discard codons after TIS
  --tts         discard codons before TTS
  --scale       zscore or minmax
  --stop        remove stop codon
```

## Example

```bash
cd ./4.ribo-seq/5.riboparser/14.codon_selection_time/

for sites in E P A
do
  rpf_CST \
    -l ../../../1.reference/norm/gene.norm.txt \
    --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
    --rpf ../05.merge/RIBO_merged.txt \
    --stop \
    -m 50 \
    -f 0 \
    -s $sites \
    --tis 10 \
    --tts 5 \
    -o "$sites"_site \
    &>> "$sites"_site.log
done
```

## Output files

```text
A_site_codon_selection_time.txt
A_site_iterative_codon_selection_time.txt
A_site_cst_corrplot.pdf
A_site_cst_corrplot.png
A_site_cst_corr.txt
A_site_cst_heatplot.pdf
A_site_cst_heatplot.png
A_site.log
```
