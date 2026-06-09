# 4.5.13 Codon decoding time

Codon decoding time uses Ribo-seq and RNA-seq jointly. RNA-seq can be used as a correction for mRNA abundance or transcript state.

## Command help

```bash
rpf_CDT -h
```

```text
usage: rpf_CDT [-h] --rpf RPF --rna RNA -l LIST -o OUTPUT
               [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN]
               [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop]

Required arguments:
  --rpf RPF     input Ribo-seq density file
  --rna RNA     input RNA-seq density file
  -l LIST       gene list
  -o OUTPUT     output prefix

Options:
  -s            E/P/A site
  -f            reading frame
  -m MIN        minimum RPF count
  --tis         discard codons after TIS
  --tts         discard codons before TTS
  --scale       zscore or minmax
  --stop        remove stop codon
```

## Example

```bash
cd ./4.ribo-seq/5.riboparser/13.codon_decoding_time/

for sites in E P A
do
  rpf_CDT \
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
A_site_cdt.txt
A_site_cdt_corrplot.pdf
A_site_cdt_corrplot.png
A_site_cdt_corr.txt
A_site_cdt_heatplot.pdf
A_site_cdt_heatplot.png
A_site.log
```
