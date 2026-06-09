# 4.5.11 Codon pausing score

Ribosome profiling can capture ribosome pausing. If a codon is translated more slowly, increased ribosome occupancy is expected around that codon.

Pause scores are calculated by normalizing codon-level read counts against each gene's mean read density.

## Command help

```bash
rpf_Pausing -h
```

```text
usage: rpf_Pausing [-h] -r RPF [-l LIST] -o OUTPUT [-s {E,P,A}]
                   [-f {0,1,2,all}] [-b BACKGROUND] [-m MIN]
                   [--tis TIS] [--tts TTS] [-n]
                   [--scale {zscore,minmax}] [--stop]
                   [--fig {none,png,pdf}] [--all]

Required arguments:
  -r RPF       input RPF density file
  -l LIST      gene list; default: whole
  -o OUTPUT    output prefix

Options:
  -s           E/P/A site
  -f           reading frame
  -b           background codon number
  -m MIN       minimum RPF count
  --tis        discard codons after TIS
  --tts        discard codons before TTS
  --scale      zscore or minmax
  --stop       remove stop codon
  --all        output all gene-level pausing scores
```

## Example

```bash
cd ./4.ribo-seq/5.riboparser/11.pausing_score/

for sites in E P A
do
  rpf_Pausing \
    -l ../../../1.reference/norm/gene.norm.txt \
    -r ../05.merge/RIBO_merged.txt \
    -b 0 \
    --stop \
    -m 30 \
    -s $sites \
    -f 0 \
    --scale minmax \
    -o "$sites"_site \
    &>> "$sites"_site.log
done
```

## Output files

```text
A_site_cds_codon_pausing_score.txt
A_site_cds_pausing_score.txt
A_site_sum_codon_pausing_score.txt
A_site_total_pausing_heatplot.pdf
A_site_total_pausing_heatplot.png
A_site_valid_pausing_heatplot.pdf
A_site_valid_pausing_heatplot.png
A_site.log
```
