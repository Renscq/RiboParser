# 4.5.12 Codon occupancy

Codon occupancy analysis assigns A/P/E-site codons and normalizes read counts at each codon against average per-codon density within ORFs.

## Command help

```bash
rpf_Occupancy -h
```

```text
usage: rpf_Occupancy [-h] -r RPF [-l LIST] -o OUTPUT [-s {E,P,A}]
                     [-f {0,1,2,all}] [-m MIN] [-n]
                     [--tis TIS] [--tts TTS]
                     [--scale {zscore,minmax}] [--stop] [--all]

Required arguments:
  -r RPF       input RPF density file
  -l LIST      gene list; default: whole
  -o OUTPUT    output prefix

Options:
  -s           E/P/A site
  -f           reading frame
  -m MIN       minimum RPF count
  -n           normalize to RPM
  --tis        discard codons after TIS
  --tts        discard codons before TTS
  --scale      zscore or minmax
  --stop       remove stop codon
  --all        output all RPF density
```

## Example

```bash
cd ./4.ribo-seq/5.riboparser/12.codon_occupancy/

for sites in E P A
do
  rpf_Occupancy \
    -l ../../../1.reference/norm/gene.norm.txt \
    -r ../05.merge/RIBO_merged.txt \
    -m 30 \
    -s "$sites" \
    -f 0 \
    --stop \
    --scale minmax \
    -o "$sites"_site \
    &>> "$sites"_site.log
done
```

## Output files

```text
A_site_codon_density.txt
A_site_codon_occupancy.txt
A_site_occupancy_corrplot.pdf
A_site_occupancy_corrplot.png
A_site_occupancy_corr.txt
A_site_occupancy_heatplot.pdf
A_site_occupancy_heatplot.png
A_site_occupancy_relative_heatplot.pdf
A_site_occupancy_relative_heatplot.png
A_site_occupancy_relative_lineplot.pdf
A_site_occupancy_relative_lineplot.png
A_site.log
```
