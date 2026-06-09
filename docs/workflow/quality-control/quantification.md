# 4.5.10 Quantification

Ribo-seq quantification differs from RNA-seq because Ribo-seq measures ribosome occupancy within CDS regions.

## Recommended principle

To reduce artifacts from translation initiation and termination, the standard workflow excludes:

- first 15 codons downstream of the start codon
- last 5 codons upstream of the stop codon

In-frame filtering can further remove out-of-frame reads that likely represent noise.

## Command help

```bash
rpf_Quant -h
```

```text
usage: rpf_Quant [-h] -r RPF -o OUTPUT [-f {0,1,2,all}]
                 [--tis TIS] [--tts TTS] [--utr5] [--utr3]

Required arguments:
  -r RPF       input RPF density file
  -o OUTPUT    output prefix

Options:
  -f           reading frame
  --tis TIS    discard codons after TIS
  --tts TTS    discard codons before TES
  --utr5       quantify 5'UTR
  --utr3       quantify 3'UTR
```

## Ribo-seq quantification

```bash
cd ./4.ribo-seq/5.riboparser/10.quantification/

rpf_Quant \
  -r ../05.merge/RIBO_merged.txt \
  --tis 15 \
  --tts 5 \
  -o RIBO \
  &>> RIBO.log
```

## Output files

```text
RIBO_cds_rpf_quant.txt
RIBO_cds_rpm_quant.txt
RIBO_cds_rpkm_quant.txt
RIBO_cds_tpm_quant.txt
RIBO_cds_rpm_bar_plot.pdf
RIBO_cds_rpm_cdf_plot.pdf
RIBO_cds_rpm_heatmap.pdf
RIBO_cds_rpm_pca_plot.pdf
RIBO_cds_rpm_pca.txt
RIBO_total.txt
RIBO.log
```
