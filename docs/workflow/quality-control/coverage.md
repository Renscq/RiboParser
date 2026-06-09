# 4.5.8 Coverage

Coverage analysis evaluates the distribution of read density along gene bodies.

## Interpretation

RiboParser provides three complementary views:

- `Aggregate density profile`: smoothed mean density across all genes
- `Gene-specific heatmap`: gene-by-bin matrix of normalized density
- `Binned coverage distribution`: proportion of genes passing coverage thresholds in each bin

## Command help

```bash
rpf_Coverage -h
```

```text
usage: rpf_Coverage [-h] -t TRANSCRIPT -r RPF [-o OUTPUT]
                    [-f {0,1,2,all}] [-m MIN] [-b BIN] [-n]
                    [--thread THREAD] [--outlier] [--set {intersect,union}]
                    [--heat] [--bar]

Required arguments:
  -t TRANSCRIPT  transcript annotation in TXT format
  -r RPF         RPF density file

Options:
  -f             reading frame
  -m MIN         minimum RPF count
  -b BIN         bins for 5'UTR, CDS, 3'UTR
  -n             normalize to RPM
  --outlier      filter outliers
  --heat         draw heatmap
  --bar          draw barplot
```

## RNA-seq coverage

```bash
cd ./3.rna-seq/5.riboparser/08.coverage/

rpf_Coverage \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RNA_merged.txt \
  -m 50 \
  --outlier \
  -b 10,100,10 \
  -n \
  --heat \
  -o RNA \
  &>> RNA.log

rpf_Percent \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RNA_merged.txt \
  -n \
  -m 50 \
  -f 0 \
  -o RNA \
  &>> RNA.log
```

## Ribo-seq coverage

```bash
cd ./4.ribo-seq/5.riboparser/08.coverage/

rpf_Coverage \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -m 50 \
  --outlier \
  -b 10,150,10 \
  -n \
  --heat \
  -o RIBO \
  &>> RIBO.log
```

## Output files

```text
RIBO_SRR1944912_10_150_10_coverage.txt
RIBO_SRR1944912_10_150_10_heat_plot.png
RIBO_SRR1944912_coverage_bar_plot.pdf
RIBO_SRR1944912_coverage_bar_plot.png
RIBO_SRR1944912_coverage_line_plot.pdf
RIBO_SRR1944912_coverage_line_plot.png
RIBO.log
```
