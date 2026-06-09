# 4.5 Quality control

RiboParser provides RNA-seq and Ribo-seq quality-control modules.

## RNA-seq analysis directory

```bash
mkdir -p ./3.rna-seq/5.riboparser/
cd ./3.rna-seq/5.riboparser/

mkdir 01.qc 02.digestion 03.offset 04.density 05.merge \
  06.periodicity 07.metaplot 08.coverage 09.correlation 10.shuffle
```

## Ribo-seq analysis directory

```bash
mkdir -p ./4.ribo-seq/5.riboparser/
cd ./4.ribo-seq/5.riboparser/

mkdir 01.qc 02.digestion 03.offset 04.density 05.merge \
  06.periodicity 07.metaplot 08.coverage 09.correlation 10.quantification \
  11.pausing_score 12.codon_occupancy 13.codon_decoding_time 14.codon_selection_time \
  15.coefficient_of_variation 16.meta_codon 17.shuffle 18.retrieve 19.frame_shift
```

## QC modules

| Section | Module | Main command | Data type |
|---|---|---|---|
| 4.5.1 | Check | `rpf_Check` | RNA-seq / Ribo-seq |
| 4.5.2 | Enzymatic bias | `rpf_Digest` | RNA-seq / Ribo-seq |
| 4.5.3 | P-site offset | `rpf_Offset`, `rna_Offset` | Ribo-seq / RNA-seq |
| 4.5.4 | RPF density | `rpf_Density`, `rna_Density` | Ribo-seq / RNA-seq |
| 4.5.5 | Merge density | `merge_dst_list`, `rpf_Merge` | RNA-seq / Ribo-seq |
| 4.5.6 | Periodicity | `rpf_Periodicity` | Ribo-seq |
| 4.5.7 | Metaplot | `rpf_Metaplot` | RNA-seq / Ribo-seq |
| 4.5.8 | Coverage | `rpf_Coverage`, `rpf_Percent` | RNA-seq / Ribo-seq |
| 4.5.9 | Correlation | `rpf_Corr` | RNA-seq / Ribo-seq |
| 4.5.10 | Quantification | `rpf_Quant` | Ribo-seq |
| 4.5.11 | Codon pausing score | `rpf_Pausing` | Ribo-seq |
| 4.5.12 | Codon occupancy | `rpf_Occupancy` | Ribo-seq |
| 4.5.13 | Codon decoding time | `rpf_CDT` | RNA-seq + Ribo-seq |
| 4.5.14 | Codon selection time | `rpf_CST` | RNA-seq + Ribo-seq |
| 4.5.15 | Coefficient of Variation | `rpf_CoV`, `rpf_Cumulative_CoV` | Ribo-seq |
| 4.5.16 | Meta-codon analysis | `rpf_Meta_Codon` | Ribo-seq |
