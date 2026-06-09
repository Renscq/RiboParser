# 4.5 Quality control overview

This section contains complete command documentation for RNA-seq and Ribo-seq quality control.

## Modules

| Section | Command | Purpose |
|---|---|---|
| 4.5.1 | `rpf_Check` | read length distribution, filtered BAM, saturation |
| 4.5.2 | `rpf_Digest` | digestion and ligation bias |
| 4.5.3 | `rna_Offset`, `rpf_Offset` | RNA/Ribo offset table |
| 4.5.4 | `rna_Density`, `rpf_Density` | RNA/Ribo density generation |
| 4.5.5 | `merge_dst_list`, `rpf_Merge` | merge density files |
| 4.5.6 | `rpf_Periodicity` | 3-nt periodicity |
| 4.5.7 | `rpf_Metaplot` | TIS/TTS metagene profiles |
| 4.5.8 | `rpf_Coverage`, `rpf_Percent` | gene-body coverage |
| 4.5.9 | `rpf_Corr` | sample correlation |

## Directory structure

```bash
mkdir 01.qc 02.digestion 03.offset 04.density 05.merge \
  06.periodicity 07.metaplot 08.coverage 09.correlation
```
