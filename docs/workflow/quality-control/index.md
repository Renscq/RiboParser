# 4.5 Quality control

The quality-control module checks aligned RNA-seq and Ribo-seq data before downstream gene-level, codon-level, smORF, and SeRP analyses.

This module covers BAM inspection, enzymatic digestion bias, offset estimation, density generation, density merging, periodicity, metagene profiles, coverage, and sample correlation.

## Recommended workflow

```text
rpf_Check
  ↓
rpf_Digest
  ↓
rna_Offset / rpf_Offset
  ↓
rna_Density / rpf_Density
  ↓
merge_dst_list → rpf_Merge
  ↓
rpf_Periodicity / rpf_Metaplot / rpf_Coverage / rpf_Percent / rpf_Corr
```

## Step summary

| Section | Main command | Purpose |
|---|---|---|
| 4.5.1 Check | `rpf_Check` | Inspect aligned BAM files, filter transcriptome reads, calculate read-length distribution, and optionally evaluate saturation. |
| 4.5.2 Enzymatic bias | `rpf_Digest` | Detect 5' and 3' end sequence preferences caused by nuclease digestion, ligation, or library-preparation bias. |
| 4.5.3 Offset table | `rna_Offset`, `rpf_Offset` | Generate RNA-seq offset tables and predict Ribo-seq P-site offset tables. |
| 4.5.4 Read density | `rna_Density`, `rpf_Density` | Convert BAM alignments to transcript-level density files. |
| 4.5.5 Merge density | `merge_dst_list`, `rpf_Merge` | Build density file lists and merge single-sample density files into multi-sample matrices. |
| 4.5.6 Periodicity | `rpf_Periodicity` | Calculate frame-specific read distribution and 3-nt periodicity. |
| 4.5.7 Metaplot | `rpf_Metaplot` | Generate TIS/TTS-centered metagene profiles. |
| 4.5.8 Coverage | `rpf_Coverage`, `rpf_Percent` | Evaluate read coverage across transcript regions and summarize coverage percentage. |
| 4.5.9 Correlation | `rpf_Corr` | Calculate sample correlation at gene level and nucleotide-position level. |

## Suggested directory structure

```bash
mkdir -p 01.qc 02.digestion 03.offset 04.density 05.merge \
         06.periodicity 07.metaplot 08.coverage 09.correlation
```

## Input dependency

| Step | Main input | Main output |
|---|---|---|
| Check | transcriptome BAM and `gene.norm.txt` | filtered BAM and read-length statistics |
| Digestion bias | filtered BAM, transcript annotation, transcript FASTA | end-sequence motif and digestion-site profiles |
| Offset table | filtered BAM or read-length range | RNA/Ribo offset tables |
| Density | filtered BAM, offset table, transcript annotation, transcript FASTA | single-sample density file |
| Merge density | single-sample density files | merged density matrix |
| Downstream QC | merged density matrix | periodicity, metaplot, coverage, and correlation summaries |

## Notes

- Run `rpf_Check` before density generation so that downstream steps use filtered and indexed BAM files.
- Use `rna_Offset` for RNA-seq and `rpf_Offset` for Ribo-seq.
- Ribo-seq codon-level analysis should only be performed after checking periodicity, metagene pattern, coverage, and replicate correlation.
- Merge helper commands are documented in the page where their outputs are most commonly used.
