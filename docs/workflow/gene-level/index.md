# 4.6 Gene-level analysis

## Function

Gene-level analysis summarizes Ribo-seq or RNA-seq density at transcript and gene scale. These modules quantify CDS-level abundance, detect in-gene reading-frame switches, and draw IGV-like density profiles for genes of interest.

## Submodules

| Section | Command | Description |
|---|---|---|
| 4.6.1 Gene quantification | `rpf_Quant` | Summarize RPF counts, RPM, RPKM, and TPM for CDS and optional UTR regions. |
| 4.6.2 Frameshift detection | `rpf_Shift` | Detect in-gene reading-frame switches (programmed ribosomal frameshifts) from the merged density file. |
| 4.6.3 Gene density plot | `rpf_Geneplot` | Draw IGV-like transcript/genome-coordinate density profiles for genes of interest. |

## Suggested directory structure

The gene-level modules are organized under `sce/4.ribo-seq/5.riboparser/` (see the project layout in [3 New project](../../new-project.md)). Sub-directories `01.qc` to `09.correlation` are created by the quality-control modules (4.5); create the gene-level ones (`10.quantification`, `20.frame_shift`, and `23.geneplot`, matching sections 4.6.1-4.6.3) as follows:

```bash
mkdir -p 10.quantification 20.frame_shift 23.geneplot
```

## Input dependency

The gene-level modules consume the merged density matrix produced by `rpf_Merge` in `05.merge` (see [4.5.5 Merge density](../quality-control/merge-density.md)), usually `<prefix>_rpf_merged.jsonl.gz` (or the plain TXT table). The per-module input and output are summarized below:

| Section | Main input | Main output |
|---|---|---|
| Gene quantification | merged density (`-r`) | raw RPF counts and RPM/RPKM/TPM tables plus QC plots |
| Frameshift detection | merged frame-resolved density (`-r`, use the JSONL file) | frameshift candidates, classification summary, and per-candidate plots |
| Gene density plot | merged density (`-r`) and a target gene or target list (`-g` / `--target-list`) | IGV-like per-target, per-sample density figures |

`rpf_Shift` needs the frame-resolved density records, so it should be run on the `*_rpf_merged.jsonl.gz` output of `rpf_Merge` rather than a frame-unresolved summary. `rpf_Geneplot` reads the same merged file and accepts both the compact JSONL and the plain TXT table.

## Notes

- For standard translation-level quantification, use CDS density and trim codons near the start and stop codons to reduce initiation and termination artifacts.
- Gene-level metrics are useful for expression overview, replicate quality control, and downstream differential analysis, but they should be interpreted together with read periodicity, metagene profiles, and coverage uniformity.
- Frameshift candidates from `rpf_Shift` are statistical detections: always validate the shift site and direction manually (e.g. with the per-candidate gene plots) before reporting biological conclusions. `rpf_Geneplot` is the general-purpose tool for visually inspecting the density of any gene of interest.
