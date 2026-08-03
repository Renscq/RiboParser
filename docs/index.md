# RiboParser documentation

**RiboParser** is a modular command-line toolkit for comprehensive RNA-seq and Ribo-seq data analysis.

This documentation keeps the repository README short and moves the complete tutorial, parameter explanations, command templates, output interpretation, and helper scripts into a structured website.

## Complete workflow

The first four steps belong to the **public pipeline** and can be handled by general bioinformatics tools (e.g. STAR, Bowtie, RSEM). The remaining steps are the **RiboParser-specific analysis**.

```text
# 1-4: Public pipeline (general tools)
Reference preparation
→ Raw data download
→ Raw data cleaning
→ Alignment and quantification

# 5-10: RiboParser-specific analysis
→ Quality control
→ Gene-level analysis
→ Codon-level analysis
→ smORF analysis
→ SeRP analysis
→ Visualization and downstream analysis
```

## Main command groups

| Group | Commands |
|---|---|
| Ribo-quality | `rpf_Reference`, `rpf_Check`, `rpf_Digest`, `rpf_Offset`, `rpf_Density`, `rpf_Merge`, `rpf_Periodicity`, `rpf_Metaplot`, `rpf_Coverage`, `rpf_Corr`, `rpf_Quant`, `rpf_Percent` |
| RNA | `rna_Offset`, `rna_Density` |
| Ribo-pausing | `rpf_Pausing`, `rpf_Occupancy`, `rpf_CoV`, `rpf_Cumulative_CoV`, `rpf_CDT`, `rpf_CST`, `rpf_Odd_Ratio`, `rpf_Meta_Codon` |
| Ribo-utils | `rpf_Shuffle`, `rpf_Shift`, `rpf_Retrieve`, `rpf_Bam2bw`, `rpf_Geneplot` |
| smORF | `smorf_scanner`, `smorf_cluster`, `smorf_evidence`, `smorf_quant` |
| SeRP | `serp_overlap`, `serp_peak`, `serp_properties` |

### Helper scripts

| Group | Commands | Purpose |
|---|---|---|
| Ribo-utils not covered above | `rpf_Shuffle`, `rpf_Shift`, `rpf_Bam2bw`, `rpf_Retrieve`, `rpf_Geneplot` | shuffled controls, translation frame shift, signal track generation, retrieve rpf density, gene plotting |
| Ribo-utils merge helpers | `merge_length`, `merge_saturation`, `merge_digestion`, `merge_offset`, `merge_offset_detail`, `merge_dst_list`, `merge_period`, `merge_metagene`, `merge_coverage`, `merge_quant`, `merge_pausing`, `merge_occupancy`, `merge_cdt`, `merge_cst`, `merge_odd_ratio` | length distribution, gene saturation, rpf digestion, rpf offset, rpf offset details, rpf density file, 3nt periodicity, rpf metagene, rpf coverage, rpf quant, translation pausing score, translation occupancy, translation codon decoding time, translation codon selection time, translation codon odd ratio |
| FASTA helpers | `fa_gc_sum`, `fa_len_flt`, `fa_len_sum`, `fa_split`, `line_feed`, `nt2aa`, `rand_seq`, `retrieve_seq`, `revs` | sequence preprocessing and extraction |
| FASTQ helpers | `fq_len_flt`, `fq_len_sum`, `fq_length`, `fq_split`, `fq_trim`, `fq2fa`, `fq2txt`, `phred_quality`, `simulate_fastq` | read preprocessing and simulation |
| bedGraph helpers | `bg2meta`, `rpm_smooth` | signal processing |
| Bowtie/RSEM helpers | `merge_bwt_log`, `merge_rsem` | mapping and quantification summary |
| RiboCode/RiboTISH helpers | `ribocode_bed_format`, `ribotish_format` | external ORF-tool output formatting |
| Unix helpers | `dos2unix` | line-ending conversion |
| oligo helpers | `get_overlap_seq`, `get_tissue_freq`, `get_win_seq` | sequence-window and oligo-related processing |
