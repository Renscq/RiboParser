# 5 Other toolkits

This section intentionally avoids duplicating commands already documented in the workflow. It summarizes only auxiliary utilities and where to find them.

## Utility groups

| Group | Commands | Purpose |
|---|---|---|
| Ribo-utils not covered above | `rpf_Shuffle`, `rpf_Bam2bw`, `rpf_Retrieve` , `rpf_Geneplot` | shuffled controls, signal track generation, retrieve rpf density, gene plotting |
| FASTA helpers | `fa_gc_sum`, `fa_len_flt`, `fa_len_sum`, `fa_split`, `line_feed`, `nt2aa`, `rand_seq`, `retrieve_seq`, `revs` | sequence preprocessing and extraction |
| FASTQ helpers | `fq_len_flt`, `fq_len_sum`, `fq_length`, `fq_split`, `fq_trim`, `fq2fa`, `fq2txt`, `phred_quality`, `simulate_fastq` | read preprocessing and simulation |
| bedGraph helpers | `bg2meta`, `rpm_smooth` | signal processing |
| Bowtie/RSEM helpers | `merge_bwt_log`, `merge_rsem` | mapping and quantification summary |
| RiboCode/RiboTISH helpers | `ribocode_bed_format`, `ribotish_format` | external ORF-tool output formatting |
| Unix helpers | `dos2unix` | line-ending conversion |
| oligo helpers | `get_overlap_seq`, `get_tissue_freq`, `get_win_seq` | sequence-window and oligo-related processing |
