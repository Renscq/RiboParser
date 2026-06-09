# 5.8 Helper scripts

## Bowtie

| Command | Function |
|---|---|
| `merge_bwt_log` | Merge Bowtie mapping log files |

## RSEM

| Command | Function |
|---|---|
| `merge_rsem` | Merge RSEM `expected_count`, `TPM`, or `FPKM` columns |

## merge_ribo

| Command | Function |
|---|---|
| `merge_length` | Merge read length distribution |
| `merge_saturation` | Merge saturation results |
| `merge_digestion` | Merge digestion bias results |
| `merge_offset` | Merge offset results |
| `merge_offset_detail` | Merge detailed offset end-distribution files |
| `merge_dst_list` | Create density file list |
| `merge_period` | Merge periodicity results |
| `merge_metagene` | Merge metagene results |
| `merge_coverage` | Merge coverage results |
| `merge_quant` | Merge quantification results |
| `merge_pausing` | Merge pausing results |
| `merge_occupancy` | Merge occupancy results |
| `merge_cdt` | Merge codon decoding time results |
| `merge_cst` | Merge codon selection time results |
| `merge_odd_ratio` | Merge odds ratio results |

## FASTA utilities

| Command | Function |
|---|---|
| `retrieve_seq` | Retrieve FASTA sequence by ID |
| `fa_gc_sum` | GC summary |
| `fa_len_flt` | Filter FASTA by length |
| `fa_len_sum` | Summarize FASTA sequence length |
| `fa_split` | Split FASTA file |
| `nt2aa` | Translate nucleotide sequence to amino acid sequence |
| `rand_seq` | Generate random sequences |
| `revs` | Reverse or reverse-complement sequences |

## FASTQ utilities

| Command | Function |
|---|---|
| `fq_len_flt` | Filter FASTQ by length |
| `fq_len_sum` | Summarize FASTQ read length |
| `fq_length` | Calculate FASTQ read length distribution |
| `fq_split` | Split FASTQ file |
| `fq_trim` | Trim FASTQ reads |
| `fq2fa` | Convert FASTQ to FASTA |
| `fq2txt` | Convert FASTQ to text |
| `phred_quality` | Evaluate Phred quality |
| `simulate_fastq` | Simulate FASTQ reads |

## Other utilities

| Command | Function |
|---|---|
| `bg2meta` | bedGraph meta-profile processing |
| `rpm_smooth` | Smooth RPM signal |
| `ribocode_bed_format` | Format RiboCode BED output |
| `ribotish_format` | Format RiboTISH output |
| `dos2unix` | Convert Windows line endings |
