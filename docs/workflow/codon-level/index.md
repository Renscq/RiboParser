# 4.7 Codon-level analysis

Codon-level analysis is the most Ribo-seq-specific part of RiboParser.

## Required QC before codon-level analysis

- reliable offset table
- strong 3-nt periodicity
- sufficient RPF count
- reproducible biological replicates
- clear CDS enrichment
- reasonable digestion bias

## Submodules

| Section | Command | Purpose |
|---|---|---|
| 4.7.1 | `rpf_Pausing` | codon pausing score |
| 4.7.2 | `rpf_Occupancy` | codon occupancy |
| 4.7.3 | `rpf_CDT` | codon decoding time |
| 4.7.4 | `rpf_CST` | codon selection time |
| 4.7.5 | `rpf_CoV`, `rpf_Cumulative_CoV` | coverage-dependent variation |
| 4.7.6 | `rpf_Meta_Codon` | meta-codon profiles |
| 4.7.7 | `rpf_Odd_Ratio` | codon odds-ratio enrichment |
