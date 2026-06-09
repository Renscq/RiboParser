# 4.7 Codon-level analysis

Codon-level analysis is the key Ribo-seq-specific downstream analysis.

## Prerequisites

Before codon-level analysis, confirm:

- reliable P-site offset
- strong 3-nt periodicity
- sufficient in-frame RPF density
- good sample reproducibility
- clear CDS enrichment

## Main analyses

| Analysis | Command |
|---|---|
| Codon pausing score | `rpf_Pausing` |
| Codon occupancy | `rpf_Occupancy` |
| Codon decoding time | `rpf_CDT` |
| Codon selection time | `rpf_CST` |
| Coefficient of variation | `rpf_CoV`, `rpf_Cumulative_CoV` |
| Meta-codon analysis | `rpf_Meta_Codon` |
| Codon odds ratio | `rpf_Odd_Ratio` |

## Typical analysis order

```text
RPF density
→ merged RPF density
→ periodicity check
→ codon pausing score
→ codon occupancy
→ codon decoding time
→ codon selection time
→ CoV
→ meta-codon analysis
```

Detailed command examples are provided in the corresponding Quality Control subsections.
