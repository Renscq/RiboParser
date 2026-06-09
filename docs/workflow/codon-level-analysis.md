# 4.7 Codon-level analysis

Codon-level analysis is the most Ribo-seq-specific downstream component of RiboParser.

## Main analyses

- codon pausing score
- codon occupancy
- codon decoding time
- codon selection time
- codon odds ratio
- coefficient of variation
- meta-codon analysis

## Related commands

```bash
rpf_Pausing -h
rpf_Occupancy -h
rpf_CDT -h
rpf_CST -h
rpf_Odd_Ratio -h
rpf_CoV -h
rpf_Meta_Codon -h
```

## Prerequisites

Before codon-level analysis, confirm:

- reliable P-site offset
- strong 3-nt periodicity
- sufficient in-frame RPF density
- good sample reproducibility
- clear CDS enrichment
