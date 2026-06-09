# 3 Quick start

The complete analysis contains reference preparation, raw data processing, RNA-seq analysis, and Ribo-seq analysis.

## Recommended project structure

```text
sce/
├── 1.reference
│   ├── genome
│   ├── mrna
│   ├── ncrna
│   ├── norm
│   ├── rrna
│   ├── rsem-index
│   ├── star-index
│   └── trna
├── 2.rawdata
│   ├── ribo-seq
│   └── rna-seq
├── 3.rna-seq
│   ├── 1.cleandata
│   ├── 2.bowtie
│   ├── 3.star
│   ├── 4.quantification
│   └── 5.riboparser
│       ├── 01.qc
│       ├── 02.digestion
│       ├── 03.offset
│       ├── 04.density
│       ├── 05.merge
│       ├── 06.periodicity
│       ├── 07.metaplot
│       ├── 08.coverage
│       ├── 09.correlation
│       ├── 10.shuffle
│       └── 11.retrieve
└── 4.ribo-seq
    ├── 1.cleandata
    ├── 2.bowtie
    ├── 3.star
    ├── 4.quantification
    └── 5.riboparser
        ├── 01.qc
        ├── 02.digestion
        ├── 03.offset
        ├── 04.density
        ├── 05.merge
        ├── 06.periodicity
        ├── 07.metaplot
        ├── 08.coverage
        ├── 09.correlation
        ├── 10.quantification
        ├── 11.pausing_score
        ├── 12.codon_occupancy
        ├── 13.codon_decoding_time
        ├── 14.codon_selection_time
        ├── 15.coefficient_of_variation
        ├── 16.meta_codon
        ├── 17.shuffle
        ├── 18.retrieve
        └── 19.frame_shift
```

## Minimal workflow

```text
Prepare references
→ Download RNA-seq/Ribo-seq raw data
→ Clean reads
→ Classify reads with Bowtie
→ Align reads with STAR
→ Quantify with RSEM
→ Run RiboParser quality control
→ Generate density files
→ Merge density files
→ Run periodicity, metagene, coverage, and correlation analysis
→ Run gene-level and codon-level analyses
```
