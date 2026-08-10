# 3 New project

## Project layout

To standardize project analysis, the following directory structure is recommended.

```text
.
sce
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
│       ├── 10.quantification
│       ├── 11.shuffle
│       └── 12.retrieve
├── 4.ribo-seq
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
│       ├── 10.quantification
│       ├── 11.pausing_score
│       ├── 12.codon_occupancy
│       ├── 13.codon_decoding_time
│       ├── 14.codon_selection_time
│       ├── 15.coefficient_of_variation
│       ├── 16.cumulative_of_cov
│       ├── 17.meta_codon
│       ├── 18.odd_ratio
│       ├── 19.pause_site_plot
│       ├── 20.frame_shift
│       ├── 21.shuffle
│       ├── 22.retrieve
│       └── 23.geneplot
├── 5.smorf
│   ├── 01.scanner
│   ├── 02.cluster
│   ├── 03.evidence
│   └── 04.quant
└── 6.serp
    ├── 01.peak
    ├── 02.overlap
    ├── 03.summary
    ├── 04.peak_plot
    ├── 05.metaplot
    └── 06.properties
```

## Basic package check

```bash
riboparser -v
riboparser -c
riboparser -d
riboparser -m
```
