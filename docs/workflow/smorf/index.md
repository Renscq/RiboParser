# 4.8 smORF analysis

RiboParser provides an integrated smORF workflow.

## Workflow

```text
smorf_scanner
→ smorf_filter
→ smorf_evidence
→ smorf_integrate
```

## Main evidence levels

| Evidence | Meaning |
|---|---|
| ORF sequence | start codon, stop codon, length, strand, category |
| Kozak context | sequence or PWM-based start codon context |
| Ribo-seq signal | RPF sum, covered codons, coverage ratio |
| Periodicity | frame-specific translation evidence |
| Start/stop behavior | initiation and termination signatures |
| Multi-sample support | consistency across samples |
| Integration | ORF-level evidence table and matrix |
