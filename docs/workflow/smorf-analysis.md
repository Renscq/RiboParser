# 4.8 smORF analysis

RiboParser provides smORF-related tools for candidate scanning, filtering, translation evidence evaluation, and integration.

## Related commands

```bash
smorf_scanner -h
smorf_filter -h
smorf_evidence -h
smorf_integrate -h
```

## Suggested workflow

```text
smorf_scanner
→ smorf_filter
→ smorf_evidence
→ smorf_integrate
```

## Evidence types

- ORF sequence features
- start codon context
- RPF density
- frame specificity
- periodicity
- multi-sample support
- conservation evidence
- mass spectrometry support if available
