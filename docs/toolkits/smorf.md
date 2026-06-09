# 5.7 smORF commands

| Command | Function |
|---|---|
| `smorf_scanner` | Scan candidate small ORFs |
| `smorf_filter` | Filter candidate smORFs |
| `smorf_evidence` | Evaluate translation evidence using Ribo-seq P-site density |
| `smorf_integrate` | Integrate ORF-level evidence across samples |

Recommended order:

```text
smorf_scanner
→ smorf_filter
→ smorf_evidence
→ smorf_integrate
```

See `Workflow → smORF analysis` for detailed examples.
