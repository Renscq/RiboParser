# 1.3 Suggested documentation layout

The original README was a complete tutorial. In this website version, the tutorial is split into focused pages.

```text
RiboParser/
├── README.md
├── mkdocs.yml
├── requirements-docs.txt
├── docs/
│   ├── index.md
│   ├── overview/
│   ├── installation/
│   ├── workflow/
│   │   ├── reference-preparation.md
│   │   ├── raw-data-download.md
│   │   ├── raw-data-cleaning.md
│   │   ├── alignment-and-quantification.md
│   │   ├── quality-control/
│   │   ├── gene-level-analysis.md
│   │   ├── codon-level-analysis.md
│   │   └── smorf-analysis.md
│   ├── toolkits/
│   ├── performance.md
│   ├── license.md
│   ├── acknowledgements.md
│   └── deployment.md
└── .github/
    └── workflows/
        └── docs.yml
```

## Principle

- `README.md`: short GitHub homepage
- `docs/`: full tutorial and command documentation
- `mkdocs.yml`: navigation and theme configuration
- `.github/workflows/docs.yml`: GitHub Pages deployment
