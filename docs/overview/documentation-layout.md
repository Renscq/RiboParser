# 1.3 Suggested documentation layout

The long README is reorganized into a website-style documentation system.

## Repository-level files

```text
RiboParser/
├── README.md
├── mkdocs.yml
├── requirements-docs.txt
├── docs/
└── .github/
    └── workflows/
        └── docs.yml
```

## Documentation layout

```text
docs/
├── index.md
├── overview/
├── installation/
├── workflow/
│   ├── reference-preparation.md
│   ├── raw-data-download.md
│   ├── raw-data-cleaning.md
│   ├── alignment-and-quantification.md
│   ├── quality-control/
│   ├── gene-level-analysis.md
│   ├── codon-level-analysis.md
│   └── smorf-analysis.md
├── toolkits/
├── performance.md
├── license.md
├── acknowledgements.md
└── deployment.md
```

## Principle

- `README.md`: GitHub repository homepage.
- `docs/`: complete documentation website source.
- `mkdocs.yml`: website navigation and theme configuration.
- `.github/workflows/docs.yml`: GitHub Pages deployment workflow.
