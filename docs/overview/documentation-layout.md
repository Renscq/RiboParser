# 1.3 Suggested documentation layout

```text
RiboParser/
├── README.md
├── mkdocs.yml
├── requirements-docs.txt
├── docs/
│   ├── overview/
│   ├── installation/
│   ├── workflow/
│   │   ├── quality-control/
│   │   ├── gene-level/
│   │   ├── codon-level/
│   │   └── smorf/
│   ├── toolkits/
│   ├── performance.md
│   ├── license.md
│   └── acknowledgements.md
└── .github/workflows/docs.yml
```

`Deployment` was removed from website navigation. Deployment notes are kept only in `UPLOAD_GUIDE.md`.
