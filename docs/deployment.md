# Deployment

This page describes how to deploy this documentation to GitHub Pages.

## Target URL

```text
https://renscq.github.io/RiboParser/
```

## Required files

```text
README.md
mkdocs.yml
requirements-docs.txt
docs/
.github/workflows/docs.yml
```

## Deployment steps

```bash
git add README.md mkdocs.yml requirements-docs.txt docs .github/workflows/docs.yml
git commit -m "Enrich RiboParser documentation from original README"
git push origin main
```

Then open GitHub:

```text
Settings → Pages → Build and deployment → Source → GitHub Actions
```

Check:

```text
Actions → Deploy documentation
```
