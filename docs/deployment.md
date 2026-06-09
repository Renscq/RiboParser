# Deployment

This page describes how to deploy the documentation website to GitHub Pages.

## Target URL

For the repository:

```text
https://github.com/Renscq/RiboParser
```

the GitHub Pages project site URL is:

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

### 1. Commit files

```bash
git add README.md mkdocs.yml requirements-docs.txt docs .github/workflows/docs.yml
git commit -m "Build GitHub Pages documentation with MkDocs"
git push origin main
```

### 2. Enable GitHub Pages

Open:

```text
Settings → Pages → Build and deployment → Source
```

Select:

```text
GitHub Actions
```

### 3. Check GitHub Actions

Open:

```text
Actions → Deploy documentation
```

After the workflow succeeds, visit:

```text
https://renscq.github.io/RiboParser/
```

## Local preview

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

Then open:

```text
http://127.0.0.1:8000/
```
