# 2.3 GitHub

## Purpose

Install the development version directly from the GitHub source repository. Recommended for developers or users who want the latest (possibly unreleased) features.

## Requirements

- Linux/POSIX environment
- Python >= 3.12
- [git](https://git-scm.com/)

## Clone and install

```bash
git clone https://github.com/Renscq/RiboParser.git
cd RiboParser

pip install build
python -m build
pip install .
```

## Editable development mode

For development, install in editable mode so that source changes take effect immediately:

```bash
git clone https://github.com/Renscq/RiboParser.git
cd RiboParser
pip install -e .
```

## Install external dependencies

As with pip, external tools required by the full workflow must be installed separately, for example via conda (see [2.2 pip](pip.md)).
