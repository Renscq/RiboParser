# 2 Installation

RiboParser requires **Python >= 3.12** on a **Linux/POSIX** system.

RiboParser can be installed in any of the following three ways:

| Method | Section | Description |
|---|---|---|
| conda / micromamba | [2.1 conda / micromamba](conda.md) | Install RiboParser together with external bioinformatics dependencies in a reproducible environment |
| pip | [2.2 pip](pip.md) | Install the released Python package into an existing environment |
| GitHub | [2.3 GitHub](github.md) | Install the development version directly from the source repository |

- **conda / micromamba** is the recommended method: it manages RiboParser and the external tools required by the full workflow (e.g. Bowtie, STAR, SAMtools, RSEM) in a single reproducible environment.
- **pip** only installs the RiboParser package itself; external tools must be installed separately.
- **GitHub** provides the latest development version and is intended for developers or users who want unreleased features.
