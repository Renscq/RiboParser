# 2.3 GitHub

Install the development version from GitHub:

```bash
git clone https://github.com/Renscq/RiboParser.git
cd RiboParser

pip install build
python -m build
pip install .
```

For editable development mode:

```bash
pip install -e .
```

## Run the test

```bash
# check the version
riboparser -v

# check the citation
riboparser -c

# check the dependency
riboparser -d

# check the module
riboparser -m
```
