# 2.1 pip

Install RiboParser with pip:

```bash
pip install riboparser
```

Check the installation:

```bash
riboparser -v
riboparser -c
riboparser -d
riboparser -m
```

## GitHub build + pip install

```bash
cd RiboParser

# install build dependency
pip install build

# build the package
python -m build

# install local package
pip install .
```
