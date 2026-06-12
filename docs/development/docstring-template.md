# RiboParser docstring template

This page defines the recommended documentation style for Python modules,
classes, functions, and command-line scripts in RiboParser.

RiboParser uses **NumPy-style docstrings** for public APIs and complex internal
functions. This format is readable in source code and works well with Sphinx,
pydocstyle, Ruff, VS Code, and other Python tooling.

## File header

Command-line entry scripts should keep the lightweight RiboParser header:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-06-12
# Version: 0.2.7
# Function: Calculate RPF saturation from read-to-gene mapping.
# Input: BAM-derived read mapping table.
# Output: RPF saturation table and summary statistics.
```

## Module docstring

Major analysis modules should include a short module-level docstring after the
header:

```python
"""RPF saturation analysis.

This module provides classes and functions for estimating Ribo-seq library
saturation by random read resampling. It is mainly used to evaluate whether the
current sequencing depth is sufficient to detect translated mRNA features.
"""
```

## Function docstring

Use NumPy-style docstrings for public functions, class methods, and complex
helper functions:

```python
def function_name(arg1: str, arg2: int = 10) -> dict[str, list[int]]:
    """Short summary of the function.

    Extended description if needed. Explain why the function exists and what
    biological or computational problem it solves.

    Parameters
    ----------
    arg1 : str
        Description of the first argument.
    arg2 : int, optional
        Description of the second argument. Default is 10.

    Returns
    -------
    dict[str, list[int]]
        Description of the returned object.

    Raises
    ------
    FileNotFoundError
        Raised when the input file does not exist.
    ValueError
        Raised when an input argument is invalid.

    Notes
    -----
    Additional workflow details, biological assumptions, and coordinate
    conventions can be described here.

    Examples
    --------
    >>> result = function_name("input.txt", arg2=20)
    >>> len(result)
    10
    """
```

## Class docstring

Class docstrings should document construction parameters and important mutable
attributes:

```python
class RpfSaturation:
    """Estimate RPF saturation from read-to-gene mappings.

    Parameters
    ----------
    bam_seq_dict : dict[str, list[str]]
        Dictionary mapping read names to candidate gene IDs.
    sample_name : str
        Sample name used in output tables.
    resample_ratios : list[float]
        Resampling ratios used for saturation estimation.

    Attributes
    ----------
    mrna_dict : dict[str, list[str]]
        Dictionary mapping gene IDs to detected RPF sites.
    saturation_result : pandas.DataFrame
        Saturation result table generated after resampling.
    """
```

## RiboParser-specific sections

### Coordinate System

Use this section whenever genomic, transcriptomic, or codon-level coordinates
are involved:

```python
Coordinate System
-----------------
Input genomic intervals are treated as 0-based half-open coordinates.
Transcript-relative codon positions are reported as 0-based codon indices.
```

### Input Format

Use this section for complex tables:

```python
Input Format
------------
The input ORF table must contain the following columns:

- ``orf_id``
- ``transcript_id``
- ``chrom``
- ``strand``
- ``orf_start``
- ``orf_end``
```

### Output Format

Use this section when a function writes a structured table:

```python
Output Format
-------------
The output table contains one row per ORF and sample. The main columns are:

- ``orf_id``
- ``sample``
- ``rpf_sum``
- ``frame0_ratio``
- ``translation_evidence``
```

### Workflow

Use this section only for non-trivial workflows. Do not repeat obvious code:

```python
Workflow
--------
The scoring procedure contains four steps:

1. Calculate ORF-level RPF density.
2. Estimate frame-specific RPF enrichment.
3. Evaluate start, stop, and release signals.
4. Assign translation evidence labels.
```

## Inline comment rules

Inline comments should explain assumptions, edge cases, biological meaning, or
coordinate conventions. Avoid comments that merely repeat the code.

Good:

```python
# Use absolute values for minus-strand bedGraph because some tools store
# reverse-strand density as negative values.
density = abs(raw_density)
```

Avoid:

```python
# Loop over genes
for gene in genes:
    ...
```

## Legacy conversion rule

Legacy RiboParser tags should be migrated as follows:

| Legacy tag | NumPy-style target |
| --- | --- |
| `@Message` | First-line summary |
| `@Input` | `Parameters` or `Notes` |
| `@Return` | `Returns` |
| `@Flow` | `Workflow` |

Do not add unsupported legacy tags to new code.
