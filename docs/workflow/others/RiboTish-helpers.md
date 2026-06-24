# RiboCode/RiboTISH helpers

RiboCode/RiboTISH helpers convert external ORF-calling results into formats that can be reused in downstream annotation, visualization, and sequence extraction steps.


## RiboTISH helper

### `ribotish_format`

#### Function

Filter RiboTISH ORF prediction output by RiboPvalue, amino-acid length, and ORF class, remove duplicate ORF blocks, and write three downstream files:

| Output | Description |
|---|---|
| `<prefix>.aa` | Amino-acid sequences in FASTA format. |
| `<prefix>.cds` | CDS/nucleotide sequences in FASTA format. |
| `<prefix>.genepred` | ORF coordinates in GenePred-like format. |

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input RiboTISH result table. |
| `-o` | Output prefix. |
| `-t` | ORF class to keep: `all`, `smorf`, or `annotated`. Default: `all`. |
| `-p` | Maximum `RiboPvalue`. Default: `0.05`. |
| `-l` | Minimum amino-acid length, based on `AALen`. Default: `8`. |
| `-c` | Clean fuzzy/overlapping classes by removing records whose `TisType` contains `Known` or `CDSFrameOverlap`. |

#### ORF class behavior

Before filtering, `TisType` strings are normalized by replacing `5'UTR` with `uORF` and `3'UTR` with `dORF`.

| `-t` value | Kept records |
|---|---|
| `all` | All records passing p-value and length filtering. |
| `smorf` | Records whose `TisType` contains `uORF`, `dORF`, or `Novel`. |
| `annotated` | Records whose `TisType` contains `Annotated`, `Extended`, `Internal`, or `Truncated`. |

#### Required input columns

The input table should contain the RiboTISH output columns used by the formatter, including:

```text
RiboPvalue    AALen    TisType    Gid    Tid    GenomePos    Blocks    AASeq    Seq
```

#### Example

```bash
ribotish_format \
  -i RiboTISH.results.txt \
  -p 0.05 \
  -l 8 \
  -t smorf \
  -c \
  -o RiboTISH.smorf
```

Output files:

```text
RiboTISH.smorf.aa
RiboTISH.smorf.cds
RiboTISH.smorf.genepred
```

### Notes

- `ribotish_format` removes duplicate ORFs based on `TisType` and `Blocks`, keeping the longest ORF after sorting by genomic position and amino-acid length.

- The generated GenePred files are headerless; add headers manually only if a downstream tool expects them.
