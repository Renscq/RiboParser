## RiboCode/RiboTISH helpers

RiboCode/RiboTISH helpers convert external ORF-calling results into formats that can be reused in downstream annotation, visualization, and sequence extraction steps.

### Command summary

| Command | Function | Main output |
|---|---|---|
| `ribocode_bed_format` | Convert RiboCode BED output to GenePred format. | `<prefix>.genepred`. |
| `ribotish_format` | Filter and convert RiboTISH output to amino-acid FASTA, CDS FASTA, and GenePred format. | `<prefix>.aa`, `<prefix>.cds`, `<prefix>.genepred`. |

## RiboCode helper

### `ribocode_bed_format`

#### Function

Convert RiboCode BED-format ORF output into GenePred format. The command can filter ORFs by class and minimum amino-acid length, then adjusts the terminal stop-codon boundary before writing the GenePred table.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input RiboCode BED file. |
| `-o` | Output prefix. The final file is `<prefix>.genepred`. |
| `-t` | ORF class to keep: `all`, `smorf`, or `annotated`. Default: `all`. |
| `-l` | Minimum amino-acid length. Default: `8`. |
| `-c` | Declared as “drop fuzzy ORFs”; currently reserved and not applied by the implementation. |

#### ORF class behavior

| `-t` value | Kept ORF types |
|---|---|
| `all` | All ORF types that pass length filtering. |
| `smorf` | `uORF`, `dORF`, `novel`, `dORF,novel`, `uORF,novel`. |
| `annotated` | `internal`, `annotated`. |

#### Example

```bash
ribocode_bed_format \
  -i RiboCode_ORFs.bed \
  -t smorf \
  -l 8 \
  -o RiboCode.smorf
```

Output file:

```text
RiboCode.smorf.genepred
```

### Output GenePred columns

The output is a headerless GenePred-like table:

```text
name    chrom    strand    txStart    txEnd    cdsStart    cdsEnd    exonCount    exonStarts    exonEnds
```

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
- `ribocode_bed_format` derives ORF length from the ORF name suffix. Make sure the RiboCode BED name field follows the expected naming convention.
- The generated GenePred files are headerless; add headers manually only if a downstream tool expects them.
