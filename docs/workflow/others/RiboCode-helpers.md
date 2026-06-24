# RiboCode/RiboTISH helpers

RiboCode/RiboTISH helpers convert external ORF-calling results into formats that can be reused in downstream annotation, visualization, and sequence extraction steps.


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


### Notes


- `ribocode_bed_format` derives ORF length from the ORF name suffix. Make sure the RiboCode BED name field follows the expected naming convention.
- The generated GenePred files are headerless; add headers manually only if a downstream tool expects them.
