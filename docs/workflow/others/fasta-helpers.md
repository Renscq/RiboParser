# FASTA helpers

FASTA helpers provide small command-line utilities for sequence statistics, format conversion, sequence extraction, translation, and random sequence generation. All commands are installed as console scripts after installing `RiboParser`.

### Command summary

| Command | Function | Main output |
|---|---|---|
| `fa_gc_sum` | Count A/T/C/G bases for each FASTA record. | Tab-delimited table with `name`, `A`, `T`, `C`, `G`. |
| `fa_len_sum` | Summarize sequence-length distribution. | Tab-delimited table with `length` and `count`. |
| `fa_len_flt` | Filter FASTA records by sequence length. | Filtered FASTA file. |
| `fa_split` | Split a FASTA file into multiple smaller FASTA files. | Multiple files named from the output prefix. |
| `line_feed` | Convert multi-line FASTA sequences to one sequence line per record. | Reformatted FASTA file. |
| `nt2aa` | Translate nucleotide/CDS FASTA sequences into amino-acid FASTA. | Amino-acid FASTA file. |
| `rand_seq` | Generate random DNA, RNA, or protein sequences. | FASTA or plain-text sequences printed to stdout. |
| `retrieve_seq` | Retrieve selected FASTA records by ID list. | Subset FASTA file and unmapped ID list. |
| `revs` | Print reverse, complement, and reverse-complement of an input sequence. | Text printed to stdout. |

### `fa_gc_sum`

#### Function

Calculate nucleotide composition for each FASTA record.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTA file. |
| `-o` | Output TXT/TSV file. |

#### Example

```bash
fa_gc_sum \
  -i transcript.fa \
  -o transcript_gc.txt
```

Output columns:

```text
name    A    T    C    G
```

### `fa_len_sum`

#### Function

Count how many FASTA records are present at each sequence length.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTA or FASTA.GZ file. |
| `-o` | Output TXT/TSV file. |

#### Example

```bash
fa_len_sum \
  -i reads.fa.gz \
  -o reads_length_distribution.txt
```

Output columns:

```text
length    count
```

### `fa_len_flt`

#### Function

Filter FASTA records by minimum and maximum length.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTA or FASTA.GZ file. |
| `-o` | Output filtered FASTA file. |
| `-m` | Minimum sequence length. Default: `8`. |
| `-M` | Maximum sequence length. Default: `1000`. |

#### Example

```bash
fa_len_flt \
  -i reads.fa.gz \
  -m 25 \
  -M 35 \
  -o reads_25_35.fa
```

### `fa_split`

#### Function

Split a FASTA file into several smaller FASTA files.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i`, `--input` | Input FASTA file. |
| `-n`, `--num` | Number of output files. Default: `10`. |
| `-o`, `--output` | Output file prefix. Default: `split`. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
fa_split \
  -i cds.fa \
  -n 8 \
  -o cds.part
```

This command writes files such as:

```text
cds.part1.fa
cds.part2.fa
...
cds.part8.fa
```

### `line_feed`

#### Function

Convert multi-line FASTA records into a one-line-per-sequence FASTA format. This is useful before running tools that expect each record to occupy exactly two lines: one header line and one sequence line.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTA file. |
| `-o` | Output FASTA file. |

#### Example

```bash
line_feed \
  -i genome_wrapped.fa \
  -o genome_one_line.fa
```

### `nt2aa`

#### Function

Translate nucleotide sequences into amino-acid sequences. RNA bases `U` are converted to `T` before translation. Sequences whose length is not divisible by 3, or that contain unexpected bases, are skipped with warning messages.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i`, `--input` | Input nucleotide/CDS FASTA file. |
| `-o`, `--output` | Output amino-acid FASTA file. Default: `results-aa.fa`. |
| `-t`, `--tis` | Require sequences to start with `ATG`; also reports a warning when the terminal codon is not `TAA`, `TAG`, or `TGA`. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
nt2aa \
  -i cds.fa \
  -o cds.aa.fa
```

With start-codon checking:

```bash
nt2aa \
  -i cds.fa \
  -t \
  -o cds.aa.checked.fa
```

### `rand_seq`

#### Function

Generate random DNA, RNA, or protein sequences.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-d` | Generate DNA sequences. |
| `-r` | Generate RNA sequences. |
| `-p` | Generate protein sequences. |
| `-l` | Sequence length. Default: `10`. |
| `-n` | Number of sequences. Default: `5`. |
| `-f` | Output format: `fa` or `txt`. Default: `fa`. |

Only one of `-d`, `-r`, or `-p` can be used in a single command.

#### Examples

```bash
rand_seq -d -l 30 -n 5 -f fa > random_dna.fa
rand_seq -r -l 25 -n 5 -f txt > random_rna.txt
rand_seq -p -l 12 -n 10 -f fa > random_peptide.fa
```

### `retrieve_seq`

#### Function

Retrieve FASTA records according to a text file containing sequence IDs. Each line of the ID file should contain one FASTA record ID. Leading `>` characters are allowed and will be stripped.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTA file. |
| `-n` | Text file containing sequence IDs. |
| `-u` | Output unmapped ID file. Default: `unmapped.ids`. |
| `-o` | Output FASTA file. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
retrieve_seq \
  -i transcript.fa \
  -n target_ids.txt \
  -u target_unmapped.ids \
  -o target.fa
```

### `revs`

#### Function

Print the reverse sequence, complement sequence, and reverse-complement sequence of a nucleotide string.

#### Parameter

| Parameter | Meaning |
|---|---|
| `seq` | Input nucleotide sequence as a positional argument. |

#### Example

```bash
revs ATGCCGTTA
```

Example output format:

```text
reversed :      ATTGCCGTA
complement:     TACGGCAAT
rev + comp:     TAACGGCAT
```

### Notes

- `fa_len_sum` and `fa_len_flt` can read `.gz` input when the file suffix is `.gz`.
- `retrieve_seq` uses FASTA record IDs. If a header contains spaces, only the ID parsed by Biopython is used for matching.
- `rand_seq` writes to stdout, so redirect the result with `>` when an output file is needed.
