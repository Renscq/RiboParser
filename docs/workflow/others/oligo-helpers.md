## oligo helpers

Oligo helpers provide utilities for overlap detection, non-redundant small-RNA/oligo frequency summarization, and sliding-window sequence generation.

### Command summary

| Command | Function | Main output |
|---|---|---|
| `get_overlap_seq` | Detect suffix-prefix overlaps between FASTA sequences. | Overlap table. |
| `get_tissue_freq` | Merge sequence frequencies across multiple FASTA files. | Non-redundant FASTA and barcode/count summary. |
| `get_win_seq` | Generate fixed-width sliding windows from each sequence. | Window table. |

### `get_overlap_seq`

#### Function

Search for overlaps where the suffix of one sequence matches the prefix of another sequence. This is useful for checking candidate oligos, adapter-like overlaps, or potentially mergeable short sequences.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTA file. |
| `-o` | Output overlap table. |
| `-l` | Minimum seed length for overlap search. Default: `6`. The reported overlap must be longer than this value. |

#### Example

```bash
get_overlap_seq \
  -i oligo.fa \
  -l 8 \
  -o oligo.overlap.txt
```

Output columns:

```text
id1    seq1    seq1_len    id2    seq2    seq2_len    overlap
```

### `get_tissue_freq`

#### Function

Combine sequence frequencies across multiple FASTA files. Each input FASTA file is treated as one sample; the sample name is inferred from the input file basename before the first period. Identical sequences are collapsed into non-redundant records, and counts across samples are stored in the FASTA header.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Text file containing one FASTA file path per line. |
| `-b` | Output barcode/read-count summary file. Default: `all_barcodes.txt`. |
| `-o` | Output non-redundant sequence FASTA file. Default: `sRNA_tissue_freq_seq.fa`. |
| `-t` | Declared as table-format output; currently parsed but not used by the implementation. |
| `-v`, `--version` | Print script version. |

#### Input list example

```text
sample1.fa
sample2.fa
sample3.fa
```

#### Command example

```bash
get_tissue_freq \
  -i fasta.list \
  -b all_barcodes.txt \
  -o sRNA_tissue_freq_seq.fa
```

Example output FASTA header structure:

```text
>t1 count_sample1 count_sample2 count_sample3 total_count
acgtacgtacgt
```

Barcode summary output:

```text
sample1    1000
sample2    1200
sample3    900
total      3100
```

### `get_win_seq`

#### Function

Generate sliding windows from each sequence in a FASTA-like file. For each non-header line, the command writes the sequence index, the original sequence, and each fixed-width window. If a sequence is shorter than the window width, the full sequence is written as the window.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i`, `--input` | Input FASTA-like file. Header lines starting with `>` are ignored. |
| `-w`, `--win` | Window width. Default: `20`. |
| `-o`, `--output` | Output TXT file. If omitted, the output is `<input>.win`. |

#### Example

```bash
get_win_seq \
  -i oligo.fa \
  -w 20 \
  -o oligo.win20.txt
```

Output columns are written without a header:

```text
sequence_index    full_sequence    window_sequence
```

### Notes

- `get_overlap_seq` performs pairwise comparisons among all FASTA records. Runtime increases quickly as the number of input sequences grows.
- `get_tissue_freq` treats sequence strings as the unit of counting. Sequence IDs are not used for collapsing.
- `get_win_seq` ignores FASTA headers and processes each sequence line independently. Use `line_feed` first if your FASTA sequences are wrapped across multiple lines.
