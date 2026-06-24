# FASTQ helpers

FASTQ helpers provide utilities for format conversion, read-length filtering and summarization, quality-score checking, trimming, splitting, and simulation.

### Command summary

| Command | Function | Main output |
|---|---|---|
| `fq2fa` | Convert FASTQ to FASTA. | FASTA file. |
| `fq2txt` | Extract read sequences from FASTQ. | Plain-text sequence file. |
| `fq_len_flt` | Filter FASTQ reads by length. | Kept FASTQ and dropped FASTQ. |
| `fq_len_sum` | Summarize read-length distribution for one FASTQ file. | Length-count table. |
| `fq_length` | Summarize read lengths for one or more FASTQ/FASTQ.GZ files. | Multi-sample length-count table. |
| `fq_split` | Split an interleaved FASTQ file into paired FASTQ files. | R1/R2 FASTQ files. |
| `fq_trim` | Trim FASTQ reads to a fixed length. | Trimmed FASTQ file. |
| `phred_quality` | Infer quality-score encoding from FASTQ quality characters. | Quality format report printed to stdout. |
| `simulate_fastq` | Generate simulated FASTQ reads. | Simulated FASTQ file. |

### `fq2fa`

#### Function

Convert FASTQ records to FASTA format. By default, the original FASTQ read names are preserved. With `-n`, records are renamed sequentially as `t1`, `t2`, and so on.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTQ file. |
| `-o` | Output FASTA file. |
| `-n` | Assign new sequential read IDs. |
| `-v`, `--version` | Print script version. |

#### Examples

```bash
fq2fa \
  -i sample.fq \
  -o sample.fa
```

```bash
fq2fa \
  -i sample.fq \
  -n \
  -o sample.renamed.fa
```

### `fq2txt`

#### Function

Extract only the sequence line from each FASTQ record and write one sequence per line.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTQ file. |
| `-o` | Output TXT file. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
fq2txt \
  -i sample.fq \
  -o sample.seq.txt
```

### `fq_len_flt`

#### Function

Filter FASTQ reads by sequence length. Reads within the specified length range are written to the output FASTQ file; reads outside the range are written to an additional dropped-read file named `<output>.drop`.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTQ file. |
| `-o` | Output FASTQ file for retained reads. |
| `-m` | Minimum read length. |
| `-M` | Maximum read length. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
fq_len_flt \
  -i sample.fq \
  -m 25 \
  -M 35 \
  -o sample_25_35.fq
```

Output files:

```text
sample_25_35.fq
sample_25_35.fq.drop
```

### `fq_len_sum`

#### Function

Summarize read-length distribution for a single FASTQ file.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTQ file. |
| `-o` | Output TXT/TSV file. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
fq_len_sum \
  -i sample.fq \
  -o sample_length_distribution.txt
```

Output format:

```text
length    reads_number
25        1200
26        1530
...
Total reads number    2730
```

### `fq_length`

#### Function

Count read lengths for one or more FASTQ or FASTQ.GZ files and merge the results into one table.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i`, `--input` | One or more FASTQ/FASTQ.GZ files. |
| `-o`, `--output` | Output file name. Default: `fq_length_distr.txt`. |

#### Example

```bash
fq_length \
  -i sample1.fq.gz sample2.fq.gz sample3.fq.gz \
  -o all_samples_length_distribution.txt
```

Output columns:

```text
Sample    Length    Read_Count
```

### `fq_split`

#### Function

Intended to split an interleaved FASTQ file into paired FASTQ files.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTQ or FASTQ.GZ file. |
| `-a` | Adapter table file. The current implementation accepts this parameter but does not use it internally. |
| `-o` | Output prefix. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
fq_split \
  -i interleaved.fq.gz \
  -a adapter.txt \
  -o sample
```

Expected output naming:

```text
sample.R1.fq
sample.R2.fq
```

### `fq_trim`

#### Function

Trim each FASTQ read to a fixed nucleotide length. Reads shorter than the target length are skipped.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input FASTQ file. |
| `-o` | Output trimmed FASTQ file. |
| `-l` | Target trim length. Default: `30`. |
| `-d` | Declared as “discard reads shorter than minimum length”; in the current implementation, short reads are skipped regardless of this flag. |

#### Example

```bash
fq_trim \
  -i sample.fq \
  -l 30 \
  -o sample.trim30.fq
```

### `phred_quality`

#### Function

Infer the likely FASTQ quality-score encoding by inspecting quality-score characters from the first `N` reads.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i`, `--input` | Input FASTQ or FASTQ.GZ file. |
| `-n`, `--number` | Number of sequences to inspect. Default: `40000`. |

#### Example

```bash
phred_quality \
  -i sample.fq.gz \
  -n 40000
```

Example output categories include:

```text
Sanger | Phred+33
Solexa | Solexa+64
Illumina1.3 | Phred+64
Illumina1.5 | Phred+64
Illumina1.8 | Phred+33
```

### `simulate_fastq`

#### Function

Generate a simulated FASTQ file with random DNA sequences and random Phred+33 quality strings.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-q`, `--qmin` | Minimum quality score. |
| `-Q`, `--qmax` | Maximum quality score. |
| `-m`, `--minlen` | Minimum read length. |
| `-M`, `--maxlen` | Maximum read length. |
| `-n`, `--number` | Number of reads to generate. |
| `-o`, `--output` | Output FASTQ file. Default: `simulated.fq`. |

#### Example

```bash
simulate_fastq \
  -q 30 \
  -Q 40 \
  -m 25 \
  -M 35 \
  -n 100000 \
  -o simulated_25_35.fq
```

### Notes

- `fq_length` supports multiple input files and `.gz` input; `fq_len_sum` is designed for a single plain FASTQ file.
- `fq_len_flt` writes both retained and dropped reads, which is useful for checking whether read-length filtering is too strict.
- For production pipelines, prefer testing each helper with a small FASTQ subset before running on full sequencing files.
