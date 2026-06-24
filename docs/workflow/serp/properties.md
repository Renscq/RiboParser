# 4.9.3 SeRP properties

## `serp_properties`

`serp_properties` calculates codon-usage statistics and translated-protein properties from nucleotide FASTA sequences.

### Function

Use this command when you want to:

- summarize codon frequency for each sequence;
- calculate RSCU and CAI-like codon-usage tables;
- summarize whole-input codon usage;
- translate nucleotide sequences and calculate protein properties such as GRAVY, flexibility, instability index, isoelectric point, and secondary-structure fractions.

Although this page is placed under SeRP analysis, the command is sequence based. It does not read a SeRP peak table directly.

### Input

The input should be a nucleotide FASTA file containing coding sequences.

Requirements:

- each sequence length should be a multiple of 3;
- sequences are translated in the current reading frame;
- a terminal stop codon is removed before protein-property calculation;
- sequences containing internal stop codons are skipped during protein-property output.

### Parameters

| Parameter | Required | Meaning |
|---|---:|---|
| `-f` | yes | Input nucleotide FASTA file. |
| `-o` | yes | Output prefix. |

### Output

| Output | Description |
|---|---|
| `<prefix>_frequency.txt` | Per-sequence codon frequency table. |
| `<prefix>_rscu.txt` | Per-sequence RSCU table. |
| `<prefix>_cai.txt` | Per-sequence CAI-like codon usage table. |
| `<prefix>_whole_codon_usage.txt` | Whole-input codon usage summary. |
| `<prefix>.Properties.txt` | Translated-protein property table. |

The protein property table contains:

| Column | Description |
|---|---|
| `ID` | FASTA record ID. |
| `Seq` | Translated amino-acid sequence. |
| `Length` | Amino-acid sequence length after removing a terminal stop codon, if present. |
| `Gravy` | Grand average of hydropathy. |
| `Flexibility` | Average predicted flexibility. |
| `Instability` | Instability index. |
| `Isoelectric_Point` | Predicted isoelectric point. |
| `Helix` | Predicted helix fraction. |
| `Turn` | Predicted turn fraction. |
| `Sheet` | Predicted sheet fraction. |

### Examples

Calculate sequence properties for CDS sequences:

```bash
serp_properties \
  -f CDS_sequences.fa \
  -o CDS_properties
```

Calculate properties for peak-associated sequences retrieved from a FASTA file:

```bash
serp_properties \
  -f SeRP_peak_sequences.fa \
  -o SeRP_peak_properties
```

### Notes

- The command expects nucleotide sequences, not amino-acid FASTA input.
- Sequence lengths that are not divisible by 3 will raise an error during codon counting.
- Internal stop codons cause the affected sequence to be skipped in the protein-property table.
- The generated `<prefix>.Properties.txt` uses a dot before `Properties`, while the codon-usage outputs use underscores.
