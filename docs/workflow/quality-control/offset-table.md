# 4.5.3 Offset table

## Purpose

Offset tables define how read positions are shifted to assign RNA-seq/Ribo-seq density to transcript coordinates. Ribo-seq codon-level analysis requires accurate P-site offset prediction. RNA-seq usually uses a constant offset.

## `rna_Offset`

### Function

Generate a simple RNA-seq offset table by assigning a constant offset to all retained read lengths.

### Input files

| Input | Description |
|---|---|
| none | this command only needs read length range and offset value |
| output prefix | sample name |

### Parameters

| Parameter | Meaning |
|---|---|
| `-o / --output` | output prefix |
| `-m / --min` | minimum read length |
| `-M / --max` | maximum read length |
| `-e / --offset` | constant offset assigned to each read length |

### Example

```bash
for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rna_Offset \
    -m 27 \
    -M 50 \
    -e 12 \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

### Output

| Output | Description |
|---|---|
| `sample_offset.txt` | RNA-seq offset table |
| `sample.log` | running log |

## `rpf_Offset`

### Function

Predict Ribo-seq P-site offsets using start/stop codon signal and/or ribosome structure-based models.

### Input files

| Input | Description |
|---|---|
| BAM | filtered transcriptome BAM |
| transcript annotation | `gene.norm.txt` |
| expected RPF length | usually the dominant RPF length |
| output prefix | sample name |

### Parameters

| Parameter | Meaning |
|---|---|
| `-b / --bam` | input BAM file |
| `-t / --transcript` | normalized transcript annotation |
| `-o / --output` | output prefix |
| `--mode` | offset model, such as `SSCBM` or `RSBM` |
| `-a` | align reads to `both`, `tis`, or `tts` |
| `-m / --min` | minimum read length |
| `-M / --max` | maximum read length |
| `-p / --exp-peak` | expected RPF length fitted to ribosome structure |
| `-s / --shift` | manually set P-site shift for different read lengths |
| `-l` | only use longest transcript per gene |
| `-d` | output detailed TIS/TTS end profiles |
| `--silence` | suppress verbose output |

### Example

```bash
for bam in ../01.qc/*.bam
do
  prefix_name=$(basename $bam .bam)

  rpf_Offset \
    -b $bam \
    -m 27 \
    -M 33 \
    -p 30 \
    -d \
    -t ../../../1.reference/norm/gene.norm.txt \
    -o $prefix_name \
    &>> $prefix_name".log"
done
```

### Outputs

| Output | Description |
|---|---|
| `sample_SSCBM_offset.txt` | start/stop codon-based offset table |
| `sample_SSCBM_offset.pdf/png` | SSCBM offset plot |
| `sample_SSCBM_offset_scale.pdf/png` | scaled SSCBM plot |
| `sample_RSBM_offset.txt` | ribosome structure-based offset table |
| `sample_RSBM_offset.pdf/png` | RSBM offset plot |
| `sample_tis_5end.txt` | TIS-aligned 5' end distribution |
| `sample_tis_3end.txt` | TIS-aligned 3' end distribution |
| `sample_tts_5end.txt` | TTS-aligned 5' end distribution |
| `sample_tts_3end.txt` | TTS-aligned 3' end distribution |

## Merge offset results

### `merge_offset_detail`

| Parameter | Meaning |
|---|---|
| `-l` | detailed end-distribution files, for example `*end.txt` |
| `-o` | output prefix |

```bash
merge_offset_detail -l *end.txt -o RIBO
```

### `merge_offset`

| Parameter | Meaning |
|---|---|
| `-l` | offset summary files, for example `*SSCBM_offset.txt` |
| `-o` | output prefix |

```bash
merge_offset -l *SSCBM_offset.txt -o RIBO_SSCBM
merge_offset -l *RSBM_offset.txt -o RIBO_RSBM
```

## Interpretation

A reliable P-site offset should show a stable dominant offset across biologically meaningful read lengths. Unstable or flat TIS/TTS profiles usually indicate weak periodicity, low read depth, or poor library quality.
