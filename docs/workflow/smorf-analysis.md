# 4.8 smORF analysis

RiboParser provides a smORF workflow covering transcriptome-wide ORF scanning, rule-based filtering, Ribo-seq evidence evaluation, and evidence integration.

## 6.1 Transcriptome-wide ORF scanning

Use the entire transcriptome sequence to scan all potential ORFs from scratch, including known and novel ORFs.

```bash
smorf_scanner -h
```

```text
usage: smorf_scanner [-h] -g GENOME -a ANNOTATION [-o OUT_PREFIX]
                     [--orf-prefix ORF_PREFIX]
                     [--start-codons START_CODONS]
                     [--min-aa MIN_AA] [--max-aa MAX_AA]
                     [--scan-strand {sense,antisense,both}]
                     [--kozak-up KOZAK_UP] [--kozak-down KOZAK_DOWN]
                     [-t THREADS] [--mark-overlap] [--remove-discarded]
                     [--include-stop]

Required arguments:
  -g, --genome       input genome FASTA
  -a, --annotation   input genePred annotation

Options:
  --start-codons     comma-separated start codons, such as ATG,CTG,GTG,TTG
  --min-aa           minimum ORF length in amino acids
  --max-aa           maximum ORF length in amino acids
  --scan-strand      sense, antisense, or both
  --kozak-up         upstream nucleotides for Kozak sequence
  --kozak-down       downstream nucleotides after start codon
  -t, --threads      number of worker processes
  --mark-overlap     mark nested or overlapping ORFs
  --remove-discarded remove same-frame internal ORFs
  --include-stop     keep stop codon symbol in peptide sequence
```

Example:

```bash
smorf_scanner \
  --genome ../genome/GCF_mine_genomic.fna \
  --annotation ../norm/mine.genepred \
  --out-prefix mine \
  --start-codons ATG \
  --min-aa 8 \
  --max-aa 10000 \
  --scan-strand both \
  --kozak-up 6 \
  --kozak-down 6 \
  --mark-overlap \
  --threads 20
```

## 6.2 Filter scanned ORFs

The transcriptome-wide scan is position-based and may generate many false positives. `smorf_filter` applies rule-based filters and Kozak PWM scoring.

```bash
smorf_filter -h
```

```text
usage: smorf_filter [-h] -i INPUT [-o OUT_PREFIX]
                    [--keep-start-codons KEEP_START_CODONS]
                    [--min-aa MIN_AA] [--max-aa MAX_AA]
                    [--keep-categories KEEP_CATEGORIES]
                    [--remove-categories REMOVE_CATEGORIES]
                    [--keep-antisense] [--keep-secondary] [--keep-partial]
                    [--kozak-mode {none,annotated,builtin,pwm,sequence}]
                    [--builtin-kozak {arabidopsis,drosophila,maize,plant,rice,terrestrial_plant,vertebrate,yeast}]
                    [--kozak-pwm KOZAK_PWM] [--kozak-seq KOZAK_SEQ]
                    [--annotated-categories ANNOTATED_CATEGORIES]
                    [--min-annotated-kozak MIN_ANNOTATED_KOZAK]
                    [--fallback-builtin-kozak {arabidopsis,drosophila,maize,plant,rice,terrestrial_plant,vertebrate,yeast}]
                    [--no-kozak-fallback]
                    [--min-kozak-pwm-score MIN_KOZAK_PWM_SCORE]
                    [--export-kozak-pwm EXPORT_KOZAK_PWM]
                    [--list-builtin-kozak]
```

Example:

```bash
smorf_filter \
  -i mine.message.txt \
  -o mine.reliable \
  --min-aa 8 \
  --max-aa 10000 \
  --kozak-mode annotated \
  --keep-categories uORF,dORF,lncORF,overlap_uORF,overlap_dORF
```

## 6.3 Evaluate smORF translation evidence

Ribo-seq can validate smORFs and detect consistently translated smORFs across samples.

```bash
smorf_evidence -h
```

Important inputs:

```text
-i, --orf-table       filtered smORF table from smorf_filter
--genepred           optional genePred file for ORF exon blocks
--density-list       TSV with sample, strand, path, and optional format
--density-plus       plus-strand P-site density file
--density-minus      minus-strand P-site density file
--density            unstranded P-site density file
--density-format     auto, wig, or bedgraph
--coord-mode         0based-half-open or 1based-closed
```

Important thresholds:

```text
--min-rpf-sum
--min-covered-codon
--min-coverage-ratio
--strong-periodicity
--moderate-periodicity
--strong-start-pause
--moderate-start-pause
--strong-stop-pause
--moderate-stop-pause
--strong-release
--moderate-release
--uniform-coverage-ratio
--uniform-gini
--uniform-max-to-mean
--skewed-max-to-mean
--skewed-top-fraction
--disperse-coverage-ratio
```

Example:

```bash
smorf_evidence \
  -i mine.reliable.passed.message.txt \
  --genepred mine.genePred \
  --density-list ribo.bedgraph.list \
  -o mine.smorf.riboseq_evidence.txt
```

Example density list:

```text
sample  strand  path                                      format
ribo1   +       /project/mine/ribo/bedgraph/ribo1_plus.rpf.bedgraph   bedgraph
ribo1   -       /project/mine/ribo/bedgraph/ribo1_minus.rpf.bedgraph  bedgraph
ribo2   +       /project/mine/ribo/bedgraph/ribo2_plus.rpf.bedgraph   bedgraph
ribo2   -       /project/mine/ribo/bedgraph/ribo2_minus.rpf.bedgraph  bedgraph
```

These bedGraph files can be generated from `rpf_Bam2bw`.

## 6.4 Integrate smORF evidence

`smorf_integrate` combines long-format smORF Ribo-seq evidence into matrices and ORF-level summary tables.

```bash
smorf_integrate -h
```

```text
usage: smorf_integrate [-h] -i INPUT
                       [--output-matrix OUTPUT_MATRIX]
                       [--output-integrated OUTPUT_INTEGRATED]
                       [--capture-labels CAPTURE_LABELS]
                       [--pass-labels PASS_LABELS]
                       [--excellent-min-samples EXCELLENT_MIN_SAMPLES]

Required arguments:
  -i, --input        long-format smORF evidence table

Options:
  --output-matrix       ORF-by-sample matrix table
  --output-integrated   integrated ORF-level evidence table
  --capture-labels      evidence labels used to define captured samples
  --pass-labels         evidence labels used to define reliable samples
  --excellent-min-samples
```

Example:

```bash
smorf_integrate \
  -i mine.smorf.riboseq_evidence.txt \
  --output-matrix mine.smorf.frame_density_matrix.txt \
  --output-integrated mine.smorf.integrated_evidence.txt
```
