# 3 Quick start

Run the basic checks:

```bash
riboparser -v
riboparser -c
riboparser -d
riboparser -m
```

## Minimal workflow skeleton

```bash
# 1. Prepare reference
rpf_Reference \
  -g genome.fa \
  -t annotation.gtf \
  -u 30 \
  -o gene

# 2. Check Ribo-seq BAM
rpf_Check \
  -b sample.bam \
  -t gene.norm.txt \
  -o sample

# 3. Generate RPF density
rpf_Density \
  -b sample.bam \
  -m 27 \
  -M 33 \
  --period 40 \
  -p sample_SSCBM_offset.txt \
  -s gene.norm.rna.fa \
  -t gene.norm.txt \
  -o sample

# 4. Merge density files
merge_dst_list -l *_rpf.txt -o RIBO.file.list
rpf_Merge -l RIBO.file.list -o RIBO
```

## Analysis order

```text
Reference preparation
→ raw data download
→ raw data cleaning
→ alignment and quantification
→ quality control
→ density generation
→ merging
→ periodicity / metaplot / coverage / correlation
→ gene-level analysis
→ codon-level analysis
→ smORF / SeRP analysis
```
