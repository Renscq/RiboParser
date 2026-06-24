# Merge-result helpers

Most merge helpers follow this pattern:

```text
merge_xxx -l input_files -o output_prefix
```

Examples:

```bash
merge_length -l *length_distribution.txt -o RIBO
merge_saturation -l *gene_saturation.txt -o RIBO
merge_digestion -l *pwm.txt -o RIBO
merge_offset -l *SSCBM_offset.txt -o RIBO_SSCBM
merge_dst_list -l ../04.density/*_rpf.txt -o RIBO.file.list
merge_period -l *periodicity.txt -o RIBO
merge_quant -l *quant.txt -o RIBO
merge_pausing -l *pausing_score.txt -o RIBO
merge_occupancy -l *occupancy.txt -o RIBO
merge_cdt -l *cdt.txt -o RIBO
merge_cst -l *cst.txt -o RIBO
merge_odd_ratio -l *odd_ratio.txt -o RIBO
```
