# 4.5.5 Merge density

```bash
merge_dst_list -h
```

```text
usage: merge_dst_list [-h] -l LIST [LIST ...] -o OUTPUT
Required arguments:
  -l LIST     density files
  -o OUTPUT   output file list
```

```bash
merge_dst_list -l ../04.density/*_rna.txt -o RNA.file.list
merge_dst_list -l ../04.density/*_rpf.txt -o RIBO.file.list
```

```bash
rpf_Merge -h
```

```text
usage: rpf_Merge [-h] -l LIST -o OUTPUT
Required arguments:
  -l LIST    sample list in TXT format
  -o OUTPUT  output prefix
```

```bash
rpf_Merge -l RNA.file.list -o RNA
rpf_Merge -l RIBO.file.list -o RIBO
```
