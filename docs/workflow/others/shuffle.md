## `rpf_Shuffle`

### Function

Generate shuffled density as random controls.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r` | input density file |
| `-o` | output prefix |
| `-l` | gene list |
| `-s` | random seed |
| `-i` | shuffle RPFs independently for each sample |

### Example

```bash
rpf_Shuffle \
  -l ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -s 0 \
  -i \
  -o RIBO
```
