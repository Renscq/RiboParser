## `rpf_Shift`

### Function

Generate shuffled density as random controls.

### Parameters

| Parameter | Meaning |
|---|---|
| `-r` | input density file |
| `-o` | the prefix of output file |
| `-t` | the name of input transcript file in TXT format |
| `-p` | the minimum in-frame value for frame shifting screen, range [0 - 1]. (default: 45) |
| `-m` | retain transcript with more than minimum RPFs. (default: 50) |
| `--tis` | the number of codons after TIS will be discarded.. (default: 0 AA) |
| `--tts` | the number of codons after TTS will be discarded.. (default: 0 AA) |

### Example

```bash
rpf_Shift \
  -t ../../../1.reference/norm/gene.norm.txt \
  -r ../05.merge/RIBO_merged.txt \
  -p 45 \
  -m 50 \
  --tis 5 \
  --tts 5 \
  -o rpf_shift
```
