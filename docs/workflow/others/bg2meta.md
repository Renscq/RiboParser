# bedGraph helpers

bedGraph helpers are used for converting interval-style signal tracks to nucleotide-level tables and for smoothing RPM density tracks.


### Input bedGraph format

Both commands expect a tab-delimited bedGraph-like file without a header:

```text
chrom    start    end    value
```

Coordinates follow the standard 0-based, half-open interval style used by bedGraph:

```text
chr1    100    103    2.5
```

represents positions `100-101`, `101-102`, and `102-103` with value `2.5`.

### `bg2meta`

#### Function

Convert bedGraph intervals to single-nucleotide rows for downstream plotting or meta-position analysis.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input bedGraph file. |
| `-o` | Output TXT/TSV file. |
| `-v`, `--version` | Print script version. |

#### Example

```bash
bg2meta \
  -i sample.plus.bedgraph \
  -o sample.plus.meta.txt
```

Example input:

```text
chr1    100    103    2.5
```

Example output:

```text
chr1    101    2.5
chr1    102    2.5
chr1    103    2.5
```
