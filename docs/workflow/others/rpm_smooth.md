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


### `rpm_smooth`

#### Function

Smooth RPM density in bedGraph files. The command first expands interval bins to single-nucleotide bins, fills missing positions with zero, smooths values chromosome by chromosome, clips negative smoothed values to zero, and writes a bedGraph-like output.

#### Parameters

| Parameter | Meaning |
|---|---|
| `-i` | Input bedGraph file. |
| `-o` | Output smoothed bedGraph file. |
| `-m` | Smoothing method: `SG` for Savitzky-Golay or `Gaus` for Gaussian. Default: `SG`. |
| `-w` | Window width for Savitzky-Golay smoothing. Default: `11`. |
| `-p` | Polynomial order for Savitzky-Golay smoothing. Default: `3`. |
| `-s` | Gaussian sigma. Default: `1`. |
| `-z` | Remove rows with smoothed RPM equal to zero. Default: keep zero rows. |

#### Savitzky-Golay example

```bash
rpm_smooth \
  -i sample.rpm.bedgraph \
  -m SG \
  -w 11 \
  -p 3 \
  -z \
  -o sample.rpm.SG.bedgraph
```

#### Gaussian example

```bash
rpm_smooth \
  -i sample.rpm.bedgraph \
  -m Gaus \
  -s 2 \
  -z \
  -o sample.rpm.Gaus.bedgraph
```

### Notes

- `rpm_smooth` can generate dense chromosome-wide single-nucleotide output. For large genomes, run it per chromosome or on filtered regions when possible.
- For Savitzky-Golay smoothing, the window width should be appropriate for the expected signal shape. Very small windows may under-smooth, while very large windows may flatten true peaks.
- Use `-z` when you want a compact bedGraph output without zero-coverage positions.
