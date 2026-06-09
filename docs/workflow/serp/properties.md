# 4.9.3 SeRP properties

## Purpose

`serp_properties` summarizes properties of SeRP peaks or signal regions.

## Input files

| Input | Description |
|---|---|
| peak table | output from `serp_peak` or external peak file |
| annotation | optional gene/transcript annotation |
| signal file | optional signal file for intensity statistics |

## Parameters

| Parameter | Meaning |
|---|---|
| `-i / --input` | input peak or region table |
| `-a / --annotation` | annotation file |
| `-s / --signal` | signal file |
| `-o / --output` | output prefix |
| `--strand` | strand-aware assignment |
| `--distance-to-feature` | calculate distance to selected gene feature |
| `--summary` | generate summary statistics |

## Example

```bash
serp_properties \
  -i serp_peak.txt \
  -a gene.norm.txt \
  -s serp_signal.bedgraph \
  -o serp_properties
```

## Output files

| Output | Description |
|---|---|
| `serp_properties.txt` | peak property table |
| `serp_properties_summary.txt` | summary statistics |
| `serp_properties.log` | running log |

## Result interpretation

Use property summaries to evaluate peak length, intensity, annotation class, positional enrichment, and sample-specific signal patterns.
