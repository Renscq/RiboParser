# 4.9 SeRP analysis

SeRP-related modules are placed in the main workflow because they are specialized Ribo-seq/RiboParser analyses rather than general helper scripts.

## Modules

| Section | Command | Purpose |
|---|---|---|
| 4.9.1 | `serp_overlap` | overlap analysis for SeRP signals or regions |
| 4.9.2 | `serp_peak` | SeRP peak detection |
| 4.9.3 | `serp_properties` | peak or signal property summary |

## General input types

- Ribo-seq or SeRP signal tracks
- BED/bedGraph-like regions
- peak tables
- transcript or gene annotations
