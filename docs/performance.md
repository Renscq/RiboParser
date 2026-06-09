# 6 Performance

The original README assessed the workflow on a CentOS 7 system using 12 threads across RNA-seq and Ribo-seq datasets from three species.

| Species | Dataset | Library | Sample number | Sample size | Index building elapsed time | Index disk usage | Preprocessing & alignment elapsed time | Preprocessing disk usage | RiboParser elapsed time | RiboParser disk usage |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| *S. cerevisiae* | GSE67387 | Ribo-seq | 6 | 32 G | 38 s | 357 M | 43 m 21 s | 30 G | 1 h 26 m 23 s | 3.6 G |
| *S. cerevisiae* | GSE67387 | RNA-seq | 6 | 17 G | 38 s | 357 M | 32 m 52 s | 26 G | 37 m 42 s | 2.8 G |
| *M. musculus* | GSE114064 | Ribo-seq | 6 | 43 G | 59 m 8 s | 36 G | 50 m 27 s | 31 G | 32 m 3 s | 7.9 G |
| *M. musculus* | GSE114064 | RNA-seq | 6 | 60 G | 59 m 8 s | 36 G | 4 h 14 m 45 s | 62 G | 29 m 30 s | 7.6 G |
| *H. sapiens* | GSE131650 | Ribo-seq | 6 | 42 G | 1 h 55 m 56 s | 44 G | 2 h 11 m 42 s | 29 G | 2 h 18 m 57 s | 14 G |
| *H. sapiens* | GSE131650 | RNA-seq | 6 | 54 G | 1 h 55 m 56 s | 44 G | 1 h 15 m 15 s | 30 G | 40 m 35 s | 11 G |

## System recommendation

For optimal performance, deploy RiboParser on Linux-based systems. The benchmark was tested on Ubuntu 20.04 LTS / CentOS 7-like environments.
