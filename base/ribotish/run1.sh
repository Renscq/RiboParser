grep -vE 'Annotated|Internal|Truncated|Extended|CDSFrameOverlap|Known' /mnt/t64/rensc/usda110/nodules-ribo/gmx/sorf/ribotish/gmx4_all.txt | sed 's/gene-//g' > test.txt
