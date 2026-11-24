# RiboParser 0.3.16

- Split incompatible smORF splice paths into separate start families. Record a shared suffix once and use branch-specific P-sites before exporting a splice-path call.
- Add a priority table and stop-excluded, per-sample evidence for the selected ORF. The reliable table and genePred retain their schemas for downstream quantification.
- Add an optional audit for ORFs overlapping annotated CDSs. The report distinguishes shared CDS signal, signal at codons outside all overlapping CDS isoforms, and support at an external initiation site. RNA junction and ORF-specific peptide evidence can be supplied as separate inputs.
