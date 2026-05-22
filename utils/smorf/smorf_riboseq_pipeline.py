#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Rensc
date: 2026-05-23

Main pipeline functions for smORF Ribo-seq evidence analysis.
"""

from typing import Dict, List, Tuple

import pandas as pd

from .smorf_riboseq_constants import DensityTrack, EvidenceThresholds
from .smorf_riboseq_density import load_chrom_density, scan_density_max_positions
from .smorf_riboseq_io import eprint, get_orf_blocks
from .smorf_riboseq_metrics import (
    calculate_pausing,
    calculate_periodicity,
    calculate_release,
    classify_coverage_shape,
    classify_translation_evidence,
    quantify_profile,
)
from .smorf_riboseq_profile import (
    extract_downstream_profile,
    extract_transcript_profile,
    get_coding_nt_length,
    make_codon_profile,
)


def process_one_track(
    orf_table: pd.DataFrame,
    track: DensityTrack,
    chrom_sizes: Dict[str, int],
    genepred_blocks: Dict[str, Tuple[List[int], List[int]]],
    thresholds: EvidenceThresholds,
    args,
) -> pd.DataFrame:
    """Process one density track and return an ORF-level evidence table."""
    eprint(f"[smORFEvidence] Processing sample={track.sample}, strand={track.strand}, file={track.path}")

    if not chrom_sizes:
        eprint("[smORFEvidence] No chromosome sizes provided. Scanning density file to infer chromosome sizes.")
        chrom_sizes = scan_density_max_positions(track.path, track.file_format)

    if not chrom_sizes:
        raise ValueError("Failed to determine chromosome sizes. Please provide --chrom-sizes.")

    if track.strand in {"+", "-"}:
        sub_table = orf_table[orf_table["strand"].eq(track.strand)].copy()
    else:
        sub_table = orf_table.copy()

    results = []
    grouped = sub_table.groupby("chrom", sort=False)

    for chrom, chrom_orfs in grouped:
        if chrom not in chrom_sizes:
            eprint(f"[Warning] Chromosome not found in chrom sizes or density: {chrom}. Skipped.")
            continue

        chrom_size = int(chrom_sizes[chrom])
        eprint(f"[smORFEvidence] Loading chromosome {chrom} ({chrom_size:,} bp), ORFs={len(chrom_orfs):,}")
        density = load_chrom_density(track.path, track.file_format, chrom, chrom_size)

        for idx, (_, row) in enumerate(chrom_orfs.iterrows(), start=1):
            if args.progress_every > 0 and idx % args.progress_every == 0:
                eprint(f"[smORFEvidence] {track.sample}:{chrom}: processed {idx:,}/{len(chrom_orfs):,} ORFs")

            result = score_one_orf(
                row=row,
                density=density,
                track=track,
                genepred_blocks=genepred_blocks,
                thresholds=thresholds,
                args=args,
            )
            results.append(result)

        del density

    return pd.DataFrame(results)


def score_one_orf(row, density, track, genepred_blocks, thresholds, args) -> dict:
    """Score one ORF using one P-site density track."""
    starts, ends = get_orf_blocks(row, genepred_blocks)

    nt_profile = extract_transcript_profile(
        density=density,
        starts=starts,
        ends=ends,
        strand=row["strand"],
    )

    if len(nt_profile) > int(row["nt_length"]):
        nt_profile = nt_profile[:int(row["nt_length"])]

    coding_nt_length = get_coding_nt_length(row, len(nt_profile))
    codon_profile = make_codon_profile(nt_profile, coding_nt_length)

    downstream_nt_profile = extract_downstream_profile(
        density=density,
        genomic_start=int(row["genomic_start"]),
        genomic_end=int(row["genomic_end"]),
        strand=str(row["strand"]),
        nt_window=args.post_stop_codons * 3,
    )

    quant = quantify_profile(nt_profile, codon_profile)
    periodicity = calculate_periodicity(nt_profile, coding_nt_length, thresholds)
    pausing = calculate_pausing(codon_profile, thresholds, args.pseudocount)
    release = calculate_release(
        codon_profile=codon_profile,
        downstream_nt_profile=downstream_nt_profile,
        post_stop_codons=args.post_stop_codons,
        thresholds=thresholds,
        pseudocount=args.pseudocount,
    )
    coverage_shape = classify_coverage_shape(codon_profile, thresholds)
    evidence = classify_translation_evidence(quant, periodicity, release, coverage_shape, thresholds)

    result = {
        "sample": track.sample,
        "orf_id": row["orf_id"],
        "gene_id": row["gene_id"],
        "transcript_id": row["transcript_id"],
        "chrom": row["chrom"],
        "strand": row["strand"],
        "category": row["category"] if "category" in row.index else ".",
        "genomic_start": int(row["genomic_start"]),
        "genomic_end": int(row["genomic_end"]),
        "nt_length": int(row["nt_length"]),
        "coding_nt_length": int(coding_nt_length),
        "coding_codon_count": int(len(codon_profile)),
    }

    result.update(quant)
    result.update(periodicity)
    result.update(pausing)
    result.update(release)
    result.update(coverage_shape)
    result["translation_evidence"] = evidence

    return result


def run_riboseq_evidence(args) -> pd.DataFrame:
    """Run the full smORF Ribo-seq evidence workflow."""
    from .smorf_riboseq_io import (
        read_chrom_sizes,
        read_density_list,
        read_genepred,
        read_orf_table,
    )

    eprint("[smORFEvidence] Reading ORF table.")
    orf_table = read_orf_table(args.orf_table, args.coord_mode)
    eprint(f"[smORFEvidence] Loaded PASS ORFs: {len(orf_table):,}")

    eprint("[smORFEvidence] Reading chromosome sizes.")
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)

    eprint("[smORFEvidence] Reading genePred blocks.")
    genepred_blocks = read_genepred(args.genepred, args.coord_mode)
    if genepred_blocks:
        eprint(f"[smORFEvidence] Loaded genePred blocks: {len(genepred_blocks):,}")

    eprint("[smORFEvidence] Reading P-site density files.")
    tracks = read_density_list(args)
    thresholds = build_thresholds(args)

    all_results = []
    for track in tracks:
        result = process_one_track(
            orf_table=orf_table,
            track=track,
            chrom_sizes=chrom_sizes,
            genepred_blocks=genepred_blocks,
            thresholds=thresholds,
            args=args,
        )
        all_results.append(result)

    if all_results:
        output_table = pd.concat(all_results, axis=0, ignore_index=True)
    else:
        output_table = pd.DataFrame()

    if not args.keep_no_evidence and not output_table.empty:
        output_table = output_table[output_table["translation_evidence"].ne("NoEvidence")].copy()

    return output_table


def build_thresholds(args) -> EvidenceThresholds:
    """Build evidence threshold object from CLI arguments."""
    return EvidenceThresholds(
        min_rpf_sum=args.min_rpf_sum,
        min_covered_codon=args.min_covered_codon,
        min_coverage_ratio=args.min_coverage_ratio,
        strong_periodicity=args.strong_periodicity,
        moderate_periodicity=args.moderate_periodicity,
        strong_start_pause=args.strong_start_pause,
        moderate_start_pause=args.moderate_start_pause,
        strong_stop_pause=args.strong_stop_pause,
        moderate_stop_pause=args.moderate_stop_pause,
        strong_release=args.strong_release,
        moderate_release=args.moderate_release,
        uniform_coverage_ratio=args.uniform_coverage_ratio,
        uniform_gini=args.uniform_gini,
        uniform_max_to_mean=args.uniform_max_to_mean,
        skewed_max_to_mean=args.skewed_max_to_mean,
        skewed_top_fraction=args.skewed_top_fraction,
        disperse_coverage_ratio=args.disperse_coverage_ratio,
    )
