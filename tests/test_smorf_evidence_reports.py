#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-11-24
# Version: dev001
# Function: Verify candidate-level smORF evidence and CDS overlap calls.
# Input: Small coding paths and synthetic P-site counts.
# Output: Assertions for priority, splice, and annotated-CDS evidence.

"""Test candidate-level evidence independently of family scoring."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from utils.smorf.cluster.annotation import AnnotatedORF, AnnotationIndex
from utils.smorf.density import ChromDensity
from utils.smorf.evidence.cds_overlap import (
    EvidenceContext,
    cds_masks,
    cds_sample_evidence,
    overlap_status,
)
from utils.smorf.evidence.reporting import (
    DensityProfiles,
    coding_positions,
    priority_sample_pass,
    splice_sample_pass,
    summarize_counts,
)


def cds(name, blocks, strand):
    """Create a genePred-style CDS with offsets aligned to genomic blocks."""
    return AnnotatedORF(
        name, name, b"gene", b"1", strand, blocks,
        sum(end - start for start, end in blocks),
        AnnotationIndex._block_offsets(blocks, strand),
    )


class CandidateEvidenceTests(unittest.TestCase):
    """Exercise stop-excluded P-sites and per-library support."""

    def test_short_candidates_need_three_distributed_frame_zero_codons(self):
        """A 19 aa ORF qualifies with three occupied codons across two segments."""
        counts = np.zeros((19, 3))
        counts[[0, 3, 10], 0] = [4, 3, 3]
        self.assertTrue(priority_sample_pass(summarize_counts(counts), 19))

        longer = np.pad(counts, ((0, 1), (0, 0)))
        self.assertFalse(priority_sample_pass(summarize_counts(longer), 20))
        longer[16, 0] = 1
        self.assertTrue(priority_sample_pass(summarize_counts(longer), 20))

    def test_separate_samples_cannot_supply_each_others_signal(self):
        """Depth, periodicity, and breadth must pass in the same sample."""
        deep_peak = np.zeros((20, 3))
        deep_peak[0, 0] = 10
        sparse = np.zeros_like(deep_peak)
        sparse[[4, 8, 14], 0] = 1
        self.assertFalse(priority_sample_pass(summarize_counts(deep_peak), 20))
        self.assertFalse(priority_sample_pass(summarize_counts(sparse), 20))
        self.assertTrue(priority_sample_pass(summarize_counts(deep_peak + sparse), 20))

        off_frame = deep_peak + sparse
        off_frame[9, 1] = 10
        self.assertFalse(priority_sample_pass(summarize_counts(off_frame), 20))

    def test_split_stop_is_excluded_on_both_strands(self):
        """The three stop bases are removed even when an intron divides them."""
        plus = coding_positions(((100, 110), (200, 202)), "+", 9)
        minus = coding_positions(((100, 102), (200, 210)), "-", 9)
        self.assertEqual(plus.tolist(), list(range(100, 109)))
        self.assertEqual(minus.tolist(), list(range(209, 200, -1)))

    def test_half_open_density_tracks_preserve_phase_and_strand_sign(self):
        """Density outside a half-open interval must not leak into an adjacent phase."""
        density = ChromDensity(
            chrom="1", starts=np.array([100, 104]), ends=np.array([102, 105]),
            values=np.array([-2, -4]), max_end=105,
        )
        positions = coding_positions(((100, 109),), "-", 6)
        counts = DensityProfiles._counts_at(density, positions)
        np.testing.assert_array_equal(counts, [[0, 0, 0], [0, 4, 0]])

    def test_shared_density_cannot_support_a_specific_splice_path(self):
        """Only diagnostic frame-zero positions from the selected ORF apply."""
        positions = coding_positions(((100, 130),), "+", 27)
        counts = np.zeros((9, 3))
        counts[[5, 8], 0] = 50
        self.assertEqual(splice_sample_pass(counts, positions, (100, 103), 3.0), (0.0, 0, False))
        counts[[0, 1], 0] = [2, 1]
        self.assertEqual(splice_sample_pass(counts, positions, (100, 103), 3.0), (3.0, 2, True))
        self.assertEqual(splice_sample_pass(counts, positions, (97,), 3.0), (0.0, 0, False))


class AnnotatedCDSTests(unittest.TestCase):
    """Distinguish candidate-specific signal from annotated-CDS signal."""

    def test_all_annotated_isoforms_mask_coding_positions(self):
        """A second CDS can explain the apparent extension of the first."""
        positions = coding_positions(((100, 130), (200, 233)), "+", 60)
        parents = (cds(b"one", ((100, 130),), b"+"), cds(b"two", ((200, 233),), b"+"))
        masks = cds_masks(positions, "+", parents)
        self.assertEqual(masks.parent_ids, ("one", "two"))
        self.assertEqual(int(masks.shared_phase0.sum()), 20)
        self.assertFalse(masks.outside_codons.any())

    def test_competing_annotated_frame_and_partial_codon(self):
        """An overlapping frame contributes host signal and masks partial codons."""
        positions = coding_positions(((100, 130),), "-", 27)
        masks = cds_masks(positions, "-", (cds(b"host", ((108, 128),), b"-"),))
        self.assertEqual(masks.parent_ids, ("host",))
        self.assertEqual(int(masks.shared_phase0.sum()), 0)
        self.assertEqual(int(masks.outside_codons.sum()), 1)
        counts = np.ones((9, 3))
        self.assertGreater(cds_sample_evidence(counts, masks).competing_phase0_rpf, 0)

    def test_shared_only_and_initiation_are_distinct_evidence(self):
        """A shared body stays ambiguous unless initiation is independently observed."""
        positions = coding_positions(((100, 136),), "+", 33)
        masks = cds_masks(positions, "+", (cds(b"host", ((97, 136),), b"+"),))
        counts = np.zeros((11, 3))
        counts[[0, 4, 8], 0] = 10
        sample = cds_sample_evidence(counts, masks)
        self.assertEqual(overlap_status([sample], 1, False), "SharedCDSSignal")
        self.assertEqual(overlap_status([sample], 1, True), "AlternativeStartSupported")
        self.assertEqual(overlap_status([sample], 2, True), "SharedCDSSignal")

    def test_extension_can_be_supported_without_initiation_data(self):
        """Three outside codons can support translated extent independently."""
        positions = coding_positions(((100, 160),), "+", 57)
        masks = cds_masks(positions, "+", (cds(b"host", ((121, 160),), b"+"),))
        self.assertEqual(int(masks.outside_codons.sum()), 7)
        counts = np.zeros((19, 3))
        counts[[1, 3, 5], 0] = 5
        sample = cds_sample_evidence(counts, masks)
        self.assertEqual(overlap_status([sample], 1, False), "IndependentExtentSupported")

    def test_external_evidence_matches_exact_starts_and_preserves_source(self):
        """Keep exact LTM sites, transcript junctions, and ORF peptides separate."""
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            annotation = root / "transcripts.genePred"
            annotation.write_text("tx\t1\t+\t100\t160\t100\t160\t1\t100,\t160,\n", encoding="utf-8")
            tis = root / "tis.tsv"
            tis.write_text(
                "chrom\tstrand\tstart0\tltm_reads\tltm_minus_chx\tgene_id\tsource\ttissue\n"
                "1\t+\t100\t30\t1\tgene\tTIS-1\tseedling\n"
                "1\t+\t103\t40\t1\tother_gene\tTIS-1\tseedling\n"
                "1\t+\t103\t50\t-1\tgene\tTIS-1\tseedling\n"
                "1\t+\t106\t60\t2\tgene\tTIS-2\tseedling\n",
                encoding="utf-8",
            )
            junctions = root / "junctions.tsv"
            junctions.write_text(
                "transcript_id\tsupporting_reads\tsource\ntx\t7\tRNA-1\n", encoding="utf-8",
            )
            peptides = root / "peptides.tsv"
            peptides.write_text("orf_id\tevidence_id\norf\tMS-1\n", encoding="utf-8")
            context = EvidenceContext(annotation, tis, junctions, peptides)

            for start, exact, alternative in ((100, 1, 0), (103, 0, 0), (106, 1, 1)):
                with self.subTest(start=start):
                    blocks = ((start, 139),)
                    positions = coding_positions(blocks, "+", 136 - start)
                    _, metadata = context.describe("1", "+", blocks, positions, {"gene"}, {"tx"}, {"orf"})
                    self.assertEqual(metadata["tis_exact_start"], exact)
                    self.assertEqual(metadata["tis_alternative_start"], alternative)
                    self.assertEqual(metadata["junction_context"], "tx|RNA-1|7")
                    self.assertEqual(metadata["peptide_evidence"], "MS-1")
            self.assertEqual(metadata["tis_source"], "tis.tsv|TIS-2|seedling")


if __name__ == "__main__":
    unittest.main()
