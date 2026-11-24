#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-11-24
# Version: dev001
# Function: Exercise splice-aware smORF family construction.
# Input: Small synthetic models and the AT1G79075 regression geometry.
# Output: Unit-test assertions for family identities and shared coding paths.

"""Test family construction on strand-aware, spliced candidate coordinates."""

import unittest

from utils.smorf.cluster.annotation import TranscriptMeta
from utils.smorf.cluster.family import _FamilyMixin
from utils.smorf.cluster.models import Candidate, ClusterSummary


class FamilyBuilder(_FamilyMixin):
    """Use the production family methods with the canonical start order."""

    start_rank = {b"ATG": 0}


def candidate(identifier, blocks, order, strand=b"+", transcript=b"AT1G79075.1"):
    """Build one complete candidate with scanner-style half-open exon bounds."""
    length = sum(end - start for start, end in blocks)
    meta = TranscriptMeta(transcript, b"AT1G79075", b"1", strand, 1000, order, False, ())
    return Candidate(
        fields=(), input_order=order, transcript_meta=meta, blocks=blocks,
        oriented_blocks=blocks if strand == b"+" else blocks[::-1], nt_length=length,
        aa_length=length // 3 - 1, start_codon=b"ATG", stop_codon=b"TAA",
        kozak_score=None, orf_id=identifier, transcript_id=transcript, gene_id=meta.gene_id,
        chrom=meta.chrom, strand=strand, source_strand=b"sense", category=b"lncORF",
        stop_boundary=blocks[-1][1] if strand == b"+" else blocks[0][0],
    )


def families(candidates):
    """Run the same exact deduplication and clustering used by the CLI."""
    builder = FamilyBuilder()
    return builder._cluster_representatives(builder._deduplicate_exact(candidates, ClusterSummary()))


class FamilyTests(unittest.TestCase):
    """Verify that alternate coding paths remain scoreable without duplicated ORFs."""

    def test_shared_suffix_does_not_merge_incompatible_scaffolds(self):
        """The 16 aa AT1G79075 suffix belongs to two incompatible longer paths."""
        long_blocks = ((29746615, 29746705),)
        shared_blocks = ((29746654, 29746705),)
        spliced_blocks = ((29746013, 29746025), (29746639, 29746705))
        for strand in (b"+", b"-"):
            with self.subTest(strand=strand):
                def orient(blocks):
                    """Mirror the observed intervals for a minus-strand regression."""
                    if strand == b"+":
                        return blocks
                    return tuple((40000000 - end, 40000000 - start) for start, end in blocks[::-1])

                result = families([
                    candidate(b"long", orient(long_blocks), 1, strand),
                    candidate(b"shared", orient(shared_blocks), 2, strand),
                    candidate(b"spliced", orient(spliced_blocks), 3, strand),
                    candidate(b"duplicate", orient(shared_blocks), 4, strand, b"AT1G79075.2"),
                ])
                by_id = {family.primary.primary.orf_id: family for family in result}
                self.assertEqual(set(by_id), {b"long", b"shared", b"spliced"})
                shared = by_id[b"shared"]
                self.assertEqual({m.orf_id for m in shared.primary.members}, {b"shared", b"duplicate"})
                self.assertEqual(shared.splice_path_role, b"SharedBody")
                self.assertEqual(set(shared.splice_scaffold_ids), {b"long", b"spliced"})
                for name in (b"long", b"spliced"):
                    self.assertEqual(by_id[name].splice_path_role, b"Branch")
                    self.assertTrue(by_id[name].splice_unique_positions)

    def test_nested_starts_still_form_one_family(self):
        """A single coding path retains the existing one-family representation."""
        result = families([candidate(b"long", ((100, 160),), 1), candidate(b"short", ((118, 160),), 2)])
        self.assertEqual(len(result), 1)
        self.assertEqual(len(result[0].representatives), 2)

    def test_split_stop_is_ordered_on_each_strand(self):
        """A split codon retains all three physical coordinates in coding order."""
        plus = candidate(b"plus", ((100, 131), (200, 202)), 1)
        minus = candidate(b"minus", ((100, 102), (200, 231)), 2, b"-")
        self.assertEqual(plus.stop_positions, (130, 200, 201))
        self.assertEqual(minus.stop_positions, (200, 101, 100))


if __name__ == "__main__":
    unittest.main()
