#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Collapse duplicates and construct alternative-start families.
# Input: Validated gene-level candidates.
# Output: Primary families and member mappings.

"""Collapse duplicates and construct alternative-start families."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from typing import Any

from .config import PRIMARY_SELECTION_RULE
from .models import Candidate, ClusterSummary, ExactRepresentative, Family, _SuffixState, _UnionFind
from .output import _ShardWriter


class _FamilyMixin:
    """Provide exact-deduplication and family-construction methods."""

    @staticmethod
    def _category_rank(category: bytes) -> int:
        if category in {b"annotated_ORF", b"annotated_mORF"}:
            return 0
        if category in {
            b"uORF",
            b"dORF",
            b"iORF",
            b"same_frame_iORF",
            b"emORF",
            b"overlap_uORF",
            b"overlap_dORF",
        }:
            return 1
        if category == b"lncORF":
            return 2
        if category == b"other_ORF":
            return 3
        return 4

    def _exact_rank(self, candidate: Candidate) -> tuple[Any, ...]:
        return (
            self._category_rank(candidate.category),
            -candidate.transcript_length,
            candidate.input_order,
            candidate.orf_id,
        )

    def _primary_rank(
        self,
        representative: ExactRepresentative,
    ) -> tuple[Any, ...]:
        candidate = representative.primary
        protected = 0 if candidate.category in {b"annotated_ORF", b"annotated_mORF"} else 1
        kozak_rank = -candidate.kozak_score if candidate.kozak_score is not None else float("inf")
        return (
            protected,
            self.start_rank.get(candidate.start_codon, len(self.start_rank)),
            -candidate.aa_length,
            kozak_rank,
            -candidate.transcript_length,
            candidate.input_order,
            candidate.orf_id,
        )

    def _deduplicate_exact(
        self,
        candidates: list[Candidate],
        summary: ClusterSummary,
    ) -> list[ExactRepresentative]:
        grouped: dict[
            tuple[bytes, bytes, bytes, tuple[tuple[int, int], ...]],
            ExactRepresentative,
        ] = {}
        for candidate in candidates:
            key = candidate.exact_key
            representative = grouped.get(key)
            if representative is None:
                grouped[key] = ExactRepresentative(
                    primary=candidate,
                    members=[candidate],
                )
                continue
            representative.members.append(candidate)
            if self._exact_rank(candidate) < self._exact_rank(representative.primary):
                representative.primary = candidate
            summary.exact_duplicates_collapsed += 1
        representatives = list(grouped.values())
        representatives.sort(key=lambda item: item.primary.input_order)
        return representatives

    @staticmethod
    def _suffix_compatible(
        candidate: Candidate,
        registered_extent: int,
    ) -> bool:
        current_extent = candidate.start_extent()
        if candidate.strand == b"+":
            return registered_extent <= current_extent
        return registered_extent >= current_extent

    def _cluster_bucket_indexed(
        self,
        bucket_indices: list[int],
        representatives: list[ExactRepresentative],
        union_find: _UnionFind,
    ) -> None:
        """Union one stop bucket with an indexed splice-suffix algorithm."""
        ordered = sorted(
            bucket_indices,
            key=lambda index: (
                -representatives[index].primary.nt_length,
                representatives[index].primary.input_order,
            ),
        )
        states: dict[tuple[Any, ...], _SuffixState] = {}

        for index in ordered:
            candidate = representatives[index].primary
            own_key = candidate.path_key(0)
            state = states.get(own_key)
            if state is None:
                state = _SuffixState()
                states[own_key] = state
            else:
                if state.anchor is not None:
                    union_find.union(index, state.anchor)
                if state.pending:
                    retained_pending: list[tuple[int, int]] = []
                    for pending_index, pending_extent in state.pending:
                        if self._suffix_compatible(candidate, pending_extent):
                            union_find.union(index, pending_index)
                        else:
                            retained_pending.append((pending_index, pending_extent))
                    state.pending = retained_pending
            if state.anchor is None:
                state.anchor = index

            for block_index in range(1, len(candidate.oriented_blocks)):
                suffix_key = candidate.path_key(block_index)
                suffix_state = states.get(suffix_key)
                if suffix_state is None:
                    suffix_state = _SuffixState()
                    states[suffix_key] = suffix_state
                suffix_state.pending.append((index, candidate.suffix_extent(block_index)))

    def _cluster_representatives(
        self,
        representatives: list[ExactRepresentative],
    ) -> list[Family]:
        if not representatives:
            return []
        union_find = _UnionFind(len(representatives))
        buckets: dict[
            tuple[bytes, bytes, bytes, int, bytes],
            list[int],
        ] = defaultdict(list)
        for index, representative in enumerate(representatives):
            candidate = representative.primary
            buckets[
                (
                    candidate.chrom,
                    candidate.strand,
                    candidate.source_strand,
                    candidate.stop_boundary,
                    candidate.stop_codon,
                )
            ].append(index)
        for bucket_indices in buckets.values():
            if len(bucket_indices) > 1:
                self._cluster_bucket_indexed(
                    bucket_indices,
                    representatives,
                    union_find,
                )

        components: dict[int, list[ExactRepresentative]] = defaultdict(list)
        for index, representative in enumerate(representatives):
            components[union_find.find(index)].append(representative)

        families: list[Family] = []
        for component in components.values():
            component.sort(key=lambda item: item.primary.input_order)
            primary = min(component, key=self._primary_rank)
            families.append(Family(primary=primary, representatives=component))
        families.sort(key=lambda item: item.primary.primary.input_order)
        return families

    def _ordered_start_codons(
        self,
        members: Sequence[Candidate],
    ) -> bytes:
        codons = {member.start_codon for member in members}
        return b",".join(
            sorted(
                codons,
                key=lambda codon: (
                    self.start_rank.get(codon, len(self.start_rank)),
                    codon,
                ),
            )
        )

    @staticmethod
    def _family_type(
        representative_count: int,
        exact_duplicate_count: int,
    ) -> bytes:
        has_alt = representative_count > 1
        has_exact = exact_duplicate_count > 0
        if has_alt and has_exact:
            return b"MIXED"
        if has_alt:
            return b"ALTERNATIVE_START"
        if has_exact:
            return b"TRANSCRIPT_DUPLICATE"
        return b"SINGLETON"

    @staticmethod
    def _all_members(family: Family) -> list[Candidate]:
        members = [
            member for representative in family.representatives for member in representative.members
        ]
        members.sort(key=lambda item: item.input_order)
        return members

    def _write_family(
        self,
        family: Family,
        local_family_number: int,
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> None:
        primary_candidate = family.primary.primary
        all_members = self._all_members(family)
        representative_count = len(family.representatives)
        exact_duplicate_count = sum(
            len(representative.members) - 1 for representative in family.representatives
        )
        alt_start_count = representative_count - 1
        collapsed_count = len(all_members) - 1

        fields = list(primary_candidate.fields)
        fields.extend([b""] * (self.plan.family_count - len(fields)))
        family_lookup = self.plan.family_lookup
        fields[family_lookup[b"family_id"]] = str(local_family_number).encode("ascii")
        fields[family_lookup[b"family_role"]] = b"PRIMARY"
        fields[family_lookup[b"family_type"]] = self._family_type(
            representative_count,
            exact_duplicate_count,
        )
        fields[family_lookup[b"family_size"]] = str(len(all_members)).encode("ascii")
        fields[family_lookup[b"family_representative_count"]] = str(representative_count).encode(
            "ascii"
        )
        fields[family_lookup[b"family_transcript_count"]] = str(
            len({member.transcript_id for member in all_members})
        ).encode("ascii")
        fields[family_lookup[b"family_exact_duplicate_count"]] = str(exact_duplicate_count).encode(
            "ascii"
        )
        fields[family_lookup[b"family_alt_start_count"]] = str(alt_start_count).encode("ascii")
        fields[family_lookup[b"family_collapsed_count"]] = str(collapsed_count).encode("ascii")
        fields[family_lookup[b"family_categories"]] = b",".join(
            sorted({member.category for member in all_members})
        )
        fields[family_lookup[b"family_start_codons"]] = self._ordered_start_codons(all_members)
        fields[family_lookup[b"primary_selection_rule"]] = PRIMARY_SELECTION_RULE
        writer.write_family(fields)

        primary_id = family.primary.primary.orf_id
        ordered_representatives = [family.primary] + [
            representative
            for representative in family.representatives
            if representative is not family.primary
        ]
        family_number_bytes = str(local_family_number).encode("ascii")
        for representative in ordered_representatives:
            representative_id = representative.primary.orf_id
            ordered_members = [representative.primary] + [
                member for member in representative.members if member is not representative.primary
            ]
            for member in ordered_members:
                if member.orf_id == primary_id:
                    family_role = b"PRIMARY"
                    collapse_reason = b"PRIMARY"
                elif member.orf_id == representative_id:
                    family_role = b"COLLAPSED"
                    collapse_reason = b"ALTERNATIVE_START"
                else:
                    family_role = b"COLLAPSED"
                    collapse_reason = b"EXACT_TRANSCRIPT_DUPLICATE"
                writer.write_member(
                    (
                        family_number_bytes,
                        primary_id,
                        representative_id,
                        member.orf_id,
                        member.transcript_id,
                        member.category,
                        member.start_codon,
                        str(member.aa_length).encode("ascii"),
                        str(member.transcript_length).encode("ascii"),
                        family_role,
                        collapse_reason,
                    )
                )

        summary.family_count += 1
        summary.primary_categories[
            primary_candidate.category.decode("utf-8", errors="replace")
        ] += 1
        summary.alt_starts_collapsed += alt_start_count
        if len(all_members) == 1:
            summary.singleton_families += 1
        else:
            summary.multi_member_families += 1
