#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev001
# Function: Parse and filter candidate ORFs.
# Input: Scanner rows, annotation index, and cluster rules.
# Output: Validated candidates with annotation-overlap filtering.

"""Parse and filter candidate ORFs."""

from __future__ import annotations

import heapq
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from .annotation import AnnotatedORF, AnnotatedOverlap, AnnotationIndex, TranscriptMeta
from .config import _BUFFER_BYTES, ANNOTATED_CATEGORIES, STOP_CODONS
from .models import Candidate, ClusterSummary
from .output import _ShardWriter


class _FilteringMixin:
    """Provide scanner parsing and annotation-aware filtering methods."""

    @staticmethod
    def _parse_int(value: bytes) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _parse_blocks(
        starts_value: bytes,
        ends_value: bytes,
        exon_count: int | None,
    ) -> tuple[tuple[int, int], ...] | None:
        try:
            starts = AnnotationIndex._parse_int_list(starts_value)
            ends = AnnotationIndex._parse_int_list(ends_value)
        except ValueError:
            return None
        if not starts or len(starts) != len(ends) or exon_count != len(starts):
            return None
        blocks = tuple(zip(starts, ends))
        if any(end <= start for start, end in blocks):
            return None
        return blocks

    @staticmethod
    def _unique_reasons(reasons: Iterable[bytes]) -> list[bytes]:
        return list(dict.fromkeys(reasons))

    def _annotate_kozak(
        self,
        fields: list[bytes],
        start_codon: bytes,
    ) -> float | None:
        scorer = self.kozak_scorer
        if scorer is None:
            fields[self.i_kozak_name] = b"none"
            fields[self.i_kozak_score] = b"NA"
            fields[self.i_kozak_level] = b"NA"
            fields[self.i_kozak_valid] = b"NA"
            return None
        sequence = fields[self.i_kozak_seq] if self.i_kozak_seq >= 0 else b""
        start_index = fields[self.i_kozak_start] if self.i_kozak_start >= 0 else b""
        score, level, valid_ratio = scorer.score(
            sequence,
            start_index,
            start_codon,
        )
        fields[self.i_kozak_name] = scorer.spec.name
        fields[self.i_kozak_score] = b"NA" if score is None else f"{score:.6f}".encode("ascii")
        fields[self.i_kozak_level] = level
        fields[self.i_kozak_valid] = f"{valid_ratio:.6f}".encode("ascii")
        return score

    def _parse_candidate(
        self,
        raw_fields: list[bytes],
        input_order: int,
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> Candidate | None:
        summary.input_orfs += 1
        if len(raw_fields) != self.plan.input_count:
            fields = raw_fields[: self.plan.input_count]
            fields.extend([b""] * (self.plan.filter_count - len(fields)))
            fields[self.i_structure] = b"FAIL"
            fields[self.i_filter_status] = b"FAIL"
            fields[self.i_filter_reason] = b"invalid_column_count"
            summary.basic_removed += 1
            summary.removal_reasons["invalid_column_count"] += 1
            writer.write_removed(fields)
            return None

        fields = raw_fields
        if len(fields) < self.plan.filter_count:
            fields.extend([b""] * (self.plan.filter_count - len(fields)))

        reasons: list[bytes] = []
        structure_reasons: list[bytes] = []
        source_strand = fields[self.i_source_strand].strip().lower()
        priority = fields[self.i_priority].strip()
        completeness = fields[self.i_completeness].strip()
        category = fields[self.i_category].strip()
        start_codon = fields[self.i_start_codon].strip().upper().replace(b"U", b"T")
        stop_codon = fields[self.i_stop_codon].strip().upper().replace(b"U", b"T")
        transcript_id = fields[self.i_transcript_id].strip()
        gene_id = fields[self.i_gene_id].strip()
        chrom = fields[self.i_chrom].strip()
        strand = fields[self.i_strand].strip()

        if source_strand not in {b"sense", b"antisense"}:
            reasons.append(b"invalid_source_strand")
        elif self.config.require_sense and source_strand != b"sense":
            reasons.append(b"non_sense_strand")
        if not priority:
            reasons.append(b"missing_priority")
        if not completeness:
            reasons.append(b"missing_completeness")
        elif completeness != b"complete":
            reasons.append(b"incomplete_orf:" + completeness)
        if not category:
            reasons.append(b"missing_category")
        else:
            if category in self.config.remove_categories:
                reasons.append(b"removed_category:" + category)
            if self.config.keep_categories and category not in self.config.keep_categories:
                reasons.append(b"category_not_allowed:" + category)
        if not start_codon:
            reasons.append(b"missing_start_codon")
        elif self.config.keep_start_codons and start_codon not in self.config.keep_start_codons:
            reasons.append(b"start_codon_not_allowed:" + start_codon)

        aa_length = self._parse_int(fields[self.i_aa_length])
        nt_length = self._parse_int(fields[self.i_nt_length])
        exon_count = self._parse_int(fields[self.i_exon_count])
        if aa_length is None:
            reasons.append(b"invalid_aa_length")
        else:
            if aa_length < self.config.min_aa:
                reasons.append(b"too_short")
            if aa_length > self.config.max_aa:
                reasons.append(b"too_long")
        if nt_length is None or nt_length <= 0:
            structure_reasons.append(b"invalid_nt_length")
        elif aa_length is not None and nt_length != aa_length * 3 + 3:
            structure_reasons.append(b"inconsistent_orf_length")
        if completeness == b"complete" and stop_codon not in STOP_CODONS:
            structure_reasons.append(b"invalid_stop_codon")

        if self.i_ambiguous >= 0:
            ambiguous = self._parse_int(fields[self.i_ambiguous])
            if ambiguous is None or ambiguous < 0:
                structure_reasons.append(b"invalid_ambiguous_codon_count")
            elif ambiguous > self.config.max_ambiguous_codons:
                structure_reasons.append(b"too_many_ambiguous_codons")

        blocks = self._parse_blocks(
            fields[self.i_exon_starts],
            fields[self.i_exon_ends],
            exon_count,
        )
        if blocks is None:
            structure_reasons.append(b"invalid_exon_blocks")
        elif nt_length is not None and sum(end - start for start, end in blocks) != nt_length:
            structure_reasons.append(b"inconsistent_block_length")

        transcript_meta = self.annotation.transcripts.get(transcript_id)
        if transcript_meta is None:
            reasons.append(b"transcript_not_in_annotation")
        else:
            if gene_id != transcript_meta.gene_id:
                reasons.append(b"annotation_gene_mismatch")
            if chrom != transcript_meta.chrom:
                reasons.append(b"annotation_chrom_mismatch")
            expected_strand = transcript_meta.strand
            if source_strand == b"antisense":
                expected_strand = b"-" if expected_strand == b"+" else b"+"
            if strand != expected_strand:
                reasons.append(b"annotation_strand_mismatch")

        kozak_score = self._annotate_kozak(fields, start_codon)
        if self.config.kozak_spec is not None and self.config.min_kozak_score > 0:
            if kozak_score is None:
                reasons.append(b"insufficient_kozak_context")
            elif kozak_score < self.config.min_kozak_score:
                reasons.append(b"weak_kozak_pwm")

        reasons.extend(structure_reasons)
        unique_reasons = self._unique_reasons(reasons)
        fields[self.i_structure] = b"FAIL" if structure_reasons else b"PASS"
        if unique_reasons:
            fields[self.i_filter_status] = b"FAIL"
            fields[self.i_filter_reason] = b";".join(unique_reasons)
            summary.basic_removed += 1
            summary.removal_reasons.update(
                reason.decode("utf-8", errors="replace") for reason in unique_reasons
            )
            writer.write_removed(fields)
            return None

        fields[self.i_filter_status] = b"PASS"
        fields[self.i_filter_reason] = b"PASS"
        summary.basic_passed += 1
        assert transcript_meta is not None
        assert blocks is not None
        assert nt_length is not None
        assert aa_length is not None
        oriented_blocks = blocks if strand == b"+" else tuple(reversed(blocks))
        terminal_block = oriented_blocks[-1]
        stop_boundary = terminal_block[1] if strand == b"+" else terminal_block[0]
        return Candidate(
            fields=tuple(fields),
            input_order=input_order,
            transcript_meta=transcript_meta,
            blocks=blocks,
            oriented_blocks=oriented_blocks,
            nt_length=nt_length,
            aa_length=aa_length,
            start_codon=start_codon,
            stop_codon=stop_codon,
            kozak_score=kozak_score,
            orf_id=fields[self.i_orf_id],
            transcript_id=transcript_id,
            gene_id=gene_id,
            chrom=chrom,
            strand=strand,
            source_strand=source_strand,
            category=category,
            stop_boundary=stop_boundary,
        )

    def _iter_transcript_groups(
        self,
        input_path: str | Path,
        start: int,
        end: int,
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> Iterator[tuple[TranscriptMeta, list[Candidate]]]:
        """Stream and filter contiguous transcript blocks in one byte range."""
        current_transcript_id: bytes | None = None
        current_meta: TranscriptMeta | None = None
        current_candidates: list[Candidate] = []
        last_annotation_order = 0
        input_order = 0

        with Path(input_path).open("rb", buffering=_BUFFER_BYTES) as handle:
            handle.seek(start)
            while handle.tell() < end:
                line = handle.readline()
                if not line:
                    break
                input_order += 1
                raw_fields = line.rstrip(b"\r\n").split(b"\t")
                transcript_id = (
                    raw_fields[self.i_transcript_id].strip()
                    if len(raw_fields) > self.i_transcript_id
                    else b""
                )

                if current_transcript_id is not None and transcript_id != current_transcript_id:
                    if current_meta is not None:
                        if current_meta.annotation_order <= last_annotation_order:
                            raise ValueError(
                                "ORF message table is not ordered like the source "
                                "genePred annotation at transcript "
                                + current_transcript_id.decode(
                                    "utf-8",
                                    errors="replace",
                                )
                            )
                        last_annotation_order = current_meta.annotation_order
                        yield current_meta, current_candidates
                    current_candidates = []
                    current_transcript_id = None
                    current_meta = None

                if current_transcript_id is None and transcript_id:
                    current_transcript_id = transcript_id
                    current_meta = self.annotation.transcripts.get(transcript_id)

                candidate = self._parse_candidate(
                    raw_fields=raw_fields,
                    input_order=input_order,
                    summary=summary,
                    writer=writer,
                )
                if candidate is not None:
                    if current_meta is None:
                        current_meta = candidate.transcript_meta
                        current_transcript_id = candidate.transcript_id
                    current_candidates.append(candidate)

            if current_transcript_id is not None and current_meta is not None:
                if current_meta.annotation_order <= last_annotation_order:
                    raise ValueError(
                        "ORF message table is not transcript-contiguous or "
                        "annotation-ordered at transcript "
                        + current_transcript_id.decode(
                            "utf-8",
                            errors="replace",
                        )
                    )
                yield current_meta, current_candidates

    def _iter_gene_groups(
        self,
        transcript_groups: Iterator[tuple[TranscriptMeta, list[Candidate]]],
    ) -> Iterator[tuple[bytes, list[Candidate]]]:
        """Aggregate active genes with a last-order min-heap."""
        pending: dict[bytes, list[Candidate]] = {}
        close_heap: list[tuple[int, int, bytes]] = []

        def flush_before(order: int) -> Iterator[tuple[bytes, list[Candidate]]]:
            while close_heap and close_heap[0][0] < order:
                _last, _first, gene_id = heapq.heappop(close_heap)
                candidates = pending.pop(gene_id, None)
                if candidates is not None:
                    yield gene_id, candidates

        for transcript_meta, candidates in transcript_groups:
            yield from flush_before(transcript_meta.annotation_order)
            gene_id = transcript_meta.gene_id
            if gene_id not in pending:
                pending[gene_id] = []
                heapq.heappush(
                    close_heap,
                    (
                        self.annotation.gene_last_order[gene_id],
                        self.annotation.gene_first_order[gene_id],
                        gene_id,
                    ),
                )
            pending[gene_id].extend(candidates)
            if transcript_meta.annotation_order == self.annotation.gene_last_order[gene_id]:
                pending_candidates = pending.pop(gene_id)
                yield gene_id, pending_candidates

        while close_heap:
            _last, _first, gene_id = heapq.heappop(close_heap)
            candidates = pending.pop(gene_id, None)
            if candidates is not None:
                yield gene_id, candidates

    @staticmethod
    def _phase_compatible_overlap(
        candidate: Candidate,
        annotated_orf: AnnotatedORF,
    ) -> AnnotatedOverlap | None:
        """Require every shared segment to use one genomic coding phase."""
        candidate_offsets = AnnotationIndex._block_offsets(
            candidate.blocks,
            candidate.strand,
        )
        candidate_index = 0
        annotated_index = 0
        overlap_nt = 0
        overlap_codon = 0

        while candidate_index < len(candidate.blocks) and annotated_index < len(
            annotated_orf.blocks
        ):
            candidate_start, candidate_end = candidate.blocks[candidate_index]
            annotated_start, annotated_end = annotated_orf.blocks[annotated_index]
            shared_start = max(candidate_start, annotated_start)
            shared_end = min(candidate_end, annotated_end)

            if shared_end > shared_start:
                shared_length = shared_end - shared_start
                if candidate.strand == b"+":
                    candidate_offset = (
                        candidate_offsets[candidate_index] + shared_start - candidate_start
                    )
                    annotated_offset = (
                        annotated_orf.block_offsets[annotated_index]
                        + shared_start
                        - annotated_start
                    )
                else:
                    candidate_offset = (
                        candidate_offsets[candidate_index] + candidate_end - shared_end
                    )
                    annotated_offset = (
                        annotated_orf.block_offsets[annotated_index] + annotated_end - shared_end
                    )

                if (candidate_offset - annotated_offset) % 3 != 0:
                    return None

                overlap_nt += shared_length
                skip_to_codon = (-candidate_offset) % 3
                available = shared_length - skip_to_codon
                if available >= 3:
                    overlap_codon += available // 3

            if candidate_end <= annotated_end:
                candidate_index += 1
            if annotated_end <= candidate_end:
                annotated_index += 1

        if overlap_nt == 0 or overlap_codon == 0:
            return None

        candidate_length = sum(end - start for start, end in candidate.blocks)
        if candidate.blocks == annotated_orf.blocks:
            relation = b"matches_annotated_ORF"
        elif overlap_nt == candidate_length:
            relation = b"contained_in_annotated_ORF_same_frame"
        elif overlap_nt == annotated_orf.nt_length:
            relation = b"extends_annotated_ORF_same_frame"
        else:
            relation = b"partial_overlap_annotated_ORF_same_frame"

        return AnnotatedOverlap(
            annotated_orf=annotated_orf,
            relation=relation,
            overlap_nt=overlap_nt,
            overlap_codon=overlap_codon,
        )

    def _best_annotated_overlap(
        self,
        candidate: Candidate,
    ) -> AnnotatedOverlap | None:
        """Return the strongest global annotated-ORF overlap match."""
        if candidate.category in ANNOTATED_CATEGORIES:
            return None

        relation_rank = {
            b"matches_annotated_ORF": 0,
            b"contained_in_annotated_ORF_same_frame": 1,
            b"extends_annotated_ORF_same_frame": 2,
            b"partial_overlap_annotated_ORF_same_frame": 3,
        }
        best_match: AnnotatedOverlap | None = None
        best_rank: tuple[Any, ...] | None = None
        for annotated_orf in self.annotation.overlapping_annotated_orfs(candidate):
            match = self._phase_compatible_overlap(candidate, annotated_orf)
            if match is None:
                continue
            rank = (
                relation_rank[match.relation],
                annotated_orf.gene_id != candidate.gene_id,
                -match.overlap_codon,
                -match.overlap_nt,
                annotated_orf.annotated_orf_id,
            )
            if best_rank is None or rank < best_rank:
                best_rank = rank
                best_match = match
        return best_match

    def _remove_annotated_overlaps(
        self,
        candidates: list[Candidate],
        summary: ClusterSummary,
        writer: _ShardWriter,
    ) -> list[Candidate]:
        """Remove non-annotated ORFs sharing translated codons with mORFs."""
        retained: list[Candidate] = []
        overlap_cache: dict[
            tuple[
                bytes,
                bytes,
                bytes,
                tuple[tuple[int, int], ...],
            ],
            AnnotatedOverlap | None,
        ] = {}
        for candidate in candidates:
            if candidate.category in ANNOTATED_CATEGORIES:
                retained.append(candidate)
                continue
            cache_key = (
                candidate.gene_id,
                candidate.chrom,
                candidate.strand,
                candidate.blocks,
            )
            if cache_key in overlap_cache:
                match = overlap_cache[cache_key]
            else:
                match = self._best_annotated_overlap(candidate)
                overlap_cache[cache_key] = match
            if match is None:
                retained.append(candidate)
                continue

            fields = list(candidate.fields)
            fields[self.i_filter_status] = b"FAIL"
            fields[self.i_filter_reason] = match.relation
            fields[self.i_matched_annotated] = match.annotated_orf.annotated_orf_id
            fields[self.i_annotated_relation] = match.relation
            fields[self.i_annotated_overlap_nt] = str(match.overlap_nt).encode("ascii")
            fields[self.i_annotated_overlap_codon] = str(match.overlap_codon).encode("ascii")
            writer.write_removed(fields)

            reason = match.relation.decode("ascii")
            summary.annotated_overlap_removed += 1
            summary.removal_reasons[reason] += 1
            if candidate.category == b"lncORF" and match.relation == b"matches_annotated_ORF":
                summary.lnc_morf_removed += 1
        return retained
