#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-06-23
# Version: 0.2.7
# Function: Core functions for RiboParser P-site offset analysis.

import os.path
import sys
from collections import OrderedDict, defaultdict
from itertools import islice

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


class Mrna(object):
    def __init__(self, record):
        self.chromosome = record[0]
        self.gene_id = record[1]
        self.transcript_id = record[2]
        self.start = record[3]
        self.end = record[4]
        self.utr5_length = int(record[5])
        self.cds_length = int(record[6])
        self.utr3_length = int(record[7])
        self.strand = record[8]
        self.rep_transcript = record[9]
        self.modified = record[10]
        self.length = self.utr5_length + self.cds_length + self.utr3_length
        self.bam = []
        self.rpf = []
        self.seq = None


class Offset(object):
    def __init__(self, args):
        self.min_length = args.min
        self.max_length = args.max
        self.nt_num = self.max_length - self.min_length + 1
        self.mrna_file = args.transcript
        self.mode = getattr(args, "mode", "both")
        self.longest = getattr(args, "longest", False)

        self.tis_offset = {
            "tis_5end": OrderedDict(),
            "tis_3end": OrderedDict(),
            "tts_5end": OrderedDict(),
            "tts_3end": OrderedDict(),
        }
        self.adj_tis_offset = OrderedDict()
        self.tis_candidate_scores = OrderedDict()
        self.tis_5end = None
        self.tis_3end = None
        self.tts_5end = None
        self.tts_3end = None
        self.detail = args.detail

        self.frame_offset = OrderedDict()
        self.frame_offset_len = {}
        self.frame_candidate_counts = OrderedDict()
        self.merge_frame_offset = pd.DataFrame()

        self.exp_offset = None
        self.shift_nt = max(1, int(args.shift))
        self.mrna_dict = OrderedDict()
        self.length_dict = {}

        self.silence = args.silence
        self.sample_file = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]
        self.pysam_input = None
        self.output_prefix = args.output

        self.peak_length = args.exp_peak
        self.peak_reads = 0

        self.screen = getattr(args, "screen", False)
        self.min_mapq = getattr(args, "min_mapq", 10)
        self.min_offset_rpfs = getattr(args, "min_offset_rpfs", 50)
        self.dp_penalty = getattr(args, "dp_penalty", 0.18)
        self.read_filter_stats = defaultdict(int)

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------
    def _candidate_window(self, length):
        """Return a plausible 1-based P-site offset window for one read length."""
        shift_nt = (length - self.peak_length) // self.shift_nt
        left_offset = int(np.floor(self.peak_length / 3)) + shift_nt
        right_offset = int(np.ceil(self.peak_length * 3 / 4)) + shift_nt + 2
        left_offset = max(1, left_offset)
        right_offset = max(left_offset + 2, right_offset)
        return left_offset, right_offset, list(range(left_offset, right_offset + 1))

    @staticmethod
    def _safe_divide(num, den):
        if den == 0:
            return 0.0
        return float(num) / float(den)

    @staticmethod
    def _warning_join(warnings):
        warnings = [str(w) for w in warnings if w]
        return ";".join(warnings) if warnings else "PASS"

    def _expected_psite(self, length):
        shift_nt = (length - self.peak_length) // self.shift_nt
        return 12 + self.peak_length - 30 + shift_nt

    def _local_score(self, row, max_rpfs):
        """Score one candidate or final row on 0-1-ish scale."""
        confidence = float(row.get("confidence", 0.0)) / 100.0
        periodicity = float(row.get("periodicity", 0.0)) / 100.0
        rpfs = float(row.get("rpfs", 0.0))
        read_score = self._safe_divide(np.log1p(rpfs), np.log1p(max_rpfs)) if max_rpfs > 0 else 0.0
        return 0.52 * confidence + 0.33 * periodicity + 0.15 * read_score

    def _transition_penalty(self, prev_length, prev_psite, now_length, now_psite):
        """Penalty for biologically implausible jumps between adjacent read lengths."""
        expected_step = (now_length - prev_length) / float(self.shift_nt)
        observed_step = now_psite - prev_psite
        jump_penalty = abs(observed_step - expected_step)
        non_monotonic_penalty = 1.0 if observed_step < -1 else 0.0
        return self.dp_penalty * (jump_penalty + non_monotonic_penalty)

    def _dp_select_candidates(self, candidates_by_length):
        """Dynamic-programming selection of one offset per read length."""
        lengths = sorted(candidates_by_length)
        if not lengths:
            return {}

        max_rpfs = 0
        for length in lengths:
            for row in candidates_by_length[length]:
                max_rpfs = max(max_rpfs, float(row.get("rpfs", 0.0)))

        dp = {}
        back = {}
        for idx, length in enumerate(lengths):
            rows = candidates_by_length[length]
            dp[length] = {}
            back[length] = {}
            for cand_idx, row in enumerate(rows):
                psite = int(row["p_site"])
                local = self._local_score(row, max_rpfs)
                if idx == 0:
                    dp[length][cand_idx] = local
                    back[length][cand_idx] = None
                    continue

                prev_length = lengths[idx - 1]
                best_score = -np.inf
                best_prev_idx = None
                for prev_idx, prev_row in enumerate(candidates_by_length[prev_length]):
                    prev_psite = int(prev_row["p_site"])
                    score = (
                        dp[prev_length][prev_idx]
                        + local
                        - self._transition_penalty(prev_length, prev_psite, length, psite)
                    )
                    if score > best_score:
                        best_score = score
                        best_prev_idx = prev_idx
                dp[length][cand_idx] = best_score
                back[length][cand_idx] = best_prev_idx

        last_length = lengths[-1]
        last_idx = max(dp[last_length], key=dp[last_length].get)
        selected = {}
        for length in reversed(lengths):
            row = dict(candidates_by_length[length][last_idx])
            row["dp_score"] = round(float(dp[length][last_idx]), 4)
            selected[length] = row
            last_idx = back[length][last_idx]
            if last_idx is None:
                break
        return selected

    @staticmethod
    def offset_scale(offset):
        offset_norm = offset.sub(offset.mean(axis=1), axis=0)
        offset_norm_row = offset_norm.div(offset_norm.std(axis=1), axis=0)
        offset_norm_row.fillna(0, inplace=True)
        return offset_norm_row

    # ------------------------------------------------------------------
    # Input parsing
    # ------------------------------------------------------------------
    def read_transcript(self):
        """Read transcript annotation and optionally retain the longest CDS per gene."""
        selected = OrderedDict()
        with open(self.mrna_file, "r") as trans_file_in:
            for line in islice(trans_file_in, 1, None):
                record = line.strip().split("\t")
                if len(record) < 11:
                    if not self.silence:
                        print("Warning: malformed transcript line skipped: " + line.strip(), flush=True)
                    continue
                now_mrna = Mrna(record)

                if not self.longest:
                    if now_mrna.transcript_id in self.mrna_dict and not self.silence:
                        print(
                            "Warning: {gene} is duplicated in transcript file; keep the longer record.".format(
                                gene=now_mrna.transcript_id
                            ),
                            flush=True,
                        )
                    old = self.mrna_dict.get(now_mrna.transcript_id)
                    if old is None or now_mrna.length > old.length:
                        self.mrna_dict[now_mrna.transcript_id] = now_mrna
                    continue

                old = selected.get(now_mrna.gene_id)
                if old is None:
                    selected[now_mrna.gene_id] = now_mrna
                elif (now_mrna.cds_length, now_mrna.length) > (old.cds_length, old.length):
                    selected[now_mrna.gene_id] = now_mrna

        if self.longest:
            self.mrna_dict = OrderedDict(
                (mrna.transcript_id, mrna) for mrna in selected.values()
            )
        print("{number} transcripts retained for offset detection.".format(number=len(self.mrna_dict)), flush=True)

    def _read_passes_screen(self, read):
        if not self.screen:
            return True
        if read.is_unmapped:
            self.read_filter_stats["unmapped"] += 1
            return False
        if read.is_secondary:
            self.read_filter_stats["secondary"] += 1
            return False
        if read.is_supplementary:
            self.read_filter_stats["supplementary"] += 1
            return False
        if read.is_duplicate:
            self.read_filter_stats["duplicate"] += 1
            return False
        if read.mapping_quality < self.min_mapq:
            self.read_filter_stats["low_mapq"] += 1
            return False
        if read.infer_read_length() is None:
            self.read_filter_stats["no_length"] += 1
            return False
        if not read.get_blocks():
            self.read_filter_stats["no_blocks"] += 1
            return False
        return True

    def get_mrna_reads(self):
        """Read BAM/SAM alignments and attach reads to transcript objects."""
        if self.sample_format.lower() == ".bam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            read_format = "rb"
        elif self.sample_format.lower() == ".sam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            read_format = "r"
        else:
            print("Unknown file format; please input BAM or SAM.", flush=True)
            sys.exit()

        self.pysam_input = pysam.AlignmentFile(self.sample_file, read_format)
        for line in self.pysam_input.fetch(until_eof=True):
            if not self._read_passes_screen(line):
                continue
            if line.reference_name in self.mrna_dict:
                read_length = line.infer_read_length()
                if line.is_reverse:
                    self.length_dict.setdefault(read_length, [0, 0])[1] += 1
                else:
                    self.length_dict.setdefault(read_length, [0, 0])[0] += 1
                self.mrna_dict[line.reference_name].bam.append(line)
            elif self.silence:
                continue
            else:
                print(line.reference_name + ": not in reference file.", flush=True)

        if not self.length_dict:
            print("No reads retained for offset detection. Please check input BAM/SAM and filters.", flush=True)
            sys.exit()

        peak_length, counts = max(self.length_dict.items(), key=lambda item: sum(item[1]))
        self.peak_reads = sum(counts)
        if self.screen:
            print("Read filter summary: {stats}".format(stats=dict(self.read_filter_stats)), flush=True)
        print(
            "Empirical peak read length: {length} nt, reads: {reads}. Expected peak length used for offset window: {peak}.".format(
                length=peak_length, reads=self.peak_reads, peak=self.peak_length
            ),
            flush=True,
        )
        self.pysam_input.close()

    # ------------------------------------------------------------------
    # SSCBM: start/stop codon based model
    # ------------------------------------------------------------------
    def get_tis_offset(self):
        """Collect 5'/3' end profiles around TIS and TTS for each read length."""
        for number in range(self.min_length, self.max_length + 1):
            self.tis_offset["tis_5end"][number] = OrderedDict({-i: 0 for i in range(self.max_length, -1, -1)})
            self.tis_offset["tis_3end"][number] = OrderedDict({i: 0 for i in range(self.max_length + 1)})
            self.tis_offset["tts_5end"][number] = OrderedDict({-i: 0 for i in range(self.max_length, -1, -1)})
            self.tis_offset["tts_3end"][number] = OrderedDict({i: 0 for i in range(self.max_length + 1)})

        for _, mrna_attr in self.mrna_dict.items():
            if len(mrna_attr.bam) == 0:
                continue
            cds_start = mrna_attr.utr5_length
            cds_end = mrna_attr.utr5_length + mrna_attr.cds_length - 6

            for line in mrna_attr.bam:
                map_start, map_end = line.get_blocks()[0]
                read_length = line.infer_read_length()
                if not (self.min_length <= read_length <= self.max_length):
                    continue

                if map_start <= cds_start <= map_end:
                    offset_5end = map_start - cds_start
                    offset_3end = map_end - cds_start
                    if offset_5end in self.tis_offset["tis_5end"][read_length]:
                        self.tis_offset["tis_5end"][read_length][offset_5end] += 1
                    if offset_3end in self.tis_offset["tis_3end"][read_length]:
                        self.tis_offset["tis_3end"][read_length][offset_3end] += 1
                elif map_start <= cds_end <= map_end:
                    offset_5end = map_start - cds_end
                    offset_3end = map_end - cds_end
                    if offset_5end in self.tis_offset["tts_5end"][read_length]:
                        self.tis_offset["tts_5end"][read_length][offset_5end] += 1
                    if offset_3end in self.tis_offset["tts_3end"][read_length]:
                        self.tis_offset["tts_3end"][read_length][offset_3end] += 1

        self.tis_5end = pd.DataFrame(self.tis_offset["tis_5end"]).T
        self.tis_3end = pd.DataFrame(self.tis_offset["tis_3end"]).T
        self.tts_5end = pd.DataFrame(self.tis_offset["tts_5end"]).T
        self.tts_3end = pd.DataFrame(self.tis_offset["tts_3end"]).T

        if self.detail:
            self.tis_5end.to_csv(self.output_prefix + "_tis_5end.txt", sep="\t", header=True, index=True)
            self.tis_3end.to_csv(self.output_prefix + "_tis_3end.txt", sep="\t", header=True, index=True)
            self.tts_5end.to_csv(self.output_prefix + "_tts_5end.txt", sep="\t", header=True, index=True)
            self.tts_3end.to_csv(self.output_prefix + "_tts_3end.txt", sep="\t", header=True, index=True)

    def shift_codon(self, now_psite, length, left_offset, right_offset):
        """Shift a P-site by 3 nt until it falls in the expected range.

        Important fix: the previous implementation used self.adj_tis_offset[length - 1][-1],
        which is periodicity, not p_site. Here the true previous p_site is used.
        """
        if pd.isna(now_psite):
            return np.nan
        now_psite = int(now_psite)
        if left_offset <= now_psite <= right_offset:
            return now_psite

        previous = self.adj_tis_offset.get(length - 1)
        if previous:
            prev_psite = previous["p_site"] if isinstance(previous, dict) else previous[7]
            if not pd.isna(prev_psite):
                prev_psite = int(prev_psite)
                if prev_psite + 3 < now_psite:
                    return self.shift_codon(now_psite - 3, length, left_offset, right_offset)
                if now_psite < prev_psite - 3:
                    return self.shift_codon(now_psite + 3, length, left_offset, right_offset)
                return now_psite

        if now_psite > right_offset:
            return self.shift_codon(now_psite - 3, length, left_offset, right_offset)
        if now_psite < left_offset:
            return self.shift_codon(now_psite + 3, length, left_offset, right_offset)
        return np.nan

    def _profile_to_positive_offset(self, profile_df):
        profile = profile_df.T.copy()
        profile.index = abs(profile.index)
        return profile.groupby(profile.index).sum()

    def _score_sscbm_site(self, profile, length, max_offset_site, source):
        left_offset, right_offset, candidates = self._candidate_window(length)
        candidate_values = profile.reindex(candidates, fill_value=0)
        peak_rpfs = float(candidate_values.get(max_offset_site, 0.0))
        bg_values = candidate_values.drop(labels=[max_offset_site], errors="ignore")
        background = float(bg_values.mean()) if not bg_values.empty else 0.0
        peak_sharpness = self._safe_divide(peak_rpfs + 1.0, background + 1.0)

        frame0_range = range(max_offset_site - 3, max_offset_site + 6, 3)
        frame1_range = range(max_offset_site - 2, max_offset_site + 6, 3)
        frame2_range = range(max_offset_site - 1, max_offset_site + 6, 3)
        frame0 = float(profile.reindex(frame0_range, fill_value=0).sum())
        frame1 = float(profile.reindex(frame1_range, fill_value=0).sum())
        frame2 = float(profile.reindex(frame2_range, fill_value=0).sum())
        rpfs = frame0 + frame1 + frame2
        max_frame = max(frame0, frame1, frame2)
        periodicity = self._safe_divide(max_frame, rpfs)

        psite = max_offset_site + 1
        psite = self.shift_codon(psite, length, left_offset, right_offset + 1)

        read_score = min(1.0, self._safe_divide(rpfs, self.min_offset_rpfs))
        sharpness_score = min(1.0, np.log1p(peak_sharpness) / np.log(6.0))
        confidence = 100.0 * (0.5 * periodicity + 0.3 * sharpness_score + 0.20 * read_score)

        warnings = []
        if rpfs < self.min_offset_rpfs:
            warnings.append("low_count")
        if peak_sharpness < 1.5:
            warnings.append("flat_peak")
        if periodicity < 0.5:
            warnings.append("weak_periodicity")
        if pd.isna(psite):
            warnings.append("invalid_psite")

        return {
            "length": int(length),
            "source": source,
            "frame0": int(max_offset_site),
            "rpfs0": frame0,
            "frame1": int(max_offset_site + 1),
            "rpfs1": frame1,
            "frame2": int(max_offset_site + 2),
            "rpfs2": frame2,
            "p_site": int(psite) if not pd.isna(psite) else np.nan,
            "rpfs": rpfs,
            "periodicity": periodicity * 100.0,
            "peak_site": int(max_offset_site),
            "peak_rpfs": peak_rpfs,
            "background": background,
            "peak_sharpness": peak_sharpness,
            "confidence": confidence,
            "warning": self._warning_join(warnings),
        }

    def _best_sscbm_for_source(self, profile_df, length, source):
        profile_matrix = self._profile_to_positive_offset(profile_df)
        profile = profile_matrix[length].copy()
        left_offset, right_offset, candidates = self._candidate_window(length)
        candidate_values = profile.reindex(candidates, fill_value=0)

        rows = []
        if candidate_values.sum() == 0:
            seed = candidates[len(candidates) // 2]
            rows.append(self._score_sscbm_site(profile, length, seed, source))
            rows[-1]["warning"] = self._warning_join([rows[-1]["warning"], "no_signal"])
            return rows[-1], rows

        for candidate in candidates:
            rows.append(self._score_sscbm_site(profile, length, candidate, source))
        rows.sort(key=lambda row: (row["confidence"], row["peak_rpfs"], row["periodicity"]), reverse=True)
        return rows[0], rows

    def _consensus_sscbm(self, length, tis_row, tts_row):
        tis_valid = not pd.isna(tis_row["p_site"])
        tts_valid = not pd.isna(tts_row["p_site"])
        warnings = []

        if tis_valid and tts_valid:
            delta = abs(int(tis_row["p_site"]) - int(tts_row["p_site"]))
            if delta == 0:
                chosen = dict(tis_row if tis_row["confidence"] >= tts_row["confidence"] else tts_row)
                source = "TIS_TTS_consensus"
            elif delta <= 3 and delta % 3 == 0:
                chosen = dict(tis_row if tis_row["confidence"] >= tts_row["confidence"] else tts_row)
                source = "TIS_TTS_shifted_consensus"
                warnings.append("tis_tts_3nt_shift")
            else:
                if tis_row["confidence"] >= tts_row["confidence"]:
                    chosen = dict(tis_row)
                    source = "TIS"
                else:
                    chosen = dict(tts_row)
                    source = "TTS"
                warnings.append("tis_tts_conflict")
        elif tis_valid:
            chosen = dict(tis_row)
            source = "TIS_only"
            warnings.append("missing_tts")
        elif tts_valid:
            chosen = dict(tts_row)
            source = "TTS_only"
            warnings.append("missing_tis")
        else:
            chosen = dict(tis_row if tis_row["confidence"] >= tts_row["confidence"] else tts_row)
            source = "invalid"
            warnings.append("invalid_tis_tts")

        chosen["source"] = source
        chosen["warning"] = self._warning_join(
            [chosen.get("warning", "PASS"), tis_row.get("warning", "PASS"), tts_row.get("warning", "PASS")] + warnings
        ).replace("PASS;", "").replace(";PASS", "")
        chosen.update(
            {
                "tis_p_site": tis_row["p_site"],
                "tis_periodicity": tis_row["periodicity"],
                "tis_rpfs": tis_row["rpfs"],
                "tis_confidence": tis_row["confidence"],
                "tis_peak_sharpness": tis_row["peak_sharpness"],
                "tts_p_site": tts_row["p_site"],
                "tts_periodicity": tts_row["periodicity"],
                "tts_rpfs": tts_row["rpfs"],
                "tts_confidence": tts_row["confidence"],
                "tts_peak_sharpness": tts_row["peak_sharpness"],
            }
        )
        return chosen

    def adjust_tis_offset(self):
        """Detect SSCBM offset by independent TIS/TTS scoring plus consensus and DP smoothing."""
        candidates_by_length = OrderedDict()
        for length in range(self.min_length, self.max_length + 1):
            tis_best, tis_candidates = self._best_sscbm_for_source(self.tis_5end, length, "TIS")
            tts_best, tts_candidates = self._best_sscbm_for_source(self.tts_5end, length, "TTS")
            consensus = self._consensus_sscbm(length, tis_best, tts_best)

            # Candidate pool for DP: all TIS/TTS candidates, not a simple TIS+TTS sum.
            pool = []
            by_psite = {}
            for row in tis_candidates + tts_candidates + [consensus]:
                if pd.isna(row["p_site"]):
                    continue
                key = int(row["p_site"])
                if key not in by_psite or row["confidence"] > by_psite[key]["confidence"]:
                    by_psite[key] = dict(row)
            for row in by_psite.values():
                pool.append(row)
            candidates_by_length[length] = pool if pool else [consensus]

        selected = self._dp_select_candidates(candidates_by_length)
        for length in range(self.min_length, self.max_length + 1):
            row = selected.get(length)
            if row is None:
                row = candidates_by_length[length][0]
            row = dict(row)
            if "dp_smoothed" not in row.get("warning", ""):
                row["warning"] = self._warning_join([row.get("warning", "PASS"), "dp_smoothed"])
            self.adj_tis_offset[length] = row
        self.tis_candidate_scores = candidates_by_length

    def write_tis_offset(self):
        adj_tis_offset = pd.DataFrame(self.adj_tis_offset).T.copy()
        if adj_tis_offset.empty:
            print("No SSCBM offset result to write.", flush=True)
            return

        column_order = [
            "length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2",
            "p_site", "rpfs", "periodicity", "confidence", "warning", "source", "dp_score",
            "peak_site", "peak_rpfs", "background", "peak_sharpness",
            "tis_p_site", "tis_periodicity", "tis_rpfs", "tis_confidence", "tis_peak_sharpness",
            "tts_p_site", "tts_periodicity", "tts_rpfs", "tts_confidence", "tts_peak_sharpness",
        ]
        for col in column_order:
            if col not in adj_tis_offset.columns:
                adj_tis_offset[col] = np.nan
        adj_tis_offset = adj_tis_offset[column_order]

        percent_cols = ["periodicity", "confidence", "tis_periodicity", "tis_confidence", "tts_periodicity", "tts_confidence"]
        for col in percent_cols:
            adj_tis_offset[col] = pd.to_numeric(adj_tis_offset[col], errors="coerce").round(2)
        float_cols = ["background", "peak_sharpness", "tis_peak_sharpness", "tts_peak_sharpness", "dp_score"]
        for col in float_cols:
            adj_tis_offset[col] = pd.to_numeric(adj_tis_offset[col], errors="coerce").round(4)
        int_cols = ["length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2", "p_site", "rpfs"]
        for col in int_cols:
            adj_tis_offset[col] = pd.to_numeric(adj_tis_offset[col], errors="coerce").fillna(0).astype(int)

        adj_tis_offset.sort_values(["length"], inplace=True)
        adj_tis_offset["ribo"] = ["first"] * adj_tis_offset.shape[0]
        adj_tis_offset.to_csv(self.output_prefix + "_SSCBM_offset.txt", sep="\t", index=False)

    def draw_tis_heatmap(self):
        def draw_figure(tis_5end, tis_3end, tts_5end, tts_3end, out_pdf, out_png):
            matplotlib.use("Agg")
            now_cmap = "Blues"
            fig = plt.figure(figsize=(12, 8), dpi=300)
            ax1 = plt.subplot(2, 2, 1)
            sns.heatmap(data=tis_5end, annot=None, linewidths=0.5, ax=ax1, cmap=now_cmap)
            ax1.set_title("RPFs 5end")
            ax1.set_ylabel("RPFs length")
            ax1.set_xlabel("from start codon")
            ax2 = plt.subplot(2, 2, 2)
            sns.heatmap(data=tis_3end, annot=None, linewidths=0.5, ax=ax2, cmap=now_cmap)
            ax2.set_title("RPFs 3end")
            ax2.set_ylabel("RPFs length")
            ax2.set_xlabel("from start codon")
            ax3 = plt.subplot(2, 2, 3)
            sns.heatmap(data=tts_5end, annot=None, linewidths=0.5, ax=ax3, cmap=now_cmap)
            ax3.set_title("RPFs 5end")
            ax3.set_ylabel("RPFs length")
            ax3.set_xlabel("from stop codon")
            ax4 = plt.subplot(2, 2, 4)
            sns.heatmap(data=tts_3end, annot=None, linewidths=0.5, ax=ax4, cmap=now_cmap)
            ax4.set_title("RPFs 3end")
            ax4.set_ylabel("RPFs length")
            ax4.set_xlabel("from stop codon")
            plt.tight_layout()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close(fig)

        out_pdf = self.output_prefix + "_SSCBM_offset.pdf"
        out_png = self.output_prefix + "_SSCBM_offset.png"
        draw_figure(self.tis_5end, self.tis_3end, self.tts_5end, self.tts_3end, out_pdf, out_png)

        out_pdf_s = self.output_prefix + "_SSCBM_offset_scale.pdf"
        out_png_s = self.output_prefix + "_SSCBM_offset_scale.png"
        draw_figure(
            self.offset_scale(self.tis_5end),
            self.offset_scale(self.tis_3end),
            self.offset_scale(self.tts_5end),
            self.offset_scale(self.tts_3end),
            out_pdf_s,
            out_png_s,
        )

    # ------------------------------------------------------------------
    # RSBM: ribosome-structure/frame based model
    # ------------------------------------------------------------------
    def make_frame_offset(self):
        """Make candidate P-site offsets for RSBM using a window, not only 3 candidates."""
        self.exp_offset = 11 + self.peak_length - 30
        for length in range(self.min_length, self.max_length + 1):
            left_offset, right_offset, candidates = self._candidate_window(length)
            # If SSCBM is available, keep its p_site inside the candidate pool but do not rely on it exclusively.
            if length in self.adj_tis_offset:
                seed = int(self.adj_tis_offset[length]["p_site"])
                extra = list(range(max(left_offset, seed - 3), min(right_offset, seed + 3) + 1))
                candidates = sorted(set(candidates + extra))
            self.frame_offset_len[length] = candidates
            self.frame_candidate_counts[length] = OrderedDict(
                (psite, [0, 0, 0]) for psite in candidates
            )

    def get_mono_frame(self, mrna_attr, cds_start, cds_end):
        """Count frame distribution after candidate P-site projection.

        P-sites outside CDS are ignored, so UTR reads cannot drive RSBM offset selection.
        """
        for line in mrna_attr.bam:
            map_start, _ = line.get_blocks()[0]
            reads_length = line.infer_read_length()
            if not (self.min_length <= reads_length <= self.max_length):
                continue
            for psite in self.frame_offset_len.get(reads_length, []):
                psite_coord = map_start + int(psite) - 1
                if not (cds_start <= psite_coord <= cds_end):
                    continue
                frame = (psite_coord - cds_start) % 3
                self.frame_candidate_counts[reads_length][psite][frame] += 1

    def get_frame_offset(self):
        self.make_frame_offset()
        for _, mrna_attr in self.mrna_dict.items():
            if len(mrna_attr.bam) == 0:
                continue
            cds_start = mrna_attr.utr5_length
            cds_end = mrna_attr.utr5_length + mrna_attr.cds_length - 3
            self.get_mono_frame(mrna_attr, cds_start, cds_end)

    def _score_rsbm_candidate(self, length, psite, counts):
        frame0, frame1, frame2 = [float(i) for i in counts]
        rpfs = frame0 + frame1 + frame2
        periodicity = self._safe_divide(frame0, rpfs)
        read_score = min(1.0, self._safe_divide(rpfs, self.min_offset_rpfs))
        prior_distance = abs(int(psite) - self._expected_psite(length))
        prior_score = max(0.0, 1.0 - prior_distance / 10.0)
        confidence = 100.0 * (0.58 * periodicity + 0.27 * read_score + 0.15 * prior_score)

        warnings = []
        if rpfs < self.min_offset_rpfs:
            warnings.append("low_count")
        if periodicity < 0.5:
            warnings.append("weak_periodicity")
        if prior_distance > 6:
            warnings.append("far_from_structural_prior")

        return {
            "length": int(length),
            "frame0": int(psite),
            "rpfs0": frame0,
            "frame1": int(psite) + 1,
            "rpfs1": frame1,
            "frame2": int(psite) + 2,
            "rpfs2": frame2,
            "p_site": int(psite),
            "rpfs": rpfs,
            "periodicity": periodicity * 100.0,
            "confidence": confidence,
            "warning": self._warning_join(warnings),
            "source": "RSBM_window",
            "prior_psite": self._expected_psite(length),
        }

    def format_frame_offset(self):
        candidates_by_length = OrderedDict()
        for length in range(self.min_length, self.max_length + 1):
            rows = []
            for psite, counts in self.frame_candidate_counts.get(length, {}).items():
                rows.append(self._score_rsbm_candidate(length, psite, counts))
            if not rows:
                left_offset, right_offset, candidates = self._candidate_window(length)
                seed = candidates[len(candidates) // 2]
                rows = [self._score_rsbm_candidate(length, seed, [0, 0, 0])]
                rows[0]["warning"] = self._warning_join([rows[0]["warning"], "no_signal"])
            candidates_by_length[length] = rows

        selected = self._dp_select_candidates(candidates_by_length)
        self.merge_frame_offset = pd.DataFrame([selected[length] for length in sorted(selected)])
        if self.merge_frame_offset.empty:
            print("No RSBM offset result generated.", flush=True)
            return

        self.merge_frame_offset["dp_score"] = pd.to_numeric(self.merge_frame_offset["dp_score"], errors="coerce").round(4)
        self.merge_frame_offset["periodicity"] = pd.to_numeric(self.merge_frame_offset["periodicity"], errors="coerce").round(2)
        self.merge_frame_offset["confidence"] = pd.to_numeric(self.merge_frame_offset["confidence"], errors="coerce").round(2)
        self.merge_frame_offset["ribo"] = ["first"] * self.merge_frame_offset.shape[0]
        for col in ["length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2", "p_site", "rpfs"]:
            self.merge_frame_offset[col] = pd.to_numeric(self.merge_frame_offset[col], errors="coerce").fillna(0).astype(int)

    def adjust_frame_offset(self):
        """Compatibility wrapper.

        Offset continuity is now handled by dynamic programming in format_frame_offset().
        The method is kept so the CLI flow remains stable.
        """
        if self.merge_frame_offset.empty:
            return
        self.merge_frame_offset.sort_values(["length"], inplace=True)
        self.merge_frame_offset.reset_index(drop=True, inplace=True)

    def write_frame_offset(self):
        if self.merge_frame_offset.empty:
            print("No RSBM offset result to write.", flush=True)
            return
        column_order = [
            "length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2",
            "p_site", "rpfs", "periodicity", "confidence", "warning", "source", "prior_psite", "dp_score", "ribo",
        ]
        for col in column_order:
            if col not in self.merge_frame_offset.columns:
                self.merge_frame_offset[col] = np.nan
        self.merge_frame_offset[column_order].to_csv(
            self.output_prefix + "_RSBM_offset.txt", sep="\t", index=False
        )

    def draw_frame_heatmap(self):
        if self.merge_frame_offset.empty:
            return
        out_pdf = self.output_prefix + "_RSBM_offset.pdf"
        out_png = self.output_prefix + "_RSBM_offset.png"
        raw_frame_offset = self.merge_frame_offset.loc[:, ["length", "p_site", "rpfs0", "rpfs1", "rpfs2"]].copy()
        offset_num = raw_frame_offset.shape[0]
        raw_frame_offset.index = [
            str(raw_frame_offset["length"].to_list()[i]) + "_" + str(raw_frame_offset["p_site"].to_list()[i])
            for i in range(0, offset_num)
        ]
        raw_frame_offset = raw_frame_offset.drop(columns=["length", "p_site"])
        raw_frame_offset.columns = ["frame0", "frame1", "frame2"]
        raw_frame_offset = raw_frame_offset.apply(pd.to_numeric)
        scale_frame_offset = self.offset_scale(raw_frame_offset).apply(pd.to_numeric)

        matplotlib.use("Agg")
        now_cmap = "Blues"
        fig = plt.figure(figsize=(12, 6), dpi=300)
        ax1 = plt.subplot(1, 2, 1)
        sns.heatmap(data=raw_frame_offset, annot=None, linewidths=0.5, ax=ax1, cmap=now_cmap)
        ax1.set_title("raw counts")
        ax1.set_ylabel("RPFs length")
        ax1.set_xlabel("open reading frame")
        ax2 = plt.subplot(1, 2, 2)
        sns.heatmap(data=scale_frame_offset, annot=None, linewidths=0.5, ax=ax2, cmap=now_cmap)
        ax2.set_title("scaled counts")
        ax2.set_ylabel("RPFs length")
        ax2.set_xlabel("open reading frame")
        plt.tight_layout()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close(fig)
