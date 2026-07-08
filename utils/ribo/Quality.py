#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Rensc, modified by local fork
# Date: 2026-07-08
# Version: 0.2.8-dev.011
# Function: Provide optimized RiboParser core functions for rpf_Check quality analysis.
# Input: RiboParser transcript annotation and BAM/SAM alignment file.
# Output: Filtered sorted BAM, BAM index, length distribution, summary JSON, and optional saturation reports.

import glob
import json
import math
import os
import random
import shutil
import subprocess
from collections import Counter
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam
from PIL import Image, ImageDraw, ImageFont


class Quality(object):
    """Quality-control engine for RPF alignment checking.

    The default workflow writes the final filtered BAM required by downstream
    RiboParser modules. It automatically prepares an indexed BAM and uses ordered multi-process
    reference fetching whenever more than one worker is available. This avoids
    the slow single Python scanner for ordinary RPF transcriptome BAM files.
    """

    def __init__(self, args):
        # opts for mRNA file import
        self.longest = args.longest
        self.mrna_file = args.transcript
        self.mrna_ids = []
        self.mrna_set = set()
        self.mrna_dict = OrderedDict()
        self.gene_to_index = {}
        self.index_to_gene = []
        self.length_dict = {}

        # opts for BAM/SAM file import
        self.thread = max(2, int(args.thread))
        self.sort_memory = "1G"
        self.input_file = args.bam
        self.sample_file = args.bam
        self.sample_format = self._detect_alignment_mode(args.bam)
        self.output_bam = None
        self.split_bam_list = []
        self.tag = args.tag
        self.reverse = args.reverse
        self.align = args.align
        self.scan_mode = "auto"
        self.progress_every = 1000000
        self.temp_output_bam = None
        self.fetch_temp_bam = None

        # opts for optional alignment-flag retention
        self.keep_secondary = getattr(args, "secondary", False)
        self.keep_supplementary = getattr(args, "supplementary", False)
        self.keep_duplicate = getattr(args, "duplicate", False)

        # opts for BAM length summary
        self.peak_length = 0
        self.peak_reads = 0
        self.mono = [0.5, 100.5]
        self.profile = None

        # opts for RPF saturation
        self.saturation_flag = args.saturation
        self.saturation = []
        self.saturation_seed = 5201314
        self.x_ticks = np.array([i * 10 for i in range(1, 10)])
        self.saturation_gene_indices = []
        self.saturation_read_gene_sets = {}
        self.saturation_matrix = None

        # opts for file output
        self.output_prefix = args.output
        self.temp_output_bam = self._make_temp_bam(".filtered.tmp.bam")

        # opts for QC reports
        self.filter_stats = self._new_filter_stats()
        self.plot_warnings = []
        self.transcript_total_input = 0
        self.transcript_total_kept = 0
        self.transcript_removed_by_longest = 0
        self.bam_reference_match_count = 0
        self.fetch_zero_weight_references = 0

    @staticmethod
    def _detect_alignment_mode(alignment_file):
        """Return the pysam open mode inferred from file extension."""
        lower_file = alignment_file.lower()
        if lower_file.endswith(".bam"):
            return "rb"
        if lower_file.endswith(".sam"):
            return "r"
        if lower_file.endswith(".cram"):
            return "rc"
        raise RuntimeError(
            "Unknown alignment format. Please provide a BAM, SAM, or CRAM file."
        )

    def _alignment_open_kwargs(self):
        """Return safe pysam.AlignmentFile keyword arguments."""
        if self.sample_format in ("rb", "rc") and self.thread > 1:
            return {"threads": self.thread}
        return {}

    def read_transcript(self):
        """Read transcript annotation and optionally keep only longest transcripts.

        Expected RiboParser gene.norm.txt columns include:
        chromosome, gene_id, transcript_id, start, end, utr5_length, cds_length,
        utr3_length, strand, rep_transcript, modified.

        Workflow
        --------
        1. Read the transcript table with or without a header.
        2. Optionally keep one representative/longest transcript per gene.
        3. Build transcript lookup tables for fast BAM scanning and saturation.
        """
        trans_file = self.mrna_file

        try:
            trans_df = pd.read_csv(trans_file, sep="\t", dtype=str)
        except Exception as exc:
            raise RuntimeError(
                f"Failed to read transcript annotation file: {trans_file}. "
                f"Original error: {exc}"
            ) from exc

        if trans_df.empty:
            raise RuntimeError(f"Transcript annotation is empty: {trans_file}")

        if "transcript_id" not in trans_df.columns:
            trans_df = pd.read_csv(trans_file, sep="\t", header=None, dtype=str)
            if trans_df.shape[1] < 3:
                raise RuntimeError(
                    f"Transcript annotation must contain at least 3 columns; "
                    f"got {trans_df.shape[1]} columns in {trans_file}"
                )

            base_cols = [
                "chromosome", "gene_id", "transcript_id", "start", "end",
                "utr5_length", "cds_length", "utr3_length", "strand",
                "rep_transcript", "modified",
            ]
            trans_df.columns = base_cols[:trans_df.shape[1]]

        self.transcript_total_input = int(len(trans_df))

        if "gene_id" not in trans_df.columns:
            trans_df["gene_id"] = trans_df["transcript_id"]

        for col in ("utr5_length", "cds_length", "utr3_length"):
            if col in trans_df.columns:
                trans_df[col] = pd.to_numeric(
                    trans_df[col],
                    errors="coerce",
                ).fillna(0).astype(int)
            else:
                trans_df[col] = 0

        if self.longest:
            selected_df = self._select_longest_transcripts(trans_df)
            self.transcript_total_kept = int(len(selected_df))
            self.transcript_removed_by_longest = int(len(trans_df) - len(selected_df))
            self._write_longest_transcripts(selected_df, trans_df)
            trans_df = selected_df
        else:
            self.transcript_total_kept = int(len(trans_df))
            self.transcript_removed_by_longest = 0
        self.bam_reference_match_count = 0
        self.fetch_zero_weight_references = 0

        transcript_ids = []
        seen_ids = set()
        for transcript_id in trans_df["transcript_id"].astype(str):
            if not transcript_id or transcript_id.lower() == "nan":
                continue
            if transcript_id in seen_ids:
                continue
            transcript_ids.append(transcript_id)
            seen_ids.add(transcript_id)

        if not transcript_ids:
            raise RuntimeError(
                "No valid transcript_id was loaded from the annotation file. "
                "Please check whether the transcript file is the RiboParser gene.norm.txt format."
            )

        self.mrna_ids = transcript_ids
        self.mrna_set = set(transcript_ids)
        self.index_to_gene = list(transcript_ids)
        self.gene_to_index = {
            transcript_id: idx for idx, transcript_id in enumerate(transcript_ids)
        }
        self.mrna_dict = OrderedDict((transcript_id, [0] * 9) for transcript_id in transcript_ids)

    @staticmethod
    def _truthy_series(series):
        """Return a boolean mask for common truthy string values."""
        return series.astype(str).str.lower().isin(["true", "1", "yes", "y", "t"])

    def _select_longest_transcripts(self, trans_df):
        """Select one representative/longest transcript per gene."""
        df = trans_df.copy()

        if "rep_transcript" in df.columns:
            rep_mask = self._truthy_series(df["rep_transcript"])
            rep_df = df.loc[rep_mask].copy()
            if not rep_df.empty and rep_df["gene_id"].nunique() == len(rep_df):
                return rep_df

        df["gene_length"] = df["utr5_length"] + df["cds_length"] + df["utr3_length"]
        df = df.sort_values(
            by=[
                "gene_id", "cds_length", "gene_length",
                "utr5_length", "utr3_length", "transcript_id",
            ],
            ascending=[True, False, False, False, False, True],
        )
        return df.groupby("gene_id", as_index=False, sort=False).head(1)

    def _write_longest_transcripts(self, selected_df, all_df):
        """Write a transcript-retention report when longest mode is enabled."""
        out_file = self.output_prefix + "_longest_transcripts.txt"
        selected_ids = set(selected_df["transcript_id"].astype(str))

        report_cols = [
            col for col in [
                "gene_id", "transcript_id", "cds_length", "utr5_length",
                "utr3_length", "rep_transcript",
            ] if col in all_df.columns
        ]

        report = all_df.loc[:, report_cols].copy()
        report["selected_by_rpf_check"] = report["transcript_id"].astype(str).isin(selected_ids)
        report.to_csv(out_file, sep="\t", index=False)

    @staticmethod
    def _require_samtools():
        """Return the external samtools executable path or fail with a clear error."""
        samtools_path = shutil.which("samtools")
        if not samtools_path:
            raise RuntimeError(
                "samtools was not found in PATH. Please install samtools and make sure it is available."
            )
        return samtools_path

    @staticmethod
    def _run_samtools_command(command, step_name):
        """Run a samtools command and keep successful runs quiet.

        Samtools linked from conda/micromamba can emit repeated ncurses/libtinfo
        compatibility messages even when the command succeeds. Those lines do not
        help QC interpretation, so successful commands stay silent. Full stdout
        and stderr are still included when samtools exits with an error.
        """
        process = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if process.returncode != 0:
            raise RuntimeError(
                f"{step_name} failed with return code {process.returncode}\n"
                f"Command: {' '.join(command)}\n"
                f"STDOUT:\n{process.stdout}\n"
                f"STDERR:\n{process.stderr}"
            )

        return process

    def _make_temp_prefix(self, label):
        """Return a samtools temporary prefix next to the output prefix."""
        return self.output_prefix + "." + label

    def _make_temp_bam(self, suffix):
        """Return a temporary BAM path next to the output prefix."""
        return self.output_prefix + suffix

    def _samtools_sort(self, input_bam, output_bam, temp_prefix=None, threads=None, memory=None):
        """Sort BAM/SAM with samtools CLI first and pysam as a fallback."""
        sort_thread = max(1, int(threads if threads is not None else self.thread))
        sort_memory = memory if memory is not None else self.sort_memory
        temp_prefix = temp_prefix if temp_prefix is not None else output_bam + ".tmp"

        try:
            command = [
                self._require_samtools(),
                "sort",
                "-@",
                str(sort_thread),
                "-m",
                str(sort_memory),
                "-T",
                temp_prefix,
                "-o",
                output_bam,
                input_bam,
            ]
            self._run_samtools_command(command, "samtools sort")
        except Exception as samtools_error:
            print(
                f"samtools sort failed; fallback to pysam.sort. Error: {samtools_error}",
                flush=True,
            )
            try:
                pysam.sort(
                    "-@",
                    str(sort_thread),
                    "-m",
                    str(sort_memory),
                    "-T",
                    temp_prefix,
                    "-o",
                    output_bam,
                    input_bam,
                )
            except Exception as pysam_error:
                raise RuntimeError(
                    "Both samtools sort and pysam.sort failed.\n"
                    f"samtools error: {samtools_error}\n"
                    f"pysam error: {pysam_error}"
                ) from pysam_error

        self._assert_non_empty_file(output_bam, "Sorted BAM")

    def _samtools_index(self, bam_file, threads=None):
        """Index BAM with samtools CLI first and pysam as a fallback."""
        index_thread = max(1, int(threads if threads is not None else self.thread))
        try:
            command = [
                self._require_samtools(),
                "index",
                "-@",
                str(index_thread),
                bam_file,
            ]
            self._run_samtools_command(command, "samtools index")
        except Exception as samtools_error:
            print(
                f"samtools index failed; fallback to pysam.index. Error: {samtools_error}",
                flush=True,
            )
            try:
                if index_thread > 1:
                    pysam.index("-@", str(index_thread), bam_file)
                else:
                    pysam.index(bam_file)
            except Exception as pysam_error:
                raise RuntimeError(
                    "Both samtools index and pysam.index failed.\n"
                    f"samtools error: {samtools_error}\n"
                    f"pysam error: {pysam_error}"
                ) from pysam_error

    def _samtools_cat(self, output_bam, bam_files):
        """Concatenate BAM files with samtools CLI first and pysam as a fallback."""
        bam_files = list(bam_files)

        if not bam_files:
            raise RuntimeError("No BAM files were provided to samtools cat.")

        if len(bam_files) == 1:
            shutil.copyfile(bam_files[0], output_bam)
            self._assert_non_empty_file(output_bam, "Concatenated BAM")
            return

        try:
            command = [
                self._require_samtools(),
                "cat",
                "-o",
                output_bam,
                *bam_files,
            ]
            self._run_samtools_command(command, "samtools cat")
        except Exception as samtools_error:
            print(
                f"samtools cat failed; fallback to pysam.cat. Error: {samtools_error}",
                flush=True,
            )
            try:
                pysam.cat("-o", output_bam, *bam_files)
            except Exception as pysam_error:
                raise RuntimeError(
                    "Both samtools cat and pysam.cat failed.\n"
                    f"samtools error: {samtools_error}\n"
                    f"pysam error: {pysam_error}"
                ) from pysam_error

        self._assert_non_empty_file(output_bam, "Concatenated BAM")

    def _samtools_merge(self, output_bam, bam_files, threads=None):
        """Merge coordinate-sorted split BAM files into one final BAM.

        This is faster than cat -> sort when each split BAM is already sorted.
        The split BAMs generated by fetch mode are written in BAM-header
        reference order, so they are valid inputs for samtools merge.
        """
        bam_files = list(bam_files)
        merge_thread = max(1, int(threads if threads is not None else self.thread))

        if not bam_files:
            raise RuntimeError("No BAM files were provided to samtools merge.")

        if len(bam_files) == 1:
            shutil.copyfile(bam_files[0], output_bam)
            self._assert_non_empty_file(output_bam, "Merged BAM")
            return

        try:
            command = [
                self._require_samtools(),
                "merge",
                "-@",
                str(merge_thread),
                "-f",
                output_bam,
                *bam_files,
            ]
            self._run_samtools_command(command, "samtools merge")
        except Exception as samtools_error:
            print(
                f"samtools merge failed; fallback to pysam.merge. Error: {samtools_error}",
                flush=True,
            )
            try:
                pysam.merge("-@", str(merge_thread), "-f", output_bam, *bam_files)
            except Exception as pysam_error:
                raise RuntimeError(
                    "Both samtools merge and pysam.merge failed.\n"
                    f"samtools error: {samtools_error}\n"
                    f"pysam error: {pysam_error}"
                ) from pysam_error

        self._assert_non_empty_file(output_bam, "Merged BAM")

    @staticmethod
    def _assert_non_empty_file(file_name, label):
        """Fail if a required file does not exist or is empty."""
        if not os.path.exists(file_name):
            raise RuntimeError(f"{label} was not generated: {file_name}")
        if os.path.getsize(file_name) == 0:
            raise RuntimeError(f"{label} is empty: {file_name}")

    @staticmethod
    def _has_bam_index(bam_file):
        """Return True when a BAM index exists next to the BAM file."""
        return (
            os.path.exists(bam_file + ".bai")
            or os.path.exists(os.path.splitext(bam_file)[0] + ".bai")
            or os.path.exists(bam_file + ".csi")
        )

    @staticmethod
    def _bam_sort_order(bam_file):
        """Return BAM header sort order when it is available."""
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam_in:
                header_dict = bam_in.header.to_dict()
            return str(header_dict.get("HD", {}).get("SO", "unknown")).lower()
        except Exception:
            return "unknown"

    def sort_index_bam(self, input_bam=None, output_bam=None, cleanup=False, try_index_first=False):
        """Finalize a BAM file as sorted and indexed output.

        When ``try_index_first`` is True, the method first assumes the input BAM
        is already coordinate-sorted. If indexing succeeds, it skips the expensive
        final sort step. If indexing fails, it falls back to samtools sort/index.
        """
        input_bam = input_bam or self.sample_file
        output_bam = output_bam or self.output_prefix + ".bam"
        sorted_tmp = self._make_temp_bam(".sorting.tmp.bam")
        temp_prefix = self._make_temp_prefix("sort.tmp")

        stale_files = [
            sorted_tmp,
            output_bam + ".bai",
            output_bam + ".csi",
            os.path.splitext(output_bam)[0] + ".bai",
            os.path.splitext(output_bam)[0] + ".csi",
        ]
        if os.path.abspath(input_bam) != os.path.abspath(output_bam):
            stale_files.append(output_bam)

        for stale_file in stale_files:
            if os.path.exists(stale_file):
                os.remove(stale_file)

        if try_index_first:
            print(f"try direct final BAM indexing: {input_bam} -> {output_bam}", flush=True)
            if os.path.abspath(input_bam) != os.path.abspath(output_bam):
                os.replace(input_bam, output_bam)
            try:
                self._samtools_index(output_bam, threads=self.thread)
                self.output_bam = output_bam
                print("direct final BAM indexing succeeded; skip final sort", flush=True)
                return output_bam
            except Exception as index_error:
                print(
                    "direct final BAM indexing failed; fallback to final sort. "
                    f"Error: {index_error}",
                    flush=True,
                )
                for stale_index in [
                    output_bam + ".bai",
                    output_bam + ".csi",
                    os.path.splitext(output_bam)[0] + ".bai",
                    os.path.splitext(output_bam)[0] + ".csi",
                ]:
                    if os.path.exists(stale_index):
                        os.remove(stale_index)
                input_bam = output_bam

        print(f"sort final BAM: {input_bam} -> {output_bam}", flush=True)
        self._samtools_sort(
            input_bam,
            sorted_tmp,
            temp_prefix=temp_prefix,
            threads=self.thread,
            memory=self.sort_memory,
        )
        os.replace(sorted_tmp, output_bam)

        print(f"index final BAM: {output_bam}", flush=True)
        self._samtools_index(output_bam, threads=self.thread)
        self.output_bam = output_bam

        if cleanup and os.path.abspath(input_bam) != os.path.abspath(output_bam):
            if os.path.exists(input_bam):
                os.remove(input_bam)

        for stale_file in glob.glob(temp_prefix + "*"):
            if os.path.exists(stale_file):
                os.remove(stale_file)

        return output_bam

    def _prepare_indexed_bam_for_fetch(self):
        """Ensure an indexed BAM is available for fetch-mode random access."""
        if self.sample_format == "rb" and self._has_bam_index(self.sample_file):
            self.fetch_temp_bam = None
            print(f"reuse existing BAM index for fetch mode: {self.sample_file}", flush=True)
            return

        sort_order = self._bam_sort_order(self.sample_file)
        if self.sample_format == "rb" and sort_order == "coordinate":
            self.fetch_temp_bam = None
            print(f"index coordinate-sorted BAM for fetch mode: {self.sample_file}", flush=True)
            self._samtools_index(self.sample_file, threads=self.thread)
            return

        temp_sorted_bam = self._make_temp_bam(".fetch.sorted.bam")
        temp_sort_prefix = self._make_temp_prefix("fetch.sort.tmp")
        self.fetch_temp_bam = temp_sorted_bam

        for stale_file in [temp_sorted_bam, temp_sorted_bam + ".bai", temp_sorted_bam + ".csi"]:
            if os.path.exists(stale_file):
                os.remove(stale_file)

        print(f"sort input BAM for fetch mode: {self.sample_file} -> {temp_sorted_bam}", flush=True)
        self._samtools_sort(
            self.sample_file,
            temp_sorted_bam,
            temp_prefix=temp_sort_prefix,
            threads=self.thread,
            memory=self.sort_memory,
        )
        print(f"index sorted fetch BAM: {temp_sorted_bam}", flush=True)
        self._samtools_index(temp_sorted_bam, threads=self.thread)
        self.sample_file = temp_sorted_bam
        self.sample_format = "rb"

    @staticmethod
    def _new_filter_stats():
        """Return an initialized filter-statistics counter."""
        return Counter({
            "total_records_scanned": 0,
            "non_mrna_reads": 0,
            "unmapped_reads": 0,
            "secondary_reads": 0,
            "supplementary_reads": 0,
            "duplicate_reads": 0,
            "qc_fail_reads": 0,
            "removed_secondary_reads": 0,
            "removed_supplementary_reads": 0,
            "removed_duplicate_reads": 0,
            "missing_length_reads": 0,
            "missing_unique_tag_reads": 0,
            "removed_multi_mapping_reads": 0,
            "removed_reverse_strand_reads": 0,
            "kept_plus_reads": 0,
            "kept_minus_reads": 0,
            "passed_filter_reads": 0,
            "written_reads": 0,
            "fetch_missing_reference": 0,
        })

    @staticmethod
    def _passes_unique_filter(reads, align, tag, stats):
        """Return True when an alignment passes unique-read filtering."""
        if tag == 0:
            return True

        align = str(align).lower()

        if reads.has_tag("NH"):
            try:
                return reads.get_tag("NH") <= 1
            except Exception:
                stats["missing_unique_tag_reads"] += 1
                return False

        if align in ("bowtie2", "bowtie"):
            if reads.has_tag("XS"):
                return False

        if reads.mapping_quality is not None and reads.mapping_quality > 0:
            return True

        stats["missing_unique_tag_reads"] += 1
        return False

    @staticmethod
    def _record_length(length_dict, read_length, is_reverse):
        """Record read length by strand."""
        if read_length not in length_dict:
            length_dict[read_length] = [0, 0]
        if is_reverse:
            length_dict[read_length][1] += 1
        else:
            length_dict[read_length][0] += 1

    @staticmethod
    def _passes_flag_filters(reads, stats, keep_secondary, keep_supplementary, keep_duplicate):
        """Return True when an alignment passes SAM flag-based filtering."""
        if reads.is_unmapped:
            stats["unmapped_reads"] += 1
            return False

        if reads.is_secondary:
            stats["secondary_reads"] += 1
            if not keep_secondary:
                stats["removed_secondary_reads"] += 1
                return False

        if reads.is_supplementary:
            stats["supplementary_reads"] += 1
            if not keep_supplementary:
                stats["removed_supplementary_reads"] += 1
                return False

        if reads.is_duplicate:
            stats["duplicate_reads"] += 1
            if not keep_duplicate:
                stats["removed_duplicate_reads"] += 1
                return False

        if reads.is_qcfail:
            stats["qc_fail_reads"] += 1

        return True

    @staticmethod
    def _update_saturation_cache(
        reads,
        gene_index,
        tag,
        saturation_gene_indices,
        saturation_read_gene_sets,
    ):
        """Update saturation cache using integer gene indices."""
        if tag == 1:
            saturation_gene_indices.append(gene_index)
            return

        try:
            saturation_read_gene_sets[reads.query_name].add(gene_index)
        except KeyError:
            saturation_read_gene_sets[reads.query_name] = {gene_index}

    @staticmethod
    def _process_alignment_record(
        reads,
        reference_name,
        length_dict,
        stats,
        gene_to_index,
        tag,
        align,
        saturation_flag,
        saturation_gene_indices,
        saturation_read_gene_sets,
        keep_secondary=True,
        keep_supplementary=True,
        keep_duplicate=True,
        bam_out=None,
    ):
        """Apply all QC filters to one alignment record and optionally write it."""
        gene_index = gene_to_index.get(reference_name)
        if gene_index is None:
            stats["non_mrna_reads"] += 1
            return

        if not Quality._passes_flag_filters(
            reads,
            stats,
            keep_secondary,
            keep_supplementary,
            keep_duplicate,
        ):
            return

        if not Quality._passes_unique_filter(reads, align, tag, stats):
            stats["removed_multi_mapping_reads"] += 1
            return

        # query_length is much cheaper than infer_read_length() for ordinary RPF reads.
        # Use infer_read_length() only as a fallback for unusual records.
        read_length = reads.query_length
        if read_length is None:
            read_length = reads.infer_read_length()
        if read_length is None:
            stats["missing_length_reads"] += 1
            return

        # Count plus and minus strand reads independently after all active filters.
        Quality._record_length(length_dict, read_length, reads.is_reverse)

        if reads.is_reverse:
            stats["kept_minus_reads"] += 1
        else:
            stats["kept_plus_reads"] += 1

        stats["passed_filter_reads"] += 1

        if bam_out is not None:
            bam_out.write(reads)
            stats["written_reads"] += 1

        if saturation_flag:
            Quality._update_saturation_cache(
                reads,
                gene_index,
                tag,
                saturation_gene_indices,
                saturation_read_gene_sets,
            )

    def _merge_length_dict(self, length_dict):
        """Merge one length dictionary into the object-level length dictionary."""
        for key, value in length_dict.items():
            if key in self.length_dict:
                self.length_dict[key][0] += value[0]
                self.length_dict[key][1] += value[1]
            else:
                self.length_dict[key] = list(value)

    def _merge_filter_stats(self, stats):
        """Merge one stats counter into the object-level stats counter."""
        for key, value in stats.items():
            self.filter_stats[key] = self.filter_stats.get(key, 0) + value

    def _merge_saturation_cache(self, gene_indices, read_gene_sets):
        """Merge worker-level saturation caches into object-level caches."""
        if gene_indices:
            self.saturation_gene_indices.extend(gene_indices)

        if read_gene_sets:
            for read_name, gene_set in read_gene_sets.items():
                self.saturation_read_gene_sets.setdefault(read_name, set()).update(gene_set)

    def _choose_scan_mode(self):
        """Choose a fast scan mode without falling back to a single Python scanner.

        The previous auto mode could fall back to a single Python scanner for unsorted BAM/SAM files.
        That was safe but slow for large transcriptome alignments. The optimized
        policy is now aggressive: the input alignment is automatically converted,
        sorted, and indexed when needed, then scanned by multi-process
        transcript-reference fetching.
        """
        if self.thread <= 1:
            self.thread = 2
            print(
                "thread <= 1 was requested; rpf_Check uses 2 workers to avoid single-process filtering",
                flush=True,
            )

        try:
            with pysam.AlignmentFile(
                self.sample_file,
                self.sample_format,
                **self._alignment_open_kwargs(),
            ) as bam_in:
                references = set(bam_in.references)
        except Exception as error:
            raise RuntimeError(
                "Cannot read alignment header. rpf_Check requires reference names "
                "to match transcript_id values from the annotation file. "
                f"Original error: {error}"
            ) from error

        matched_reference_count = sum(
            1 for transcript_id in self.mrna_ids if transcript_id in references
        )
        self.bam_reference_match_count = int(matched_reference_count)
        print(
            "BAM/transcript reference match: {matched}/{total}".format(
                matched=matched_reference_count,
                total=len(self.mrna_ids),
            ),
            flush=True,
        )

        if matched_reference_count < 1:
            raise RuntimeError(
                "None of the transcript_id values from the annotation file matched "
                "the alignment reference names. rpf_Check expects transcriptome BAM/SAM "
                "references generated by rpf_Reference."
            )

        return "fetch"

    def _cleanup_final_bam_intermediates(self):
        """Remove temporary BAM and samtools sort files after successful finalization."""
        stale_files = [
            self.temp_output_bam,
            self.temp_output_bam + ".bai",
            self.temp_output_bam + ".csi",
            self.fetch_temp_bam,
            self.fetch_temp_bam + ".bai" if self.fetch_temp_bam else None,
            self.fetch_temp_bam + ".csi" if self.fetch_temp_bam else None,
        ]

        for stale_file in stale_files:
            if stale_file and os.path.exists(stale_file):
                os.remove(stale_file)

        for pattern in [
            self.output_prefix + ".bam.sort.tmp*",
            self.output_prefix + ".filtered.tmp.bam.sort.tmp*",
            self.output_prefix + ".fetch.sort.tmp*",
            self.output_prefix + ".merged.unsorted.bam.sort.tmp*",
        ]:
            for stale_file in glob.glob(pattern):
                if os.path.exists(stale_file):
                    os.remove(stale_file)

    def scan_mrna_reads(self, mode=None):
        """Scan mRNA-aligned reads with multi-process fetch mode and write final BAM."""
        if mode is not None and str(mode).lower() not in {"auto", "fetch"}:
            raise RuntimeError("scan_mrna_reads() only supports automatic multi-process fetch mode.")

        resolved_mode = self._choose_scan_mode()
        self.scan_mode = resolved_mode

        self.scan_mrna_reads_fetch()
        self.ensure_non_empty_results()
        self.merge_sort_index_bam()

        self._cleanup_final_bam_intermediates()
        print("rpf_Check parallel scan and final BAM generation done", flush=True)

    @staticmethod
    def _contiguous_split_references(reference_weights, worker_count):
        """Split references into ordered contiguous chunks by mapped-read weight.

        The previous bin-packing splitter produced well-balanced workers, but the
        reference sets were interleaved across workers. That forced a slow
        coordinate merge after all ``*_split_*.bam`` files were generated. This
        splitter keeps each worker as a contiguous BAM-header reference interval,
        so the final BAM can be produced by fast concatenation instead of a full
        coordinate merge.
        """
        worker_count = max(1, int(worker_count))
        ordered_refs = sorted(reference_weights, key=lambda item: int(item[2]))
        if worker_count >= len(ordered_refs):
            return [[(item[0], int(item[3]))] for item in ordered_refs], [int(item[1]) for item in ordered_refs]

        total_weight = sum(max(1, int(item[1])) for item in ordered_refs)
        target_weight = max(1, int(np.ceil(total_weight / worker_count)))
        splits = []
        split_weights = []
        current_refs = []
        current_weight = 0
        remaining_workers = worker_count

        for idx, (reference_name, weight, _reference_order, gene_index) in enumerate(ordered_refs):
            remaining_refs = len(ordered_refs) - idx
            weight = max(1, int(weight))
            should_close = (
                current_refs
                and current_weight >= target_weight
                and remaining_workers > 1
                and remaining_refs >= remaining_workers
            )
            if should_close:
                splits.append(current_refs)
                split_weights.append(current_weight)
                current_refs = []
                current_weight = 0
                remaining_workers -= 1

            current_refs.append((reference_name, int(gene_index)))
            current_weight += weight

        if current_refs:
            splits.append(current_refs)
            split_weights.append(current_weight)

        return splits, split_weights

    def _get_reference_weights_for_fetch(self):
        """Return mRNA reference weights from BAM index statistics.

        If an index statistic is unavailable or reports zero mapped reads for a
        valid transcript reference, the weight falls back to 1. This keeps the program in multi-process fetch mode instead of
        silently reverting to a single-process scanner.
        """
        with pysam.AlignmentFile(self.sample_file, self.sample_format) as bam_in:
            references = list(bam_in.references)
            reference_order = {name: idx for idx, name in enumerate(references)}
            try:
                index_stats = bam_in.get_index_statistics()
            except Exception:
                index_stats = []

        reference_set = set(references)
        mapped_by_reference = {
            item.contig: int(item.mapped) for item in index_stats
        }

        reference_weights = []
        missing_reference_count = 0
        zero_weight_count = 0
        for transcript_id in self.mrna_ids:
            if transcript_id not in reference_set:
                missing_reference_count += 1
                continue

            weight = int(mapped_by_reference.get(transcript_id, 0))
            if weight <= 0:
                zero_weight_count += 1
                weight = 1

            reference_weights.append((
                transcript_id,
                weight,
                reference_order.get(transcript_id, 0),
                self.gene_to_index[transcript_id],
            ))

        self.filter_stats["fetch_missing_reference"] += missing_reference_count

        if not reference_weights:
            raise RuntimeError(
                "None of the transcript_id values from the annotation file were found in the BAM header."
            )

        self.fetch_zero_weight_references = int(zero_weight_count)

        return reference_weights

    @staticmethod
    def _scan_fetch_worker(args):
        """Worker for ordered fetch-mode scanning and split-BAM writing.

        The hot loop is intentionally inlined to reduce per-read Python function
        calls. Each worker receives transcript references as ``(name, gene_index)``
        pairs, so no dictionary lookup is needed for each alignment record.
        """
        (
            bam_in_file,
            sample_format,
            bam_out_file,
            reference_list,
            tag,
            align,
            saturation_flag,
            keep_secondary,
            keep_supplementary,
            keep_duplicate,
        ) = args

        length_dict = {}
        stats = Quality._new_filter_stats()
        saturation_gene_indices = []
        saturation_read_gene_sets = {}
        align = str(align).lower()
        use_unique_filter = int(tag) == 1

        for stale_file in [bam_out_file, bam_out_file + ".bai", bam_out_file + ".csi"]:
            if os.path.exists(stale_file):
                os.remove(stale_file)

        with pysam.AlignmentFile(bam_in_file, sample_format) as bam_in:
            with pysam.AlignmentFile(bam_out_file, "wb", template=bam_in) as bam_out:
                write_read = bam_out.write

                for reference_name, gene_index in reference_list:
                    try:
                        read_iter = bam_in.fetch(reference_name)
                    except ValueError:
                        stats["fetch_missing_reference"] += 1
                        continue

                    for reads in read_iter:
                        stats["total_records_scanned"] += 1

                        if reads.is_unmapped:
                            stats["unmapped_reads"] += 1
                            continue

                        if reads.is_secondary:
                            stats["secondary_reads"] += 1
                            if not keep_secondary:
                                stats["removed_secondary_reads"] += 1
                                continue

                        if reads.is_supplementary:
                            stats["supplementary_reads"] += 1
                            if not keep_supplementary:
                                stats["removed_supplementary_reads"] += 1
                                continue

                        if reads.is_duplicate:
                            stats["duplicate_reads"] += 1
                            if not keep_duplicate:
                                stats["removed_duplicate_reads"] += 1
                                continue

                        if reads.is_qcfail:
                            stats["qc_fail_reads"] += 1

                        if use_unique_filter:
                            keep_read = True
                            if reads.has_tag("NH"):
                                try:
                                    keep_read = reads.get_tag("NH") <= 1
                                except Exception:
                                    stats["missing_unique_tag_reads"] += 1
                                    keep_read = False
                            elif align in ("bowtie2", "bowtie") and reads.has_tag("XS"):
                                keep_read = False
                            elif reads.mapping_quality is not None and reads.mapping_quality > 0:
                                keep_read = True
                            else:
                                stats["missing_unique_tag_reads"] += 1
                                keep_read = False

                            if not keep_read:
                                stats["removed_multi_mapping_reads"] += 1
                                continue

                        read_length = reads.query_length
                        if read_length is None:
                            read_length = reads.infer_read_length()
                        if read_length is None:
                            stats["missing_length_reads"] += 1
                            continue

                        length_counts = length_dict.get(read_length)
                        if length_counts is None:
                            length_counts = [0, 0]
                            length_dict[read_length] = length_counts

                        if reads.is_reverse:
                            length_counts[1] += 1
                            stats["kept_minus_reads"] += 1
                        else:
                            length_counts[0] += 1
                            stats["kept_plus_reads"] += 1

                        stats["passed_filter_reads"] += 1
                        write_read(reads)
                        stats["written_reads"] += 1

                        if saturation_flag:
                            if use_unique_filter:
                                saturation_gene_indices.append(gene_index)
                            else:
                                try:
                                    saturation_read_gene_sets[reads.query_name].add(gene_index)
                                except KeyError:
                                    saturation_read_gene_sets[reads.query_name] = {gene_index}

        return length_dict, saturation_gene_indices, saturation_read_gene_sets, stats, bam_out_file

    def scan_mrna_reads_fetch(self):
        """Scan reads by indexed transcript references with ordered contiguous workloads."""
        if not self.mrna_set:
            raise RuntimeError("No mRNA records were loaded. Cannot scan reads.")

        print("ordered fetch scanner start", flush=True)
        self.filter_stats = self._new_filter_stats()
        self.length_dict = {}
        self.saturation_gene_indices = []
        self.saturation_read_gene_sets = {}

        self._prepare_indexed_bam_for_fetch()
        reference_weights = self._get_reference_weights_for_fetch()
        worker_count = min(max(1, int(self.thread)), len(reference_weights))
        reference_splits, worker_weights = self._contiguous_split_references(reference_weights, worker_count)
        print(
            "parallel fetch workers: {workers}; transcript references: {refs}".format(
                workers=len(reference_splits),
                refs=len(reference_weights),
            ),
            flush=True,
        )

        self.split_bam_list = [
            self.output_prefix + "_split_" + str(i) + ".bam"
            for i in range(len(reference_splits))
        ]

        args = [
            (
                self.sample_file,
                self.sample_format,
                self.split_bam_list[i],
                reference_list,
                self.tag,
                self.align,
                self.saturation_flag,
                self.keep_secondary,
                self.keep_supplementary,
                self.keep_duplicate,
            )
            for i, reference_list in enumerate(reference_splits)
        ]

        if len(args) == 1:
            results = [self._scan_fetch_worker(args[0])]
        else:
            with Pool(processes=len(args)) as pool:
                results = pool.map(self._scan_fetch_worker, args)

        self.split_bam_list = []
        for length_dict, gene_indices, read_gene_sets, stats, bam_out_file in results:
            self._merge_length_dict(length_dict)
            self._merge_filter_stats(stats)
            self._merge_saturation_cache(gene_indices, read_gene_sets)
            self.split_bam_list.append(bam_out_file)

        print("ordered fetch scanner finished", flush=True)

    def filter_mrna_reads(self):
        """Preferred public method for automatic read scanning and final BAM output."""
        return self.scan_mrna_reads()

    def fliter_mrna_reads(self):
        """Deprecated spelling-compatible alias for filter_mrna_reads()."""
        return self.filter_mrna_reads()

    def merge_sort_index_bam(self):
        """Concatenate ordered split BAM files and index the final BAM.

        Fetch-mode workers now receive contiguous reference intervals in BAM-header
        order. Therefore, ``SRR*_split_0.bam``, ``SRR*_split_1.bam``, ... are
        already globally ordered. The finalization step can use ``samtools cat``
        instead of ``samtools merge``. This avoids the long post-filtering stall
        caused by coordinate merging of many split BAM files.
        """
        bam_in_list = list(getattr(self, "split_bam_list", []))
        if not bam_in_list:
            raise RuntimeError("No split BAM files were found for final BAM generation.")

        missing_bams = [bam_file for bam_file in bam_in_list if not os.path.exists(bam_file)]
        if missing_bams:
            raise FileNotFoundError("Missing split BAM files:\n" + "\n".join(missing_bams))

        empty_bams = [bam_file for bam_file in bam_in_list if os.path.getsize(bam_file) == 0]
        if empty_bams:
            raise RuntimeError("Empty split BAM files:\n" + "\n".join(empty_bams))

        final_bam = self.output_prefix + ".bam"
        concatenated_bam = self._make_temp_bam(".cat.tmp.bam")

        for stale_file in [
            final_bam,
            final_bam + ".bai",
            final_bam + ".csi",
            os.path.splitext(final_bam)[0] + ".bai",
            os.path.splitext(final_bam)[0] + ".csi",
            concatenated_bam,
        ]:
            if os.path.exists(stale_file):
                os.remove(stale_file)

        try:
            print(f"concatenate {len(bam_in_list)} ordered split BAM file(s) with samtools cat", flush=True)
            self._samtools_cat(concatenated_bam, bam_in_list)
            os.replace(concatenated_bam, final_bam)
            print(f"index final BAM: {final_bam}", flush=True)
            self._samtools_index(final_bam, threads=self.thread)
            self.output_bam = final_bam
        except Exception as cat_error:
            print(
                "fast cat/index failed; fallback to sort -> index. "
                f"Error: {cat_error}",
                flush=True,
            )
            if os.path.exists(final_bam):
                os.remove(final_bam)
            if os.path.exists(concatenated_bam):
                self.sort_index_bam(input_bam=concatenated_bam, output_bam=final_bam, cleanup=True)
            else:
                self._samtools_cat(concatenated_bam, bam_in_list)
                self.sort_index_bam(input_bam=concatenated_bam, output_bam=final_bam, cleanup=True)

        for bam_file in bam_in_list:
            if os.path.exists(bam_file):
                os.remove(bam_file)
            for index_file in [bam_file + ".bai", bam_file + ".csi"]:
                if os.path.exists(index_file):
                    os.remove(index_file)

        print("final BAM concatenation/index done", flush=True)

    def ensure_non_empty_results(self):
        """Fail early with a clear error if no reads survived rpf_Check filtering."""
        total_length_records = sum(sum(v) for v in self.length_dict.values()) if self.length_dict else 0
        passed_filter_reads = int(self.filter_stats.get("passed_filter_reads", 0)) if self.filter_stats else 0

        if total_length_records == 0 or passed_filter_reads == 0:
            try:
                self.write_summary()
            except Exception:
                pass

            raise RuntimeError(
                "No reads remained after rpf_Check filtering. "
                "Please check: "
                "1) whether transcript_id in gene.norm.txt matches BAM reference names; "
                "2) whether -a/--align matches the aligner used to generate the BAM; "
                "3) whether -g 1 removed all reads because unique-mapping tags are missing; "
                "4) whether BAM/SAM references are transcript IDs from rpf_Reference."
            )

    def detect_seq_type(self):
        """Auto-detect the type of sequence profile from plus+minus length counts."""
        if not self.length_dict:
            self.ensure_non_empty_results()

        self.peak_length, strand_counts = max(
            self.length_dict.items(),
            key=lambda item: sum(item[1]),
        )
        self.peak_reads = int(sum(strand_counts))
        display_file = self.output_bam if self.output_bam else self.input_file

        if 23 < self.peak_length < 35:
            print(f"{display_file} is detected to be monosome-seq.\n", flush=True)
            self.profile = "monosome"
            self.mono = [19.5, 40.5]
        elif 53 < self.peak_length < 65:
            print(f"{display_file} is detected to be disome-seq.\n", flush=True)
            self.profile = "disome"
            self.mono = [49.5, 70.5]
        elif 83 < self.peak_length < 95:
            print(f"{display_file} is detected to be trisome-seq.\n", flush=True)
            self.profile = "trisome"
            self.mono = [79.5, 100.5]
        elif 35 < self.peak_length < 53 or 65 < self.peak_length < 83 or 95 < self.peak_length:
            print(
                f"Warning! {display_file} does not fit the empirical length distribution.\n",
                flush=True,
            )
            print(
                "\n".join([
                    "Monosome RPFs peak length is usually ~30 nt.",
                    "Disome RPFs peak length is usually ~60 nt.",
                    "Trisome RPFs peak length is usually ~90 nt.",
                    "Please check the files and run detect_offset.py with specified peak_length.",
                ]),
                flush=True,
            )
        else:
            print(
                f"Warning! Cannot classify profile from peak length: {self.peak_length} nt.\n",
                flush=True,
            )

    def _record_plot_warning(self, message):
        """Record non-fatal plotting warnings without interrupting the QC workflow."""
        message = str(message)
        self.plot_warnings.append(message)
        warning_file = self.output_prefix + "_plot.warning.log"
        with open(warning_file, "a") as warning_out:
            warning_out.write(message + "\n")

    @staticmethod
    def _format_axis_value(value):
        """Format axis values without abbreviated labels.

        The QC report should be directly readable from log files and figures.
        Therefore labels such as 50k are intentionally avoided; integer-like
        values are printed as full numbers.
        """
        value = float(value)
        if not math.isfinite(value):
            return "NA"
        if abs(value) < 1e-12:
            return "0"
        if abs(value - round(value)) < 1e-6:
            return str(int(round(value)))
        if abs(value) >= 100:
            return f"{value:.0f}"
        if abs(value) >= 10:
            return f"{value:.1f}".rstrip("0").rstrip(".")
        return f"{value:.2f}".rstrip("0").rstrip(".")

    @staticmethod
    def _nice_tick_step(value):
        """Return a clean 1/2/5/10-based tick step."""
        value = float(value)
        if value <= 0 or not math.isfinite(value):
            return 1.0
        exponent = math.floor(math.log10(value))
        fraction = value / (10 ** exponent)
        if fraction <= 1:
            nice_fraction = 1
        elif fraction <= 2:
            nice_fraction = 2
        elif fraction <= 5:
            nice_fraction = 5
        else:
            nice_fraction = 10
        return nice_fraction * (10 ** exponent)

    @staticmethod
    def _nice_ticks(max_value, target_ticks=4):
        """Return clean y-axis tick values and a rounded y maximum."""
        max_value = float(max_value)
        if max_value <= 0 or not math.isfinite(max_value):
            return [0.0, 1.0], 1.0
        step = Quality._nice_tick_step(max_value / max(1, int(target_ticks)))
        y_max = math.ceil(max_value / step) * step
        tick_count = int(round(y_max / step))
        ticks = [i * step for i in range(tick_count + 1)]
        return ticks, y_max

    @staticmethod
    def _nice_max(value):
        """Return a rounded upper plotting limit."""
        value = float(value)
        if value <= 0 or not math.isfinite(value):
            return 1.0
        exponent = math.floor(math.log10(value))
        fraction = value / (10 ** exponent)
        if fraction <= 1:
            nice_fraction = 1
        elif fraction <= 2:
            nice_fraction = 2
        elif fraction <= 5:
            nice_fraction = 5
        else:
            nice_fraction = 10
        return nice_fraction * (10 ** exponent)

    @staticmethod
    def _linear_map(value, src_min, src_max, dst_min, dst_max):
        """Map a numeric value from one interval to another interval."""
        if src_max == src_min:
            return (dst_min + dst_max) / 2.0
        return dst_min + (float(value) - src_min) * (dst_max - dst_min) / (src_max - src_min)

    @staticmethod
    def _pdf_escape(text):
        """Escape text for raw PDF content streams."""
        return str(text).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    class _VectorCanvas(object):
        """Small pure-Python vector PDF canvas used as a stable plotting backend.

        The plot style follows the original QC figures, but the PDF is written
        directly to avoid matplotlib/fontTools backend crashes observed on some
        server environments. No auxiliary vector image files are generated.
        """

        _NAMED_COLORS = {
            "black": (0.0, 0.0, 0.0),
            "white": (1.0, 1.0, 1.0),
            "none": None,
            "gray": (0.5, 0.5, 0.5),
            "grey": (0.5, 0.5, 0.5),
        }

        def __init__(self, width, height):
            self.width = float(width)
            self.height = float(height)
            self.pdf_cmds = []
            self.ops = []
            self.rect(0, 0, self.width, self.height, fill="white", stroke="none", line_width=0)

        @classmethod
        def _rgb(cls, color):
            if color is None:
                return None
            color = str(color).strip()
            lower = color.lower()
            if lower in cls._NAMED_COLORS:
                return cls._NAMED_COLORS[lower]
            if color.startswith("#") and len(color) == 7:
                try:
                    return (
                        int(color[1:3], 16) / 255.0,
                        int(color[3:5], 16) / 255.0,
                        int(color[5:7], 16) / 255.0,
                    )
                except ValueError:
                    return (0.0, 0.0, 0.0)
            return (0.0, 0.0, 0.0)

        @classmethod
        def _stroke_cmd(cls, color):
            rgb = cls._rgb(color)
            if rgb is None:
                return ""
            return f"{rgb[0]:.4f} {rgb[1]:.4f} {rgb[2]:.4f} RG"

        @classmethod
        def _fill_cmd(cls, color):
            rgb = cls._rgb(color)
            if rgb is None:
                return ""
            return f"{rgb[0]:.4f} {rgb[1]:.4f} {rgb[2]:.4f} rg"

        def line(self, x1, y1, x2, y2, width=0.8, color="#333333"):
            if color in (None, "none") or width <= 0:
                return
            x1, y1, x2, y2 = map(float, [x1, y1, x2, y2])
            self.ops.append(("line", x1, y1, x2, y2, float(width), color))
            self.pdf_cmds.append(
                f"{width:.3f} w {self._stroke_cmd(color)} {x1:.3f} {y1:.3f} m {x2:.3f} {y2:.3f} l S"
            )

        def polyline(self, points, width=0.9, color="#1F77B4"):
            points = [(float(x), float(y)) for x, y in points]
            if len(points) < 2:
                return
            self.ops.append(("polyline", points, float(width), color))
            pdf_parts = [
                f"{width:.3f} w {self._stroke_cmd(color)} {points[0][0]:.3f} {points[0][1]:.3f} m"
            ]
            pdf_parts.extend(f"{x:.3f} {y:.3f} l" for x, y in points[1:])
            pdf_parts.append("S")
            self.pdf_cmds.append(" ".join(pdf_parts))

        def rect(self, x, y, width, height, fill="none", stroke="#333333", line_width=0.8):
            x, y, width, height = map(float, [x, y, width, height])
            if height < 0:
                y += height
                height = abs(height)
            self.ops.append(("rect", x, y, width, height, fill, stroke, float(line_width)))
            cmds = []
            if fill not in (None, "none"):
                cmds.append(
                    f"q {self._fill_cmd(fill)} {x:.3f} {y:.3f} {width:.3f} {height:.3f} re f Q"
                )
            if stroke not in (None, "none") and line_width > 0:
                cmds.append(
                    f"{line_width:.3f} w {self._stroke_cmd(stroke)} {x:.3f} {y:.3f} {width:.3f} {height:.3f} re S"
                )
            if cmds:
                self.pdf_cmds.append(" ".join(cmds))

        def circle(self, x, y, radius=2.0, fill="#1F77B4", stroke="none", line_width=0.4):
            x, y, radius = float(x), float(y), float(radius)
            self.ops.append(("circle", x, y, float(radius), fill, stroke, float(line_width)))
            c = 0.5522847498 * radius
            path = (
                f"{x + radius:.3f} {y:.3f} m "
                f"{x + radius:.3f} {y + c:.3f} {x + c:.3f} {y + radius:.3f} {x:.3f} {y + radius:.3f} c "
                f"{x - c:.3f} {y + radius:.3f} {x - radius:.3f} {y + c:.3f} {x - radius:.3f} {y:.3f} c "
                f"{x - radius:.3f} {y - c:.3f} {x - c:.3f} {y - radius:.3f} {x:.3f} {y - radius:.3f} c "
                f"{x + c:.3f} {y - radius:.3f} {x + radius:.3f} {y - c:.3f} {x + radius:.3f} {y:.3f} c"
            )
            cmds = []
            if fill not in (None, "none"):
                cmds.append(f"q {self._fill_cmd(fill)} {path} f Q")
            if stroke not in (None, "none") and line_width > 0:
                cmds.append(f"{line_width:.3f} w {self._stroke_cmd(stroke)} {path} S")
            if cmds:
                self.pdf_cmds.append(" ".join(cmds))

        def text(self, x, y, text, size=8, anchor="middle", rotate=0, color="#222222"):
            x, y, size = float(x), float(y), float(size)
            text = str(text)
            self.ops.append(("text", x, y, text, float(size), anchor, float(rotate), color))
            escaped_pdf = Quality._pdf_escape(text)
            if anchor == "middle":
                tx = -0.25 * size * len(text)
            elif anchor == "end":
                tx = -0.5 * size * len(text)
            else:
                tx = 0
            angle = math.radians(float(rotate))
            cos_a = math.cos(angle)
            sin_a = math.sin(angle)
            self.pdf_cmds.append(
                "q "
                f"{self._fill_cmd(color)} "
                f"{cos_a:.6f} {sin_a:.6f} {-sin_a:.6f} {cos_a:.6f} {x:.3f} {y:.3f} cm "
                f"BT /F1 {size:.3f} Tf {tx:.3f} 0 Td ({escaped_pdf}) Tj ET Q"
            )

        def save_pdf(self, pdf_file):
            content = "\n".join(self.pdf_cmds).encode("latin-1", "replace")
            objects = []
            objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
            objects.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
            page = (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.width:.3f} {self.height:.3f}] "
                "/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
            ).encode("latin-1")
            objects.append(page)
            objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
            stream = b"<< /Length " + str(len(content)).encode("ascii") + b" >>\nstream\n" + content + b"\nendstream"
            objects.append(stream)

            with open(pdf_file, "wb") as out:
                out.write(b"%PDF-1.4\n")
                offsets = [0]
                for idx, obj in enumerate(objects, start=1):
                    offsets.append(out.tell())
                    out.write(f"{idx} 0 obj\n".encode("ascii"))
                    out.write(obj)
                    out.write(b"\nendobj\n")
                xref_pos = out.tell()
                out.write(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
                out.write(b"0000000000 65535 f \n")
                for offset in offsets[1:]:
                    out.write(f"{offset:010d} 00000 n \n".encode("ascii"))
                out.write(
                    f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n".encode("ascii")
                )

        @classmethod
        def _png_rgb(cls, color):
            """Convert a PDF-style color value to an RGB tuple for PNG output."""
            rgb = cls._rgb(color)
            if rgb is None:
                return None
            return tuple(int(round(max(0.0, min(1.0, channel)) * 255)) for channel in rgb)

        @staticmethod
        def _load_font(size_px):
            """Load a readable TrueType font for PNG rendering."""
            size_px = max(8, int(round(size_px)))
            candidates = [
                "DejaVuSans.ttf",
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                "/usr/share/fonts/dejavu/DejaVuSans.ttf",
                "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
                "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
            ]
            for candidate in candidates:
                try:
                    return ImageFont.truetype(candidate, size_px)
                except Exception:
                    continue
            return ImageFont.load_default()

        def save_png(self, png_file, scale=2):
            """Rasterize the recorded vector drawing commands to PNG."""
            scale = max(1, int(scale))
            width_px = int(round(self.width * scale))
            height_px = int(round(self.height * scale))
            image = Image.new("RGB", (width_px, height_px), (255, 255, 255))
            draw = ImageDraw.Draw(image)

            def px(value):
                return float(value) * scale

            def py(value):
                return float(self.height - float(value)) * scale

            for op in self.ops:
                kind = op[0]

                if kind == "line":
                    _, x1, y1, x2, y2, line_width, color = op
                    color_rgb = self._png_rgb(color)
                    if color_rgb is not None:
                        draw.line(
                            [(px(x1), py(y1)), (px(x2), py(y2))],
                            fill=color_rgb,
                            width=max(1, int(round(line_width * scale))),
                        )

                elif kind == "polyline":
                    _, points, line_width, color = op
                    color_rgb = self._png_rgb(color)
                    if color_rgb is not None and len(points) >= 2:
                        draw.line(
                            [(px(x), py(y)) for x, y in points],
                            fill=color_rgb,
                            width=max(1, int(round(line_width * scale))),
                            joint="curve",
                        )

                elif kind == "rect":
                    _, x, y, width, height, fill, stroke, line_width = op
                    fill_rgb = self._png_rgb(fill) if fill not in (None, "none") else None
                    stroke_rgb = self._png_rgb(stroke) if stroke not in (None, "none") else None
                    xy = [px(x), py(y + height), px(x + width), py(y)]
                    draw.rectangle(
                        xy,
                        fill=fill_rgb,
                        outline=stroke_rgb,
                        width=max(1, int(round(line_width * scale))) if stroke_rgb else 1,
                    )

                elif kind == "circle":
                    _, x, y, radius, fill, stroke, line_width = op
                    fill_rgb = self._png_rgb(fill) if fill not in (None, "none") else None
                    stroke_rgb = self._png_rgb(stroke) if stroke not in (None, "none") else None
                    xy = [px(x - radius), py(y + radius), px(x + radius), py(y - radius)]
                    draw.ellipse(
                        xy,
                        fill=fill_rgb,
                        outline=stroke_rgb,
                        width=max(1, int(round(line_width * scale))) if stroke_rgb else 1,
                    )

                elif kind == "text":
                    _, x, y, text, size, anchor, rotate, color = op
                    font = self._load_font(size * scale)
                    color_rgb = self._png_rgb(color) or (0, 0, 0)
                    text = str(text)
                    bbox = font.getbbox(text)
                    text_width = bbox[2] - bbox[0]
                    text_height = bbox[3] - bbox[1]
                    padding = max(4, int(round(3 * scale)))
                    text_image = Image.new(
                        "RGBA",
                        (text_width + padding * 2, text_height + padding * 2),
                        (255, 255, 255, 0),
                    )
                    text_draw = ImageDraw.Draw(text_image)
                    text_draw.text(
                        (padding - bbox[0], padding - bbox[1]),
                        text,
                        fill=color_rgb + (255,),
                        font=font,
                    )
                    if abs(float(rotate)) > 1e-6:
                        text_image = text_image.rotate(
                            -float(rotate),
                            expand=True,
                            resample=Image.Resampling.BICUBIC,
                        )

                    anchor_x = px(x)
                    anchor_y = py(y)
                    if anchor == "middle":
                        paste_x = int(round(anchor_x - text_image.width / 2))
                    elif anchor == "end":
                        paste_x = int(round(anchor_x - text_image.width))
                    else:
                        paste_x = int(round(anchor_x))
                    paste_y = int(round(anchor_y - text_image.height / 2))
                    image.paste(text_image, (paste_x, paste_y), text_image)

            image.save(png_file, format="PNG")

    def _draw_panel_axes(self, canvas, left, bottom, width, height, title, xlabel, ylabel, x_ticks, x_labels, y_max):
        """Draw a clean XY axis panel with full numeric y-axis labels."""
        y_ticks, y_max = self._nice_ticks(y_max, target_ticks=4)
        for value in y_ticks:
            y = self._linear_map(value, 0, y_max, bottom, bottom + height)
            if value > 0:
                canvas.line(left, y, left + width, y, width=0.25, color="#D9D9D9")
            canvas.line(left - 3, y, left, y, width=0.6, color="#333333")
            canvas.text(left - 9, y - 2, self._format_axis_value(value), size=13, anchor="end", color="#222222")

        canvas.line(left, bottom, left + width, bottom, width=0.8, color="#333333")
        canvas.line(left, bottom, left, bottom + height, width=0.8, color="#333333")
        canvas.text(left + width / 2, bottom + height + 22, title, size=14, anchor="middle")
        canvas.text(left + width / 2, bottom - 46, xlabel, size=14, anchor="middle")
        canvas.text(left - 70, bottom + height / 2, ylabel, size=14, anchor="middle", rotate=90)

        for x, label in zip(x_ticks, x_labels):
            canvas.line(x, bottom, x, bottom - 3, width=0.6, color="#333333")
            canvas.text(x, bottom - 20, label, size=13, anchor="middle", rotate=90 if len(str(label)) > 2 else 0)

        return y_max

    def _save_vector_canvas(self, canvas, pdf_file, png_file=None):
        """Save one vector canvas as PDF and optional PNG."""
        canvas.save_pdf(pdf_file)
        if png_file:
            canvas.save_png(png_file, scale=2)

    def _draw_horizontal_dashed_reference(self, canvas, left, right, y, color="#666666", line_width=0.8, dash=6.0, gap=4.0):
        """Draw a horizontal dashed reference line on the vector canvas."""
        left = float(left)
        right = float(right)
        y = float(y)
        dash = max(1.0, float(dash))
        gap = max(1.0, float(gap))
        x = left
        while x < right:
            x2 = min(x + dash, right)
            canvas.line(x, y, x2, y, width=line_width, color=color)
            x += dash + gap

    def draw_the_length_distr(self, sorted_length):
        """Draw plus/minus RPF length distribution as vector PDF and PNG."""
        out_pdf = self.output_prefix + "_length_distribution.pdf"
        out_png = self.output_prefix + "_length_distribution.png"

        if not sorted_length:
            raise RuntimeError("No length distribution data are available for plotting.")

        length_df = pd.DataFrame.from_dict(
            sorted_length,
            orient="index",
            columns=["Plus", "Minus"],
        ).sort_index()
        x_values = length_df.index.to_numpy(dtype=float)
        plus_values = length_df["Plus"].to_numpy(dtype=float)
        minus_values = length_df["Minus"].to_numpy(dtype=float)

        canvas = self._VectorCanvas(660, 560)
        canvas.text(330, 535, "RPFs length distribution", size=17, anchor="middle")

        panels = [
            ("Plus strand", plus_values, 315, "#1F77B4"),
            ("Minus strand", minus_values, 95, "#D62728"),
        ]
        x_min, x_max = float(self.mono[0]), float(self.mono[1])
        visible_mask = (x_values >= x_min) & (x_values <= x_max)
        visible_x = x_values[visible_mask] if np.any(visible_mask) else x_values
        tick_step = max(1, int(math.ceil(len(visible_x) / 18)))
        tick_values = visible_x[::tick_step]

        for title, y_values, bottom, line_color in panels:
            left, width, height = 95, 470, 160
            y_max = self._nice_max(np.nanmax(y_values) if len(y_values) else 1)
            x_ticks = [self._linear_map(x, x_min, x_max, left, left + width) for x in tick_values]
            x_labels = [str(int(x)) for x in tick_values]
            y_max = self._draw_panel_axes(
                canvas, left, bottom, width, height, title,
                "RPFs length (nt)", "RPFs number", x_ticks, x_labels, y_max,
            )
            points = []
            for x, y in zip(x_values, y_values):
                if x < x_min or x > x_max:
                    continue
                px = self._linear_map(x, x_min, x_max, left, left + width)
                py = self._linear_map(y, 0, y_max, bottom, bottom + height)
                points.append((px, py))
            canvas.polyline(points, width=1.4, color=line_color)
            for px, py in points:
                canvas.circle(px, py, radius=1.8, fill=line_color)

        self._save_vector_canvas(canvas, out_pdf, out_png)

    def write_length_distr(self):
        """Write and plot plus/minus read-length distribution."""
        sorted_length = dict(sorted(self.length_dict.items(), key=lambda length: length[0]))
        with open(self.output_prefix + "_length_distribution.txt", "w") as length_out:
            length_out.write("\t".join(["Length", "Plus", "Minus", "Total"]) + "\n")
            for reads_length, reads_num in sorted_length.items():
                total_num = int(reads_num[0]) + int(reads_num[1])
                length_out.write(
                    "\t".join([
                        str(reads_length),
                        str(reads_num[0]),
                        str(reads_num[1]),
                        str(total_num),
                    ]) + "\n"
                )

        self.draw_the_length_distr(sorted_length)


    def _count_saturation_reads(self):
        """Return the number of read units cached for saturation analysis.

        Unique-read mode stores one integer transcript index per alignment.
        Multi-mapping mode stores one transcript-index set per query name.
        This method is intentionally lightweight because it may be called before
        the full saturation matrix is calculated.
        """
        if self.saturation_gene_indices:
            return int(len(self.saturation_gene_indices))
        if self.saturation_read_gene_sets:
            return int(len(self.saturation_read_gene_sets))
        return 0

    def _summary_dict(self):
        """Return a flat, non-redundant rpf_Check summary."""
        plus_reads = int(sum(v[0] for v in self.length_dict.values())) if self.length_dict else 0
        minus_reads = int(sum(v[1] for v in self.length_dict.values())) if self.length_dict else 0
        total_length_reads = plus_reads + minus_reads

        return OrderedDict([
            ("sample", os.path.basename(self.input_file)),
            ("input_alignment", self.input_file),
            ("output_bam", self.output_bam if self.output_bam else "NA"),
            ("parallel_workers", int(self.thread)),
            ("aligner", self.align),
            ("unique_mapped_only", bool(self.tag == 1)),
            ("longest_transcript_only", bool(self.longest)),
            ("saturation", bool(self.saturation_flag)),
            ("transcripts_input", int(self.transcript_total_input)),
            ("transcripts_kept", int(self.transcript_total_kept)),
            ("transcripts_removed", int(self.transcript_removed_by_longest)),
            ("bam_reference_matched", int(self.bam_reference_match_count)),
            ("bam_reference_missing", int(self.filter_stats.get("fetch_missing_reference", 0))),
            ("total_records_scanned", int(self.filter_stats.get("total_records_scanned", 0))),
            ("passed_filter_reads", int(self.filter_stats.get("passed_filter_reads", 0))),
            ("written_reads", int(self.filter_stats.get("written_reads", 0))),
            ("plus_reads", plus_reads),
            ("minus_reads", minus_reads),
            ("total_length_counted_reads", total_length_reads),
            ("unmapped_reads", int(self.filter_stats.get("unmapped_reads", 0))),
            ("qc_fail_reads", int(self.filter_stats.get("qc_fail_reads", 0))),
            ("secondary_reads", int(self.filter_stats.get("secondary_reads", 0))),
            ("supplementary_reads", int(self.filter_stats.get("supplementary_reads", 0))),
            ("duplicate_reads", int(self.filter_stats.get("duplicate_reads", 0))),
            ("removed_secondary_reads", int(self.filter_stats.get("removed_secondary_reads", 0))),
            ("removed_supplementary_reads", int(self.filter_stats.get("removed_supplementary_reads", 0))),
            ("removed_duplicate_reads", int(self.filter_stats.get("removed_duplicate_reads", 0))),
            ("removed_multi_mapping_reads", int(self.filter_stats.get("removed_multi_mapping_reads", 0))),
            ("missing_unique_tag_reads", int(self.filter_stats.get("missing_unique_tag_reads", 0))),
            ("missing_length_reads", int(self.filter_stats.get("missing_length_reads", 0))),
            ("peak_length", int(self.peak_length) if self.peak_length else 0),
            ("peak_reads", int(self.peak_reads) if self.peak_reads else 0),
            ("profile", self.profile if self.profile else "undetermined"),
            ("saturation_reads", self._count_saturation_reads() if self.saturation_flag else "NA"),
        ])

    def write_summary(self):
        """Write rpf_Check summary as JSON."""
        summary = self._summary_dict()

        out_json = self.output_prefix + "_rpf_check.summary.json"
        with open(out_json, "w") as out:
            json.dump(summary, out, indent=2)
            out.write("\n")

    def _get_saturation_units(self):
        """Return read-level saturation units represented by integer gene indices."""
        if self.saturation_gene_indices:
            return list(self.saturation_gene_indices)
        if self.saturation_read_gene_sets:
            return [tuple(sorted(gene_set)) for gene_set in self.saturation_read_gene_sets.values()]
        return []

    def rpf_saturation(self):
        """Calculate RPF saturation using integer gene indices.

        Workflow
        --------
        1. Use fixed seed 5201314.
        2. Shuffle read units once.
        3. Use nested prefixes at 10%, 20%, ..., 90%.
        4. Incrementally update per-transcript read counts.
        """
        units = self._get_saturation_units()
        mapped_reads = len(units)
        if mapped_reads == 0:
            raise RuntimeError(
                "Saturation was requested, but no reads are available for saturation. "
                "Please check whether reads survived filtering."
            )

        rng = random.Random(self.saturation_seed)
        rng.shuffle(units)

        steps = [int(i * 0.10 * mapped_reads) for i in range(1, 10)]
        steps = [max(1, min(step, mapped_reads)) for step in steps]

        running_counts = np.zeros(len(self.index_to_gene), dtype=np.uint32)
        covered_genes = set()
        self.saturation_matrix = np.zeros((len(self.index_to_gene), 9), dtype=np.uint32)
        self.saturation = []
        cursor = 0

        for site, step in enumerate(steps):
            for unit in units[cursor:step]:
                if isinstance(unit, int):
                    running_counts[unit] += 1
                    covered_genes.add(unit)
                else:
                    for gene_index in unit:
                        running_counts[gene_index] += 1
                        covered_genes.add(gene_index)
            cursor = step

            self.saturation.append(len(covered_genes))
            self.saturation_matrix[:, site] = running_counts

        self.mrna_dict = OrderedDict(
            (
                transcript_id,
                self.saturation_matrix[idx, :].astype(int).tolist(),
            )
            for idx, transcript_id in enumerate(self.index_to_gene)
        )

    def rpf_saturation_thread(self):
        """Compatibility alias for the optimized deterministic saturation method."""
        return self.rpf_saturation()

    def draw_gene_saturation(self):
        """Draw gene-level saturation curves as vector PDF and PNG."""
        out_pdf = self.output_prefix + "_gene_saturation.pdf"
        out_png = self.output_prefix + "_gene_saturation.png"
        out_gene = self.output_prefix + "_gene_saturation.txt"

        total_gene_num = len(self.index_to_gene)
        gene_df = pd.DataFrame(self.saturation + [total_gene_num], columns=["Count"])
        gene_df["Part"] = self.x_ticks.tolist() + [0]
        gene_df = gene_df[["Part", "Count"]]
        gene_df.to_csv(out_gene, sep="\t", index=False)

        canvas = self._VectorCanvas(780, 400)
        canvas.text(390, 372, "Gene saturation", size=17, anchor="middle")
        panels = [
            ("covered gene saturation", self.saturation, 85, "gene number", "#4C78A8"),
            ("uncovered gene saturation", [total_gene_num - i for i in self.saturation], 410, "gene number", "#F58518"),
        ]

        for title, values, left, ylabel, plot_color in panels:
            bottom, width, height = 82, 245, 205
            y_max = self._nice_max(max(max(values) if values else 1, total_gene_num))
            x_min, x_max = 5.0, 95.0
            x_positions = [self._linear_map(x, x_min, x_max, left, left + width) for x in self.x_ticks]
            y_max = self._draw_panel_axes(
                canvas, left, bottom, width, height, title,
                "reads proportion (%)", ylabel,
                x_positions, [str(i) for i in self.x_ticks], y_max,
            )

            ref_y = self._linear_map(total_gene_num, 0, y_max, bottom, bottom + height)
            self._draw_horizontal_dashed_reference(
                canvas,
                left,
                left + width,
                ref_y,
                color="#7F7F7F",
                line_width=0.8,
                dash=6.0,
                gap=4.0,
            )
            canvas.text(
                left + width - 2,
                min(bottom + height + 10, ref_y + 8),
                f"Total genes: {int(total_gene_num)}",
                size=11,
                anchor="end",
                color="#555555",
            )

            bar_width = max(4.0, width / 42.0)
            points = []
            for x, value in zip(x_positions, values):
                y = self._linear_map(value, 0, y_max, bottom, bottom + height)
                canvas.rect(x - bar_width / 2, bottom, bar_width, y - bottom, fill=plot_color, stroke=plot_color, line_width=0.4)
                points.append((x, y))
            canvas.polyline(points, width=1.2, color="#333333")
            for x, y in points:
                canvas.circle(x, y, radius=1.6, fill="#333333")

        self._save_vector_canvas(canvas, out_pdf, out_png)

    @staticmethod
    def _boxplot_stats(values):
        """Return boxplot statistics for a numeric vector."""
        arr = np.asarray(values, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return None
        q1, median, q3 = np.percentile(arr, [25, 50, 75])
        iqr = q3 - q1
        low_bound = q1 - 1.5 * iqr
        high_bound = q3 + 1.5 * iqr
        lower_values = arr[arr >= low_bound]
        upper_values = arr[arr <= high_bound]
        whisker_low = float(np.min(lower_values)) if lower_values.size else float(np.min(arr))
        whisker_high = float(np.max(upper_values)) if upper_values.size else float(np.max(arr))
        outliers = arr[(arr < whisker_low) | (arr > whisker_high)]
        if outliers.size > 120:
            outliers = np.sort(outliers)[::max(1, int(outliers.size / 120))]
        return {
            "q1": float(q1),
            "median": float(median),
            "q3": float(q3),
            "whisker_low": whisker_low,
            "whisker_high": whisker_high,
            "outliers": outliers.astype(float),
        }

    def _draw_boxplot_panel(self, canvas, data_frame, left, bottom, width, height, title):
        """Draw one reads-saturation boxplot panel with readable count-scale labels."""
        canvas.text(left + width / 2, bottom + height + 24, title, size=14, anchor="middle")
        canvas.text(left + width / 2, bottom - 48, "reads proportion (%)", size=13, anchor="middle")
        canvas.text(left - 74, bottom + height / 2, "Reads count", size=13, anchor="middle", rotate=90)

        if data_frame.empty:
            canvas.line(left, bottom, left + width, bottom, width=0.8, color="#333333")
            canvas.line(left, bottom, left, bottom + height, width=0.8, color="#333333")
            canvas.text(left + width / 2, bottom + height / 2, "No data", size=13, anchor="middle")
            return

        stats_by_col = []
        all_values = []
        for col in data_frame.columns[:9]:
            values = pd.to_numeric(data_frame[col], errors="coerce").dropna().to_numpy(dtype=float)
            values = np.maximum(values, 0.0)
            logged = np.log10(values + 1.0)
            stats = self._boxplot_stats(logged)
            stats_by_col.append(stats)
            if values.size:
                all_values.extend(values[np.isfinite(values)].tolist())

        max_count = max(all_values) if all_values else 1.0
        max_log = max(1.0, math.ceil(math.log10(max_count + 1.0)))

        y_tick_counts = [0]
        power = 0
        while 10 ** power <= max_count:
            y_tick_counts.append(10 ** power)
            power += 1
        y_tick_counts = sorted(set(y_tick_counts))

        for count_value in y_tick_counts:
            y = self._linear_map(math.log10(count_value + 1.0), 0, max_log, bottom, bottom + height)
            if count_value > 0:
                canvas.line(left, y, left + width, y, width=0.25, color="#D9D9D9")
            canvas.line(left - 3, y, left, y, width=0.6, color="#333333")
            canvas.text(left - 9, y - 2, self._format_axis_value(count_value), size=10, anchor="end")

        canvas.line(left, bottom, left + width, bottom, width=0.8, color="#333333")
        canvas.line(left, bottom, left, bottom + height, width=0.8, color="#333333")

        n = len(stats_by_col)
        box_width = max(5.0, width / (n * 2.8))
        for idx, stats in enumerate(stats_by_col):
            x = left + (idx + 0.5) * width / n
            tick_label = str(self.x_ticks[idx]) if idx < len(self.x_ticks) else str(idx + 1)
            canvas.line(x, bottom, x, bottom - 3, width=0.5, color="#333333")
            canvas.text(x, bottom - 20, tick_label, size=10, anchor="middle", rotate=90)
            if stats is None:
                continue

            def map_y(value):
                return self._linear_map(value, 0, max_log, bottom, bottom + height)

            y_q1 = map_y(stats["q1"])
            y_med = map_y(stats["median"])
            y_q3 = map_y(stats["q3"])
            y_low = map_y(stats["whisker_low"])
            y_high = map_y(stats["whisker_high"])
            canvas.line(x, y_low, x, y_q1, width=0.6, color="#333333")
            canvas.line(x, y_q3, x, y_high, width=0.6, color="#333333")
            canvas.line(x - box_width / 3, y_low, x + box_width / 3, y_low, width=0.6, color="#333333")
            canvas.line(x - box_width / 3, y_high, x + box_width / 3, y_high, width=0.6, color="#333333")
            canvas.rect(x - box_width / 2, y_q1, box_width, max(0.5, y_q3 - y_q1), fill="#DDEAF7", stroke="#1F77B4", line_width=0.6)
            canvas.line(x - box_width / 2, y_med, x + box_width / 2, y_med, width=0.8, color="#D62728")
            for outlier in stats["outliers"]:
                canvas.circle(x, map_y(outlier), radius=0.75, fill="#555555")

    def draw_rpf_saturation(self):
        """Draw read-count saturation boxplots by expression quartile as vector PDF and PNG."""
        out_pdf = self.output_prefix + "_reads_saturation.pdf"
        out_png = self.output_prefix + "_reads_saturation.png"
        out_rpf = self.output_prefix + "_reads_saturation.txt"

        if self.saturation_matrix is None:
            raise RuntimeError("Saturation matrix is not available. Run rpf_saturation() first.")

        mrna_df = pd.DataFrame(
            self.saturation_matrix,
            index=self.index_to_gene,
            columns=[str(i) for i in self.x_ticks],
        )
        mrna_df.loc[:, "mean"] = mrna_df.mean(axis=1)
        mrna_df = mrna_df.sort_values(["mean"], ascending=True)
        mrna_df.to_csv(out_rpf, sep="\t")

        quantile_list = mrna_df["mean"].quantile([0.25, 0.5, 0.75])
        mrna_0_25 = mrna_df[mrna_df["mean"] < quantile_list.iloc[0]]
        mrna_25_50 = mrna_df[
            (quantile_list.iloc[0] <= mrna_df["mean"])
            & (mrna_df["mean"] < quantile_list.iloc[1])
        ]
        mrna_50_75 = mrna_df[
            (quantile_list.iloc[1] <= mrna_df["mean"])
            & (mrna_df["mean"] < quantile_list.iloc[2])
        ]
        mrna_75_100 = mrna_df[quantile_list.iloc[2] <= mrna_df["mean"]]

        canvas = self._VectorCanvas(1320, 450)
        canvas.text(660, 418, "Reads saturation", size=17, anchor="middle")
        groups = [
            (mrna_0_25, "0-25%"),
            (mrna_25_50, "25-50%"),
            (mrna_50_75, "50-75%"),
            (mrna_75_100, "75-100%"),
        ]
        for idx, (group_df, title) in enumerate(groups):
            self._draw_boxplot_panel(
                canvas,
                group_df.iloc[:, 0:9],
                left=90 + idx * 305,
                bottom=92,
                width=220,
                height=250,
                title=title,
            )

        self._save_vector_canvas(canvas, out_pdf, out_png)
