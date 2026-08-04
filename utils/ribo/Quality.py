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

# matplotlib is imported lazily (see ``_ensure_plotting_backend``) so that the
# forked scan workers used by ``scan_mrna_reads_fetch`` never inherit
# matplotlib's C-level state (ft2font/freetype/fontTools) through fork. Loading
# it too early makes the process pool non-deterministically crash with SIGSEGV
# on some servers.
matplotlib = None  # type: ignore[assignment]
plt = None  # type: ignore[assignment]


def _ensure_plotting_backend():
    """Import matplotlib (AGG backend) on first figure drawing, in place."""
    global matplotlib, plt
    if plt is not None:
        return
    import matplotlib as _matplotlib

    _matplotlib.use("AGG")
    import matplotlib.pyplot as _plt

    matplotlib = _matplotlib
    plt = _plt


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
            # Open single-threaded here: this handle lives in the parent before
            # the process pool is forked, and an active htslib thread pool at
            # fork time is a known cause of non-deterministic SIGSEGV in the
            # worker processes.
            with pysam.AlignmentFile(self.sample_file, self.sample_format) as bam_in:
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


    # Shared academic style for the QC figures, kept next to the plotting
    # code (rather than at the top of the module) so the style parameters
    # stay with the figures they affect.
    _QC_PLOT_STYLE = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "text.color": "black",
        "axes.edgecolor": "black",
        "axes.linewidth": 0.9,
        "axes.facecolor": "white",
        "axes.labelcolor": "black",
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "legend.frameon": True,
        "legend.edgecolor": "#BBBBBB",
        "legend.fancybox": False,
        "legend.fontsize": 9,
        "figure.facecolor": "white",
        "figure.dpi": 120,
        "savefig.dpi": 300,
    }

    # Blue-orange contrast palette shared by the QC figures.
    _QC_BLUE = "#2166AC"
    _QC_BLUE_LIGHT = "#4393C3"
    _QC_ORANGE = "#E08214"
    _QC_ORANGE_DARK = "#B84507"

    def draw_the_length_distr(self, sorted_length):
        """Draw plus/minus RPF length distribution as PDF and PNG."""
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

        x_min, x_max = float(self.mono[0]), float(self.mono[1])
        mask = (x_values >= x_min) & (x_values <= x_max)
        if not np.any(mask):
            mask = np.ones_like(x_values, dtype=bool)
        x = x_values[mask]
        tick_step = max(1, int(math.ceil(len(x) / 10)))
        tick_values = x[::tick_step].astype(int)

        _ensure_plotting_backend()
        matplotlib.rcParams.update(self._QC_PLOT_STYLE)
        fig, axes = plt.subplots(
            nrows=1, ncols=2, figsize=(7.4, 3.5), sharey=True, constrained_layout=True
        )
        fig.suptitle("RPFs length distribution", fontsize=13)
        series = [
            (axes[0], "Plus strand", plus_values[mask], self._QC_BLUE),
            (axes[1], "Minus strand", minus_values[mask], self._QC_ORANGE),
        ]
        # The two panels share one y-axis, so the limit must be computed from
        # both strands together; otherwise the last panel's limit overwrites
        # the other and clips its peak outside the axes.
        y_axis_max = max(
            1.0,
            float(np.max(plus_values[mask])),
            float(np.max(minus_values[mask])),
        )
        for ax, title, y, color in series:
            ax.plot(x, y, color=color, linewidth=1.4, marker="o", markersize=3.0)
            ax.fill_between(x, 0, y, color=color, alpha=0.06, linewidth=0)
            ax.set_title(title, fontsize=11, pad=6)
            ax.set_xlabel("Read length (nt)", fontsize=10)
            ax.set_ylabel("Number of RPFs", fontsize=10)
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(0, y_axis_max * 1.18)
            ax.set_xticks(tick_values)

            # mark the most abundant read length
            peak_idx = int(np.argmax(y))
            peak_x = float(x[peak_idx])
            ax.axvline(peak_x, color="#888888", linestyle="--", linewidth=0.9)
            ax.annotate(
                f"{peak_x:.0f} nt",
                xy=(peak_x, float(y[peak_idx])),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                fontsize=9,
            )

        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.06)
        fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.06)
        plt.close(fig)

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
        """Draw gene-level saturation curves as PDF and PNG."""
        out_pdf = self.output_prefix + "_gene_saturation.pdf"
        out_png = self.output_prefix + "_gene_saturation.png"
        out_gene = self.output_prefix + "_gene_saturation.txt"

        total_gene_num = len(self.index_to_gene)
        gene_df = pd.DataFrame(self.saturation + [total_gene_num], columns=["Count"])
        gene_df["Part"] = self.x_ticks.tolist() + [0]
        gene_df = gene_df[["Part", "Count"]]
        gene_df.to_csv(out_gene, sep="\t", index=False)

        _ensure_plotting_backend()
        matplotlib.rcParams.update(self._QC_PLOT_STYLE)
        fig, axes = plt.subplots(
            nrows=1, ncols=2, figsize=(7.6, 3.5), constrained_layout=True
        )
        fig.suptitle("Gene saturation", fontsize=13)
        x_positions = np.arange(len(self.saturation))
        covered = np.asarray(self.saturation, dtype=float)
        uncovered = total_gene_num - covered
        series = [
            (axes[0], "Covered genes", covered, self._QC_BLUE),
            (axes[1], "Uncovered genes", uncovered, self._QC_ORANGE),
        ]
        for ax, title, y, color in series:
            ax.axhline(total_gene_num, color="#888888", linestyle="--", linewidth=0.9)
            ax.plot(x_positions, y, color=color, linewidth=1.5, marker="o", markersize=3.6)
            ax.fill_between(x_positions, 0, y, color=color, alpha=0.05, linewidth=0)
            ax.text(
                0.98,
                0.96,
                f"Total genes: {int(total_gene_num):,}",
                transform=ax.transAxes,
                fontsize=8.5,
                color="black",
                ha="right",
                va="top",
            )
            ax.set_title(title, fontsize=11, pad=6)
            ax.set_xlabel("Reads proportion (%)", fontsize=10)
            ax.set_ylabel("Gene number", fontsize=10)
            ax.set_xticks(x_positions, [str(i) for i in self.x_ticks])
            ax.set_ylim(bottom=0)

        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.06)
        fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.06)
        plt.close(fig)

    def draw_rpf_saturation(self):
        """Draw read-count saturation boxplots by expression quartile as PDF and PNG."""
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

        _ensure_plotting_backend()
        matplotlib.rcParams.update(self._QC_PLOT_STYLE)
        fig, axes = plt.subplots(
            nrows=1, ncols=4, figsize=(11.8, 3.3), sharey=True, constrained_layout=True
        )
        fig.suptitle("Reads saturation", fontsize=13)
        x_positions = np.arange(len(self.x_ticks))
        groups = [
            (mrna_0_25, "0-25%"),
            (mrna_25_50, "25-50%"),
            (mrna_50_75, "50-75%"),
            (mrna_75_100, "75-100%"),
        ]

        # One box per reads proportion, split by expression-quartile group.
        logged_groups = []
        max_count = 1.0
        max_log = 1.0
        for group_df, _ in groups:
            matrix = group_df.iloc[:, 0:9].to_numpy(dtype=float)
            cols = []
            for col_idx in range(matrix.shape[1]):
                col_values = matrix[:, col_idx]
                col_values = np.maximum(col_values[np.isfinite(col_values)], 0.0)
                if col_values.size == 0:
                    continue
                cols.append(np.log10(col_values + 1.0))
                max_count = max(max_count, float(col_values.max()))
                max_log = max(max_log, float(np.log10(col_values.max() + 1.0)))
            logged_groups.append(cols)

        # Shared log10(count + 1) y-axis so all four panels are comparable.
        tick_counts = [0] + [10 ** p for p in range(0, int(math.ceil(math.log10(max_count))) + 1)]
        y_ticks_log = [math.log10(c + 1.0) for c in tick_counts]
        y_tick_labels = [self._format_axis_value(c) for c in tick_counts]

        box_colors = [self._QC_BLUE, self._QC_BLUE_LIGHT, self._QC_ORANGE, self._QC_ORANGE_DARK]
        for idx, (logged, title) in enumerate(zip(logged_groups, [title for _, title in groups])):
            ax = axes[idx]
            ax.set_title(title, fontsize=11, pad=6)
            ax.set_xlabel("Reads proportion (%)", fontsize=10)
            if idx == 0:
                ax.set_ylabel("Reads count", fontsize=10)
            ax.set_xticks(x_positions)
            ax.set_xticklabels([str(i) for i in self.x_ticks])
            ax.set_yticks(y_ticks_log)
            ax.set_yticklabels(y_tick_labels)
            ax.set_ylim(0, max_log * 1.05)
            if not logged:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", fontsize=10)
                continue
            color = box_colors[idx]
            ax.boxplot(
                logged,
                positions=x_positions[: len(logged)],
                widths=0.55,
                patch_artist=True,
                medianprops=dict(color="black", linewidth=1.1),
                whiskerprops=dict(color=color, linewidth=0.9),
                capprops=dict(color=color, linewidth=0.9),
                boxprops=dict(
                    facecolor=matplotlib.colors.to_rgba(color, alpha=0.10),
                    edgecolor=color,
                    linewidth=0.9,
                ),
                flierprops=dict(
                    marker="o",
                    markersize=1.6,
                    markerfacecolor="#666666",
                    markeredgecolor="none",
                    alpha=0.5,
                ),
            )

        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.06)
        fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.06)
        plt.close(fig)
