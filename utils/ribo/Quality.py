#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Rensc, modified by local fork
# Date: 2026-06-22
# Version: 0.2.7
# Function: Provide RiboParser core functions for Quality analysis.

import json
import os
import os.path
import random
import sys
from collections import Counter
from collections import OrderedDict
from itertools import islice
from multiprocessing import Pool

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


class Quality(object):
    def __init__(self, args):
        # opts for mrna file import
        self.longest = args.longest
        self.mrna_file = args.transcript
        self.mrna_dict = OrderedDict()
        self.length_dict = {}

        # opts for bam file import
        self.thread = max(1, int(args.thread))
        self.pysam_input = None
        self.pysam_output = None
        self.sample_file = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]
        self.bam_seq_dict = {}
        self.sample_dict = None
        self.tag = args.tag
        self.reverse = args.reverse
        self.align = args.align

        # opts for bam length summary
        self.peak_length = 0
        self.peak_reads = 0
        self.mono = 1
        self.profile = None

        # opts for rpf saturation
        self.saturation_flag = args.saturation
        self.saturation = []
        self.saturation_seed = 5201314
        self.x_ticks = np.array([i * 10 for i in range(1, 10)])

        # opts for file output
        self.output_prefix = args.output

        # opts for QC reports
        self.filter_stats = Counter()
        self.transcript_total_input = 0
        self.transcript_total_kept = 0
        self.transcript_removed_by_longest = 0

    def read_transcript(self):
        """Read transcript annotation and optionally keep only longest transcripts.

        Expected RiboParser gene.norm.txt columns include:
        chromosome, gene_id, transcript_id, start, end, utr5_length, cds_length,
        utr3_length, strand, rep_transcript, modified.

        If -l/--longest is enabled:
        1. Prefer rows marked as representative transcript if available and unique per gene.
        2. Otherwise choose the largest CDS length per gene.
        3. Tie-break by total transcript length and transcript_id.
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
            # Backward-compatible fallback for files without a header.
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
                trans_df[col] = pd.to_numeric(trans_df[col], errors="coerce").fillna(0).astype(int)
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

        self.mrna_dict = OrderedDict()
        for transcript_id in trans_df["transcript_id"].astype(str):
            if transcript_id and transcript_id.lower() != "nan":
                self.mrna_dict[transcript_id] = [0] * 9

        if not self.mrna_dict:
            raise RuntimeError(
                "No valid transcript_id was loaded from the annotation file. "
                "Please check whether the transcript file is the RiboParser gene.norm.txt format."
            )

    @staticmethod
    def _truthy_series(series):
        return series.astype(str).str.lower().isin(["true", "1", "yes", "y", "t"])

    def _select_longest_transcripts(self, trans_df):
        df = trans_df.copy()

        # Prefer representative transcripts generated by rpf_Reference, only when unique per gene.
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

    def sort_index_bam(self):
        # check the data type, to support the SAM/BAM file format both.
        if self.sample_format.lower() == ".bam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            self.sample_format = "rb"
        elif self.sample_format.lower() == ".sam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            self.sample_format = "r"
        else:
            print("Unknown file format, please input the correct bam or sam file.", flush=True)
            sys.exit()

        # check the bam index file
        if self.sample_format == "rb" and (
            os.path.exists(self.sample_file + ".bai")
            or os.path.exists(os.path.splitext(self.sample_file)[0] + ".bai")
        ):
            return

        print(
            "index file: {bai} doesn't exist, create the index file.\n".format(
                bai=self.sample_file + ".bai"
            ),
            flush=True,
        )
        pysam.sort("-o", self.output_prefix + ".temp.sorted.bam", self.sample_file, "-@", str(self.thread))
        pysam.index(self.output_prefix + ".temp.sorted.bam", "-@", str(self.thread))
        self.sample_file = self.output_prefix + ".temp.sorted.bam"
        self.sample_format = "rb"

    @staticmethod
    def _new_filter_stats():
        return Counter({
            "total_records_scanned": 0,
            "unmapped_reads": 0,
            "secondary_reads": 0,
            "supplementary_reads": 0,
            "duplicate_reads": 0,
            "qc_fail_reads": 0,
            "missing_length_reads": 0,
            "missing_unique_tag_reads": 0,
            "removed_multi_mapping_reads": 0,
            "removed_reverse_strand_reads": 0,
            "kept_plus_reads": 0,
            "kept_minus_reads": 0,
            "written_reads": 0,
            "fetch_missing_reference": 0,
        })

    @staticmethod
    def _passes_unique_filter(reads, align, tag, stats):
        """Unified unique-read filter.

        tag == 0:
            keep all mapped reads.
        tag == 1:
            prefer NH tag; for Bowtie/Bowtie2 fall back to XS; otherwise use MAPQ > 0.
        """
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
            # Bowtie2 usually reports XS:i for reads with a valid second-best alignment.
            # If XS exists, treat it as non-unique.
            if reads.has_tag("XS"):
                return False

        if reads.mapping_quality is not None and reads.mapping_quality > 0:
            return True

        stats["missing_unique_tag_reads"] += 1
        return False

    @staticmethod
    def _record_length(length_dict, read_length, is_reverse):
        if read_length not in length_dict:
            length_dict[read_length] = [0, 0]
        if is_reverse:
            length_dict[read_length][1] += 1
        else:
            length_dict[read_length][0] += 1

    @staticmethod
    def filter_alignments(args):
        """Filter BAM/SAM reads with a single unified implementation."""
        (
            bam_in_file,
            bam_out_file,
            gene_list,
            tag,
            reverse_flag,
            sample_format,
            align,
            saturation_flag,
        ) = args

        bam_seq_dict = {}
        length_dict = {}
        stats = Quality._new_filter_stats()

        with pysam.AlignmentFile(bam_in_file, sample_format) as bam_in:
            with pysam.AlignmentFile(bam_out_file, "wb", template=bam_in) as bam_out:
                for gene in gene_list:
                    try:
                        gene_reads = bam_in.fetch(gene)
                    except ValueError:
                        stats["fetch_missing_reference"] += 1
                        continue

                    for reads in gene_reads:
                        stats["total_records_scanned"] += 1

                        if reads.is_unmapped:
                            stats["unmapped_reads"] += 1
                            continue

                        if reads.is_secondary:
                            stats["secondary_reads"] += 1
                        if reads.is_supplementary:
                            stats["supplementary_reads"] += 1
                        if reads.is_duplicate:
                            stats["duplicate_reads"] += 1
                        if reads.is_qcfail:
                            stats["qc_fail_reads"] += 1

                        if not Quality._passes_unique_filter(reads, align, tag, stats):
                            stats["removed_multi_mapping_reads"] += 1
                            continue

                        read_length = reads.infer_read_length()
                        if read_length is None:
                            stats["missing_length_reads"] += 1
                            continue

                        # Keep legacy behavior: count minus-strand reads in length table,
                        # but only write them to output BAM when -r/--reverse is specified.
                        Quality._record_length(length_dict, read_length, reads.is_reverse)

                        if reads.is_reverse:
                            if not reverse_flag:
                                stats["removed_reverse_strand_reads"] += 1
                                continue
                            stats["kept_minus_reads"] += 1
                        else:
                            stats["kept_plus_reads"] += 1

                        bam_out.write(reads)
                        stats["written_reads"] += 1

                        if saturation_flag:
                            try:
                                bam_seq_dict[reads.query_name].add(reads.reference_name)
                            except KeyError:
                                bam_seq_dict[reads.query_name] = {reads.reference_name}

        return length_dict, bam_seq_dict, stats

    # Compatibility wrappers. External code that imports these names will still work.
    @staticmethod
    def flt_star_results(args):
        bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format = args
        return Quality.filter_alignments(
            (bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format, "star", True)
        )[:2]

    @staticmethod
    def flt_hisat2_results(args):
        bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format = args
        return Quality.filter_alignments(
            (bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format, "hisat2", True)
        )[:2]

    @staticmethod
    def flt_bowtie2_results(args):
        bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format = args
        return Quality.filter_alignments(
            (bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format, "bowtie2", True)
        )[:2]

    def fliter_mrna_reads(self):
        """Retrieve filtered RPF reads from BAM/SAM using a unified read filter."""
        pool = Pool(processes=self.thread)

        splits = np.array_split(list(self.mrna_dict.keys()), self.thread)
        splits = [list(split) for split in splits]

        args = [
            (
                self.sample_file,
                self.output_prefix + "_split_" + str(i) + ".bam",
                split,
                self.tag,
                self.reverse,
                self.sample_format,
                self.align,
                self.saturation_flag,
            )
            for i, split in enumerate(splits)
        ]

        try:
            results = pool.map(self.filter_alignments, args)
        finally:
            pool.close()
            pool.join()

        self.filter_stats = self._new_filter_stats()

        for length_dict, bam_seq_dict, stats in results:
            for key, value in length_dict.items():
                try:
                    self.length_dict[key][0] += value[0]
                    self.length_dict[key][1] += value[1]
                except KeyError:
                    self.length_dict[key] = value

            for key, value in stats.items():
                self.filter_stats[key] += value

            if self.saturation_flag:
                for key, value in bam_seq_dict.items():
                    self.bam_seq_dict.setdefault(key, set()).update(value)

        del results
        self.ensure_non_empty_results()

    def merge_sort_index_bam(self):
        """Merge split BAM files, sort safely using different input/output names, then index."""
        bam_in_list = [
            self.output_prefix + "_split_" + str(i) + ".bam"
            for i in range(self.thread)
        ]

        merged_bam = self.output_prefix + ".merged.unsorted.bam"
        sorted_bam = self.output_prefix + ".sorted.bam"
        final_bam = self.output_prefix + ".bam"

        merge_parameters = ["-f", "-@", str(self.thread), merged_bam] + bam_in_list
        pysam.merge(*merge_parameters)

        if os.path.exists(self.output_prefix + ".temp.sorted.bam"):
            os.remove(self.output_prefix + ".temp.sorted.bam")
        if os.path.exists(self.output_prefix + ".temp.sorted.bam.bai"):
            os.remove(self.output_prefix + ".temp.sorted.bam.bai")

        for bam_file in bam_in_list:
            if os.path.exists(bam_file):
                os.remove(bam_file)

        pysam.sort("-o", sorted_bam, merged_bam, "-@", str(self.thread))
        os.replace(sorted_bam, final_bam)

        if os.path.exists(merged_bam):
            os.remove(merged_bam)

        pysam.index(final_bam, "-@", str(self.thread))

    def ensure_non_empty_results(self):
        """Fail early with a clear error if no reads survived rpf_Check filtering."""
        total_length_records = sum(sum(v) for v in self.length_dict.values()) if self.length_dict else 0
        written_reads = int(self.filter_stats.get("written_reads", 0)) if self.filter_stats else 0

        if total_length_records == 0 or written_reads == 0:
            try:
                self.write_filter_stats()
                self.write_summary()
            except Exception:
                pass

            raise RuntimeError(
                "No reads remained after rpf_Check filtering. "
                "Please check: "
                "1) whether transcript_id in gene.norm.txt matches BAM reference names; "
                "2) whether -a/--align matches the aligner used to generate the BAM; "
                "3) whether -g 1 removed all reads because unique-mapping tags are missing; "
                "4) whether minus-strand reads require -r/--reverse; "
                "5) whether the input BAM/SAM is sorted and indexed."
            )

    def detect_seq_type(self):
        """Auto-detect the type of sequence profile from total plus+minus length counts."""
        if not self.length_dict:
            self.ensure_non_empty_results()

        self.peak_length, strand_counts = max(
            self.length_dict.items(),
            key=lambda item: sum(item[1])
        )
        self.peak_reads = int(sum(strand_counts))

        if 23 < self.peak_length < 35:
            print("{bam} is detected to be monosome-seq.\n".format(bam=self.sample_file), flush=True)
            self.profile = "monosome"
            self.mono = [19.5, 40.5]
        elif 53 < self.peak_length < 65:
            print("{bam} is detected to be disome-seq.\n".format(bam=self.sample_file), flush=True)
            self.profile = "disome"
            self.mono = [49.5, 70.5]
        elif 83 < self.peak_length < 95:
            print("{bam} is detected to be trisome-seq.\n".format(bam=self.sample_file), flush=True)
            self.profile = "trisome"
            self.mono = [79.5, 100.5]
        elif 35 < self.peak_length < 53 or 65 < self.peak_length < 83 or 95 < self.peak_length:
            print(
                "Warning! {bam} doesn't fit the empirical length distribution.!\n".format(
                    bam=self.sample_file
                ),
                flush=True,
            )
            print(
                """
Monosome RPFs peak length is usually ~ 30 nt.
Disome RPFs peak length is usually ~ 60 nt.
Trisome RPFs peak length is usually ~ 90 nt.
Please check the files and run the detect_offset.py with specified peak_length!""",
                flush=True,
            )
        else:
            print(
                "Warning! Cannot classify profile from peak length: {length} nt.\n".format(
                    length=self.peak_length
                ),
                flush=True,
            )

    def draw_the_length_distr(self, sorted_length):
        # output figure name
        out_pdf = self.output_prefix + "_length_distribution.pdf"
        out_png = self.output_prefix + "_length_distribution.png"

        matplotlib.use("AGG")

        # draw the figure
        length_df = pd.DataFrame(sorted_length).T
        fig = plt.figure(figsize=(6, 6), dpi=300)
        ax1 = fig.add_subplot(2, 1, 1)
        sns.lineplot(x=length_df.index, y=length_df[0], color="#FF9900")
        ax1.set_xticks(length_df.index)
        ax1.set_xticklabels(length_df.index)
        ax1.set_xlim(self.mono)
        ax1.set_ylabel("RPFs number")
        ax1.set_xlabel("RPFs length (nt)")
        ax1.set_title("plus strand")
        plt.xticks(rotation=90)

        # draw the minus strand RPFs length distribution
        ax2 = fig.add_subplot(2, 1, 2)
        sns.lineplot(x=length_df.index, y=length_df[1], color="#0099FF")
        ax2.set_xticks(length_df.index)
        ax2.set_xticklabels(length_df.index)
        ax2.set_xlim(self.mono)
        ax2.set_ylabel("RPFs number")
        ax2.set_xlabel("RPFs length (nt)")
        ax2.set_title("Minus strand")
        plt.xticks(rotation=90)
        plt.suptitle("RPFs length distribution")
        plt.tight_layout()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()

    def write_length_distr(self):
        sorted_length = dict(sorted(self.length_dict.items(), key=lambda length: length[0], reverse=False))
        with open(self.output_prefix + "_length_distribution.txt", "w") as length_out:
            length_out.writelines("\t".join(["Length", "Plus", "Minus"]) + "\n")
            for reads_length, reads_num in sorted_length.items():
                length_out.writelines(
                    "\t".join([str(reads_length), str(reads_num[0]), str(reads_num[1])]) + "\n"
                )

        # draw the RPFs length distribution
        self.draw_the_length_distr(sorted_length)

    def write_filter_stats(self):
        """Write filter statistics as a two-column TSV."""
        out_file = self.output_prefix + "_rpf_check.filter_stats.tsv"

        ordered_keys = [
            "total_records_scanned",
            "unmapped_reads",
            "secondary_reads",
            "supplementary_reads",
            "duplicate_reads",
            "qc_fail_reads",
            "missing_length_reads",
            "missing_unique_tag_reads",
            "removed_multi_mapping_reads",
            "removed_reverse_strand_reads",
            "kept_plus_reads",
            "kept_minus_reads",
            "written_reads",
            "fetch_missing_reference",
        ]

        with open(out_file, "w") as out:
            out.write("metric\tvalue\n")
            for key in ordered_keys:
                out.write(f"{key}\t{int(self.filter_stats.get(key, 0))}\n")

    def _summary_dict(self):
        plus_reads = int(sum(v[0] for v in self.length_dict.values())) if self.length_dict else 0
        minus_reads = int(sum(v[1] for v in self.length_dict.values())) if self.length_dict else 0

        return OrderedDict([
            ("sample", os.path.basename(self.sample_file)),
            ("input_bam_or_sam", self.sample_file),
            ("output_bam", self.output_prefix + ".bam"),
            ("output_prefix", self.output_prefix),
            ("aligner", self.align),
            ("unique_tag_mode", self.tag),
            ("reverse_counted", self.reverse),
            ("longest_enabled", self.longest),
            ("thread", self.thread),
            ("saturation_enabled", self.saturation_flag),
            ("saturation_seed", self.saturation_seed),
            ("transcripts_input", int(self.transcript_total_input)),
            ("transcripts_kept", int(self.transcript_total_kept)),
            ("transcripts_removed_by_longest", int(self.transcript_removed_by_longest)),
            ("total_records_scanned", int(self.filter_stats.get("total_records_scanned", 0))),
            ("written_reads", int(self.filter_stats.get("written_reads", 0))),
            ("plus_length_counted_reads", plus_reads),
            ("minus_length_counted_reads", minus_reads),
            ("removed_multi_mapping_reads", int(self.filter_stats.get("removed_multi_mapping_reads", 0))),
            ("removed_reverse_strand_reads", int(self.filter_stats.get("removed_reverse_strand_reads", 0))),
            ("missing_unique_tag_reads", int(self.filter_stats.get("missing_unique_tag_reads", 0))),
            ("fetch_missing_reference", int(self.filter_stats.get("fetch_missing_reference", 0))),
            ("peak_length", int(self.peak_length) if self.peak_length else 0),
            ("peak_reads", int(self.peak_reads) if self.peak_reads else 0),
            ("detected_profile", self.profile if self.profile else "undetermined"),
            ("covered_reads_for_saturation", len(self.bam_seq_dict) if self.saturation_flag else "NA"),
        ])

    def write_summary(self):
        """Write rpf_Check one-row summary as TSV and JSON."""
        summary = self._summary_dict()

        out_tsv = self.output_prefix + "_rpf_check.summary.tsv"
        out_json = self.output_prefix + "_rpf_check.summary.json"

        with open(out_tsv, "w") as out:
            out.write("\t".join(summary.keys()) + "\n")
            out.write("\t".join(str(v) for v in summary.values()) + "\n")

        with open(out_json, "w") as out:
            json.dump(summary, out, indent=2)
            out.write("\n")

    @staticmethod
    def sample_keys(step, keys):
        keys_list = list(keys)
        return random.sample(keys_list, step)

    def rpf_saturation(self):
        """Calculate RPF saturation with deterministic and faster incremental sampling.

        Algorithm:
        1. Use fixed seed 5201314.
        2. Shuffle read names once.
        3. Use nested prefixes at 10%, 20%, ..., 90%.
        4. Incrementally update one Counter instead of rebuilding counters 9 times.
        """
        mapped_reads = len(self.bam_seq_dict)
        if mapped_reads == 0:
            raise RuntimeError(
                "Saturation was requested, but no reads are available for saturation. "
                "Please check whether reads survived filtering."
            )

        rng = random.Random(self.saturation_seed)
        read_names = list(self.bam_seq_dict.keys())
        rng.shuffle(read_names)

        steps = [int(i * 0.10 * mapped_reads) for i in range(1, 10)]
        steps = [max(1, min(step, mapped_reads)) for step in steps]

        running_gene_counter = Counter()
        cursor = 0
        self.saturation = []

        for site, step in enumerate(steps):
            for read_name in read_names[cursor:step]:
                running_gene_counter.update(self.bam_seq_dict[read_name])
            cursor = step

            self.saturation.append(len(running_gene_counter))

            for gene, num in running_gene_counter.items():
                if gene in self.mrna_dict:
                    self.mrna_dict[gene][site] = num

    @staticmethod
    def process_step(step_data):
        """Compatibility helper for older external code."""
        step, bam_seq_dict, seed = step_data
        rng = random.Random(seed)
        keys = list(bam_seq_dict.keys())
        temp = rng.sample(keys, step)
        temp_gene = []
        for i in temp:
            temp_gene.extend(bam_seq_dict[i])
        temp_gene_dict = dict(Counter(temp_gene))
        gene_num = len(temp_gene_dict)
        result = (gene_num, temp_gene_dict)
        return result

    def rpf_saturation_thread(self):
        """Compatibility alias for the optimized deterministic saturation method."""
        return self.rpf_saturation()

    def draw_gene_saturation(self):
        out_pdf = self.output_prefix + "_gene_saturation.pdf"
        out_png = self.output_prefix + "_gene_saturation.png"
        out_gene = self.output_prefix + "_gene_saturation.txt"

        total_gene_num = len(self.mrna_dict)

        gene_df = pd.DataFrame(self.saturation + [total_gene_num], columns=["Count"])
        gene_df["Part"] = self.x_ticks.tolist() + [0]
        gene_df = gene_df[["Part", "Count"]]
        gene_df.to_csv(out_gene, sep="\t", index=False)

        matplotlib.use("AGG")
        fig = plt.figure(figsize=(8, 4), dpi=300)
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.bar(self.x_ticks, self.saturation, width=3, color="#FF9900")
        ax1.plot(self.x_ticks, self.saturation, linewidth=0.8)
        plt.xticks(rotation=90)
        ax1.set_xlabel("reads proportion (%)")
        ax1.set_ylabel("gene number")
        plt.title("covered gene saturation")

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.bar(self.x_ticks, [total_gene_num - i for i in self.saturation], width=3, color="#FF9900")
        ax2.plot(self.x_ticks, [total_gene_num - i for i in self.saturation], linewidth=0.8)
        plt.xticks(rotation=90)
        ax2.set_xlabel("reads proportion (%)")
        ax2.set_ylabel("gene number")
        plt.title("uncovered gene saturation")

        plt.tight_layout()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()

    def draw_rpf_saturation(self):
        # output figure name
        out_pdf = self.output_prefix + "_reads_saturation.pdf"
        out_png = self.output_prefix + "_reads_saturation.png"
        out_rpf = self.output_prefix + "_reads_saturation.txt"

        mrna_df = pd.DataFrame(self.mrna_dict).T
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

        # draw the figure
        matplotlib.use("AGG")
        flierprops = dict(
            marker="o",
            markersize=2,
            markerfacecolor="#c95859",
            markeredgecolor="none",
        )

        fig = plt.figure(figsize=(17, 4), dpi=300)

        fig.add_subplot(1, 4, 1)
        ax1 = sns.boxplot(data=mrna_0_25.iloc[:, 0:9], color="#30adfe", flierprops=flierprops)
        ax1.set(xlabel="reads proportion (%)", ylabel="Reads count")
        ax1.set_yscale("log")
        ax1.set_xticklabels(self.x_ticks, rotation=90)

        fig.add_subplot(1, 4, 2)
        ax2 = sns.boxplot(data=mrna_25_50.iloc[:, 0:9], color="#30adfe", flierprops=flierprops)
        ax2.set(xlabel="reads proportion (%)", ylabel="Reads count")
        ax2.set_yscale("log")
        ax2.set_xticklabels(self.x_ticks, rotation=90)

        fig.add_subplot(1, 4, 3)
        ax3 = sns.boxplot(data=mrna_50_75.iloc[:, 0:9], color="#30adfe", flierprops=flierprops)
        ax3.set(xlabel="reads proportion (%)", ylabel="Reads count")
        ax3.set_yscale("log")
        ax3.set_xticklabels(self.x_ticks, rotation=90)

        fig.add_subplot(1, 4, 4)
        ax4 = sns.boxplot(data=mrna_75_100.iloc[:, 0:9], color="#30adfe", flierprops=flierprops)
        ax4.set(xlabel="reads proportion (%)", ylabel="Reads count")
        ax4.set_yscale("log")
        ax4.set_xticklabels(self.x_ticks, rotation=90)

        plt.suptitle("Reads saturation")
        plt.tight_layout()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()
