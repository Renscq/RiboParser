#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-08
# Version: dev004
# Function: Run legacy or replicate-consensus SeRP peak calling with unified density input.
# Input: RPF density data, control/IP sample definitions, normalization totals, and annotation.
# Output: Peak tables, BED files, peak sequences, and per-codon enrichment profiles.

"""Workflow controller for selective ribosome profiling peak calling."""

from __future__ import annotations

from collections import OrderedDict

import pandas as pd

from utils.ribo.RPFs import BASE_COLUMNS, RPFData
from utils.serp.SeRPOutput import SeRPOutput
from utils.serp import Peak as peak_module

PeakCaller = peak_module.PeakCaller


class SeRP(PeakCaller):
    """Coordinate SeRP density import, peak calling, and result export.

    Notes
    -----
    Peak-calling algorithms are implemented in :class:`utils.serp.Peak.PeakCaller`.
    Figure generation is intentionally excluded from this class and handled by
    the independent ``serp_plot`` workflow.
    """

    def __init__(self, args):
        """Initialize the SeRP peak-calling workflow from parsed arguments."""
        self.rpf_file = args.rpf
        self.total_rpf_file = args.norm
        self.annotation_file = args.anno

        self.control_samples = self._parse_sample_names(args.control)
        self.ip_samples = self._parse_sample_names(args.ip)
        self.sample_names = self.control_samples + self.ip_samples

        self.background = args.background
        self.scale = args.scale
        self.fill = args.fill
        self.correlation = args.corr
        self.threshold = args.min

        self.size = args.size
        self.polyorder = args.k
        self.back_fold = args.back_fold
        self.enrich = args.enrich
        self.gaps = args.gaps
        self.proportion = args.proportion
        self.width = args.width
        self.collision = args.collision
        self.keep_all = args.keep_all

        self.method = args.method
        self._validate_peak_caller_api()
        self.consensus_window = args.consensus_window
        self.pseudocount = args.pseudocount
        self.min_support = args.min_support
        self.min_overlap = args.min_overlap
        self.stop_trim = args.stop_trim
        self.min_codon_rpf = args.min_codon_rpf
        self.max_edge_extension = args.max_edge_extension

        self.output = SeRPOutput(
            output_prefix=args.output,
            sample_names=self.sample_names,
            upstream=args.upstream,
            downstream=args.downstream,
            ratio_out=args.ratio,
        )

        self.gene_list: list[str] = []
        self.gene_dict: OrderedDict[str, str] = OrderedDict()
        self.all_rpf = pd.DataFrame()
        self.rpf_seq = pd.DataFrame()
        self.total_rpf_num = pd.Series(dtype="float64")
        self.mean_aa_rpm = pd.Series(dtype="float64")
        self.peak_merge: list[list[str]] | pd.DataFrame = []

    def _validate_peak_caller_api(self) -> None:
        """Validate that SeRP.py and Peak.py belong to the same caller API.

        Raises
        ------
        RuntimeError
            If the installed ``utils.serp.Peak`` module is too old for the
            selected caller. This converts otherwise late ``AttributeError``
            failures into an immediate, actionable package-version error.
        """
        required_methods = [
            "filter_background_max",
            "fill_control_zeros",
            "calculate_enrichment_ratio",
            "smooth_enrichment",
            "call_peak_regions",
            "extend_collision_regions",
            "summarize_peak_enrichment",
            "ttest_peaks",
        ]
        if self.method == "consensus":
            required_methods.extend(
                [
                    "calculate_consensus_enrichment",
                    "consensus_analysis_positions",
                    "call_consensus_peak_regions",
                    "extend_peak_edges",
                    "consensus_peak_statistics",
                ]
            )

        missing_methods = [
            method_name
            for method_name in required_methods
            if not hasattr(PeakCaller, method_name)
        ]
        if not missing_methods:
            return

        api_version = getattr(peak_module, "PEAK_CALLER_API_VERSION", "unknown")
        raise RuntimeError(
            "Incompatible SeRP/Peak installation: Peak caller API {api} is "
            "missing required method(s): {methods}. Replace SeRP.py and Peak.py "
            "together from the same RiboParser update, then reinstall the package."
            .format(
                api=api_version,
                methods=", ".join(missing_methods),
            )
        )

    @staticmethod
    def _parse_sample_names(value: str) -> list[str]:
        """Parse a comma-separated sample list while preserving order."""
        return [item.strip() for item in value.split(",") if item.strip()]

    @staticmethod
    def _safe_to_numeric(series: pd.Series) -> pd.Series:
        """Convert a Series to numeric values when conversion is possible."""
        try:
            return pd.to_numeric(series)
        except (TypeError, ValueError):
            return series

    def prepare_output(self) -> None:
        """Prepare output files for a new peak-calling run."""
        self.output.prepare()

    def import_rpf(self) -> None:
        """Import TXT or JSON RPF density data through the common RPF reader.

        The three reading-frame columns for each requested sample are merged so
        downstream SeRP calculations receive the same sample-level codon table
        regardless of whether the source density file is legacy TXT or current
        JSON/JSONL format.
        """
        rpf_data = RPFData.from_file(
            rpf_file=self.rpf_file,
            sample_name=self.sample_names,
        )
        merged_rpf = rpf_data.get_frame(frame="all")

        missing_samples = [
            sample for sample in self.sample_names if sample not in merged_rpf.columns
        ]
        if missing_samples:
            raise ValueError(
                "RPF density data are missing requested sample(s): {samples}".format(
                    samples=", ".join(missing_samples)
                )
            )

        self.rpf_seq = merged_rpf.loc[:, BASE_COLUMNS].copy()
        self.all_rpf = merged_rpf.loc[
            :, ["name", "from_tis", "region"] + self.sample_names
        ].copy()
        self.gene_list = self.all_rpf["name"].drop_duplicates().astype(str).tolist()

        if self.total_rpf_file:
            total_table = pd.read_csv(
                self.total_rpf_file,
                sep="\t",
                header=None,
                usecols=[0, 1],
                names=["sample", "total_rpf"],
            )
            total_series = pd.Series(
                pd.to_numeric(total_table["total_rpf"], errors="raise").values,
                index=total_table["sample"].astype(str),
                dtype="float64",
            )
            missing_totals = [
                sample for sample in self.sample_names if sample not in total_series.index
            ]
            if missing_totals:
                raise ValueError(
                    "Normalization file is missing sample(s): {samples}".format(
                        samples=", ".join(missing_totals)
                    )
                )
            self.total_rpf_num = total_series.reindex(self.sample_names)
        else:
            self.total_rpf_num = rpf_data.total_rpf_num.reindex(self.sample_names)

        if self.total_rpf_num.isna().any():
            missing_totals = self.total_rpf_num[self.total_rpf_num.isna()].index.tolist()
            raise ValueError(
                "Cannot determine total RPF counts for sample(s): {samples}".format(
                    samples=", ".join(missing_totals)
                )
            )
        if (self.total_rpf_num <= 0).any():
            invalid_samples = self.total_rpf_num[self.total_rpf_num <= 0].index.tolist()
            raise ValueError(
                "Total RPF counts must be positive for sample(s): {samples}".format(
                    samples=", ".join(invalid_samples)
                )
            )

        print("Input density format: {fmt}".format(fmt=rpf_data.file_format), flush=True)
        print("Transcripts: {num:,}".format(num=len(self.gene_list)), flush=True)
        print("Total RPFs number:", flush=True)
        print(self.total_rpf_num, flush=True)

        self._calculate_global_background_rpm()

    def _calculate_global_background_rpm(self) -> None:
        """Calculate the global control RPM value used by fill mode 0.

        Notes
        -----
        The legacy implementation divided by ``background`` directly and
        therefore produced division-by-zero values when ``--back 0``. For that
        mode, the denominator is now the mean CDS length in codons, which keeps
        the quantity interpretable as mean control RPM per CDS codon. The
        historical ``--back 30`` behavior is otherwise retained.
        """
        cds_rpf = self.all_rpf.loc[
            self.all_rpf["region"] == "cds", self.control_samples
        ].astype(float)
        if cds_rpf.empty:
            self.mean_aa_rpm = pd.Series(0.0, index=self.control_samples)
            return

        cds_rpm = cds_rpf.div(
            self.total_rpf_num.loc[self.control_samples].values,
            axis=1,
        ) * self.scale

        if self.background > 0:
            denominator = float(self.background * max(len(self.gene_list), 1))
            self.mean_aa_rpm = cds_rpm.sum(axis=0) / denominator
        else:
            mean_cds_length = float(len(cds_rpm)) / float(max(len(self.gene_list), 1))
            denominator = max(mean_cds_length * max(len(self.gene_list), 1), 1.0)
            self.mean_aa_rpm = cds_rpm.sum(axis=0) / denominator

    def import_annotation(self) -> None:
        """Import transcript-to-gene names or assign ``-`` when unavailable."""
        self.gene_dict = OrderedDict()

        if self.annotation_file:
            with open(self.annotation_file, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    self.gene_dict[str(fields[2])] = str(fields[1])

        for transcript in self.gene_list:
            if transcript not in self.gene_dict:
                self.gene_dict[transcript] = "-"

    def _calculate_gene_rpf(
        self,
        gene_rpf: pd.DataFrame,
    ) -> tuple[list[float], pd.DataFrame, list[float]]:
        """Calculate raw gene counts and RPM values for all SeRP samples."""
        gene_rpf_sum = gene_rpf[self.sample_names].sum().tolist()
        raw_gene_rpm = gene_rpf[self.sample_names].div(
            self.total_rpf_num.loc[self.sample_names].values,
            axis=1,
        ) * self.scale
        raw_gene_rpm_sum = raw_gene_rpm[self.sample_names].sum().tolist()
        return gene_rpf_sum, raw_gene_rpm, raw_gene_rpm_sum

    def _append_peak_rows(
        self,
        transcript: str,
        gene_rpf_sum: list[float],
        raw_gene_rpm_sum: list[float],
        ck_min_corr: str,
        ip_min_corr: str,
        peak: OrderedDict,
        peak_left_start: list[int],
        peak_right_end: list[int],
        max_peak_enrich: dict,
        max_peak_site: dict,
        mean_peak_enrich: dict,
        ck_rpm: dict,
        ip_rpm: dict,
        ttest_results: dict,
        peak_results,
        utr5_region: list[int],
        cds_region: list[int],
        utr3_region: list[int],
        gap_counts: dict[int, int] | None = None,
    ) -> None:
        """Format and append called peaks to the current result collection."""
        for peak_num, peak_dict in peak.items():
            peak_len = list(peak_dict.keys())[0]
            peak_start = peak[peak_num][peak_len][0]
            peak_end = peak[peak_num][peak_len][-1]
            if gap_counts is None:
                gaps = (peak_end - peak_start + 1) - peak_len
            else:
                gaps = int(gap_counts.get(peak_num, 0))
            peak_loci = self.annotate_peak_region(
                utr5_region,
                cds_region,
                utr3_region,
                peak_start,
                peak_end,
            )

            row = self.output.build_peak_row(
                transcript=transcript,
                gene_name=self.gene_dict[transcript],
                gene_rpf_sum=gene_rpf_sum,
                raw_gene_rpm_sum=raw_gene_rpm_sum,
                ck_min_corr=ck_min_corr,
                ip_min_corr=ip_min_corr,
                peak_num=peak_num,
                peak_len=peak_len,
                gaps=gaps,
                collision_start=peak_left_start[peak_num],
                peak_start=peak_start,
                peak_end=peak_end,
                collision_end=peak_right_end[peak_num],
                max_site=max_peak_site[peak_num],
                max_fold=max_peak_enrich[peak_num],
                mean_fold=mean_peak_enrich[peak_num],
                ck_peak_rpm=ck_rpm[peak_num],
                ip_peak_rpm=ip_rpm[peak_num],
                pvalue=ttest_results[peak_num],
                peak_loci=peak_loci,
            )
            self.peak_merge.append(row)
            self.output.write_log_row(peak_results, row)
        peak_results.flush()

    def _replicate_correlation(
        self,
        raw_gene_rpm: pd.DataFrame,
        positions: list[int] | None = None,
    ) -> tuple[str, str, bool, bool]:
        """Calculate group-wise replicate correlations and pass/fail flags."""
        if positions is None:
            profile = raw_gene_rpm
        else:
            profile = raw_gene_rpm.loc[positions]

        ck_corr = profile[self.control_samples].corr()
        ip_corr = profile[self.ip_samples].corr()
        ck_min_corr = str(round(float(ck_corr.values.min()), 4))
        ip_min_corr = str(round(float(ip_corr.values.min()), 4))

        ck_corr_pairs, ck_total_pairs = self.count_correlation_pairs(
            ck_corr,
            self.control_samples,
        )
        ip_corr_pairs, ip_total_pairs = self.count_correlation_pairs(
            ip_corr,
            self.ip_samples,
        )
        return (
            ck_min_corr,
            ip_min_corr,
            ck_corr_pairs >= ck_total_pairs,
            ip_corr_pairs >= ip_total_pairs,
        )

    def _append_no_peak_row(
        self,
        peak_results,
        transcript: str,
        gene_rpf_sum: list[float],
        raw_gene_rpm_sum: list[float],
        ck_min_corr: str,
        ip_min_corr: str,
        comment: str = "Qualified",
    ) -> None:
        """Append one transcript-level status row without a called peak."""
        row = self.output.build_status_row(
            transcript,
            self.gene_dict[transcript],
            gene_rpf_sum,
            raw_gene_rpm_sum,
            ck_min_corr,
            ip_min_corr,
            comment,
            peak_num="0" if comment == "Qualified" else "-",
        )
        self.peak_merge.append(row)
        self.output.write_log_row(peak_results, row)

    def _call_legacy_transcript(
        self,
        peak_results,
        transcript: str,
        gene_rpf: pd.DataFrame,
        gene_rpf_sum: list[float],
        raw_gene_rpm: pd.DataFrame,
        raw_gene_rpm_sum: list[float],
        utr5_region: list[int],
        cds_region: list[int],
        utr3_region: list[int],
    ) -> None:
        """Run the corrected historical RiboParser SeRP peak-calling workflow."""
        background_pass = self.filter_background_max(
            self.ip_samples,
            raw_gene_rpm,
            cds_region,
        )

        gene_rpm = self.fill_control_zeros(
            self.control_samples,
            raw_gene_rpm,
            cds_region,
        )
        gene_ratio = self.calculate_enrichment_ratio(
            gene_rpm,
            self.control_samples,
            self.ip_samples,
        )
        _, gene_ratio_smooth = self.smooth_enrichment(gene_ratio, gene_rpf)
        self.output.write_raw_ratio(transcript, gene_ratio)

        if not background_pass:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                "-",
                "-",
                "Background value too large",
            )
            self.output.write_profile(gene_rpf, gene_ratio_smooth)
            return

        ck_min_corr, ip_min_corr, ck_pass, ip_pass = self._replicate_correlation(
            raw_gene_rpm
        )
        if not ck_pass:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
                "Poor repeatability of ck samples",
            )
            self.output.write_profile(gene_rpf, gene_ratio_smooth)
            return
        if not ip_pass:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
                "Poor repeatability of ip samples",
            )
            self.output.write_profile(gene_rpf, gene_ratio_smooth)
            return

        eligible_index = gene_ratio_smooth.loc[
            gene_ratio_smooth["enrich"] >= self.enrich
        ].index.tolist()
        if len(eligible_index) < self.width:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
            )
            self.output.write_profile(gene_rpf, gene_ratio_smooth)
            return

        peak = self.call_peak_regions(
            eligible_index,
            gene_ratio_smooth,
            raw_gene_rpm,
            self.control_samples,
            self.ip_samples,
        )
        if not peak:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
            )
            self.output.write_profile(gene_rpf, gene_ratio_smooth)
            return

        peak_left_start, peak_right_end = self.extend_collision_regions(
            gene_ratio_smooth,
            peak,
        )
        max_peak_enrich, max_peak_site, mean_peak_enrich = (
            self.summarize_peak_enrichment(peak, gene_ratio_smooth)
        )
        ttest_results, ck_rpm, ip_rpm, _ = self.ttest_peaks(
            peak,
            raw_gene_rpm,
            self.control_samples,
            self.ip_samples,
            transcript,
        )

        self._append_peak_rows(
            transcript,
            gene_rpf_sum,
            raw_gene_rpm_sum,
            ck_min_corr,
            ip_min_corr,
            peak,
            peak_left_start,
            peak_right_end,
            max_peak_enrich,
            max_peak_site,
            mean_peak_enrich,
            ck_rpm,
            ip_rpm,
            ttest_results,
            peak_results,
            utr5_region,
            cds_region,
            utr3_region,
        )
        self.output.write_profile(
            gene_rpf,
            gene_ratio_smooth,
            peak=peak,
            peak_left_start=peak_left_start,
            peak_right_end=peak_right_end,
        )

    def _call_consensus_transcript(
        self,
        peak_results,
        transcript: str,
        gene_rpf: pd.DataFrame,
        gene_rpf_sum: list[float],
        raw_gene_rpm: pd.DataFrame,
        raw_gene_rpm_sum: list[float],
        utr5_region: list[int],
        cds_region: list[int],
        utr3_region: list[int],
    ) -> None:
        """Run replicate-consensus SeRP peak calling for one transcript."""
        replicate_ratio, consensus_profile = self.calculate_consensus_enrichment(
            raw_gene_rpm,
            self.control_samples,
            self.ip_samples,
        )
        self.output.write_raw_ratio(transcript, replicate_ratio)

        analysis_positions = self.consensus_analysis_positions(cds_region)
        if len(analysis_positions) < self.width:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                "-",
                "-",
            )
            self.output.write_profile(gene_rpf, consensus_profile)
            return

        if self.min_codon_rpf > 0:
            mean_codon_rpf = gene_rpf.loc[analysis_positions, self.sample_names].mean(axis=0)
            if bool((mean_codon_rpf < self.min_codon_rpf).any()):
                self._append_no_peak_row(
                    peak_results,
                    transcript,
                    gene_rpf_sum,
                    raw_gene_rpm_sum,
                    "-",
                    "-",
                    "Low CDS coverage",
                )
                self.output.write_profile(gene_rpf, consensus_profile)
                return

        ck_min_corr, ip_min_corr, ck_pass, ip_pass = self._replicate_correlation(
            raw_gene_rpm,
            positions=analysis_positions,
        )
        if not ck_pass:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
                "Poor repeatability of ck samples",
            )
            self.output.write_profile(gene_rpf, consensus_profile)
            return
        if not ip_pass:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
                "Poor repeatability of ip samples",
            )
            self.output.write_profile(gene_rpf, consensus_profile)
            return

        peak, gap_counts = self.call_consensus_peak_regions(
            replicate_ratio,
            consensus_profile,
            cds_region,
        )
        if not peak:
            self._append_no_peak_row(
                peak_results,
                transcript,
                gene_rpf_sum,
                raw_gene_rpm_sum,
                ck_min_corr,
                ip_min_corr,
            )
            self.output.write_profile(gene_rpf, consensus_profile)
            return

        peak_left_start, peak_right_end = self.extend_peak_edges(
            consensus_profile,
            peak,
            max_extension=self.max_edge_extension,
            allowed_positions=analysis_positions,
        )
        max_peak_enrich, max_peak_site, mean_peak_enrich = (
            self.summarize_peak_enrichment(peak, consensus_profile)
        )
        pvalues, ck_rpm, ip_rpm, _ = self.consensus_peak_statistics(
            peak,
            raw_gene_rpm,
            self.control_samples,
            self.ip_samples,
            transcript,
        )

        self._append_peak_rows(
            transcript,
            gene_rpf_sum,
            raw_gene_rpm_sum,
            ck_min_corr,
            ip_min_corr,
            peak,
            peak_left_start,
            peak_right_end,
            max_peak_enrich,
            max_peak_site,
            mean_peak_enrich,
            ck_rpm,
            ip_rpm,
            pvalues,
            peak_results,
            utr5_region,
            cds_region,
            utr3_region,
            gap_counts=gap_counts,
        )
        self.output.write_profile(
            gene_rpf,
            consensus_profile,
            peak=peak,
            peak_left_start=peak_left_start,
            peak_right_end=peak_right_end,
        )

    def call_peaks(self) -> None:
        """Call SeRP peaks transcript by transcript without drawing figures."""
        columns = self.output.result_columns
        sample_num = len(self.sample_names)

        with self.output.open_peak_log() as peak_results:
            for transcript, gene_rpf in self.all_rpf.groupby("name", sort=False):
                transcript = str(transcript)
                print(transcript, flush=True)

                gene_rpf = gene_rpf.copy().apply(self._safe_to_numeric)
                gene_rpf.index = gene_rpf["from_tis"].astype(int)

                gene_rpf_sum, raw_gene_rpm, raw_gene_rpm_sum = self._calculate_gene_rpf(
                    gene_rpf
                )
                if sum(value > self.threshold for value in gene_rpf_sum) < sample_num:
                    row = self.output.build_status_row(
                        transcript,
                        self.gene_dict[transcript],
                        gene_rpf_sum,
                        raw_gene_rpm_sum,
                        "-",
                        "-",
                        "Too few RPFs",
                    )
                    self.peak_merge.append(row)
                    self.output.write_log_row(peak_results, row)
                    continue

                utr5_region = gene_rpf.loc[gene_rpf["region"] == "5utr"].index.tolist()
                cds_region = gene_rpf.loc[gene_rpf["region"] == "cds"].index.tolist()
                utr3_region = gene_rpf.loc[gene_rpf["region"] == "3utr"].index.tolist()
                if not cds_region:
                    row = self.output.build_status_row(
                        transcript,
                        self.gene_dict[transcript],
                        gene_rpf_sum,
                        raw_gene_rpm_sum,
                        "-",
                        "-",
                        "No CDS region",
                    )
                    self.peak_merge.append(row)
                    self.output.write_log_row(peak_results, row)
                    continue

                if self.method == "consensus":
                    self._call_consensus_transcript(
                        peak_results,
                        transcript,
                        gene_rpf,
                        gene_rpf_sum,
                        raw_gene_rpm,
                        raw_gene_rpm_sum,
                        utr5_region,
                        cds_region,
                        utr3_region,
                    )
                else:
                    self._call_legacy_transcript(
                        peak_results,
                        transcript,
                        gene_rpf,
                        gene_rpf_sum,
                        raw_gene_rpm,
                        raw_gene_rpm_sum,
                        utr5_region,
                        cds_region,
                        utr3_region,
                    )

        self.peak_merge = pd.DataFrame(self.peak_merge, columns=columns)
        print("Transcripts number: {num:,}".format(num=len(self.gene_list)), flush=True)

    def write_results(self) -> None:
        """Write final SeRP result tables through the dedicated output layer."""
        self.peak_merge = self.output.write_results(
            peak_merge=self.peak_merge,
            rpf_seq=self.rpf_seq,
            gene_dict=self.gene_dict,
        )