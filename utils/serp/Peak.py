#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-08
# Version: dev003
# Function: Provide legacy and replicate-consensus SeRP peak-calling algorithms.
# Input: Gene-level RPF/RPM profiles and peak-calling parameters from SeRP workflow.
# Output: Enrichment profiles, consensus peaks, edge regions, and peak statistics.

"""Peak-calling algorithms for selective ribosome profiling analysis."""

from __future__ import annotations

import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import ranksums, ttest_rel

from utils.serp.Enrichment import calculate_matched_local_enrichment


class PeakCaller:
    """Provide legacy and replicate-consensus SeRP peak-calling algorithms.

    Notes
    -----
    The ``legacy`` methods preserve the historical RiboParser implementation
    for reproducibility. The ``consensus`` methods calculate enrichment for
    matched biological replicate pairs, call regions independently per pair,
    and retain regions supported by replicate-level peak overlap.
    """
    def filter_background_max(self, ip_name, raw_gene_rpm, cds_region):
        """Filter genes using the maximum background signal."""
        # condition 1: bound to ribosome
        if self.background == 0:
            return True
        # condition 2: bound to nascent polypeptides
        elif self.background == 30:

            max_back_value = raw_gene_rpm.loc[0:self.background - 1].max()
            # Exclude outliers of stop codon, "cds_region[-1]-1"
            max_bound_value = raw_gene_rpm.loc[self.background:cds_region[-1] - 1].max()

            # filter the strong binding score
            if self.back_fold:
                delta_value = max_bound_value - max_back_value * self.enrich
            else:
                delta_value = max_bound_value - max_back_value * self.collision

            # flag-ip samples after background need 2 fold change higher than mock-ip
            if any(delta_value[ip_name] < 0):
                return False
            else:
                return True

    def filter_background_mean(self, ip_name, raw_gene_rpm, cds_region):
        """Filter genes using the mean background signal."""
        # condition 1: bound to ribosome
        if self.background == 0:
            return True

        # condition 2: bound to nascent polypeptides
        elif self.background == 30:

            mean_back_value = raw_gene_rpm.loc[0:self.background - 1].mean()
            # Exclude outliers of stop codon, "cds_region[-1]-1"
            mean_bound_value = raw_gene_rpm.loc[self.background:cds_region[-1] - 1].mean()

            # filter the strong binding score
            if self.back_fold:
                delta_value = mean_bound_value - mean_back_value * self.enrich
            else:
                delta_value = mean_bound_value - mean_back_value * self.collision

            # flag-ip samples after background need 2 fold change higher than mock-ip
            if any(delta_value[ip_name] < 0):
                return False
            else:
                return True

    def count_correlation_pairs(self, sp_gene_rpf_corr, sp_name):
        """Count replicate pairs passing the correlation threshold."""
        sp_corr_pairs = 0
        # more than one sample
        if len(sp_name) > 1:
            for rows in range(len(sp_name) - 1):
                for cols in range(rows, len(sp_name) - 1):
                    if sp_gene_rpf_corr.iloc[rows, cols + 1] >= self.correlation:
                        sp_corr_pairs += 1
        # only one sample in the group
        else:
            sp_corr_pairs = 1

        total_pairs = 0
        for counter in range(len(sp_name)):
            total_pairs += counter

        return sp_corr_pairs, total_pairs

    def fill_control_zeros(self, ck_name, raw_gene_rpm, cds_region):
        """Replace zero control RPM values using the configured background rule."""
        gene_rpm = raw_gene_rpm.copy()
        gene_rpm_mean = gene_rpm.loc[cds_region].mean()

        for sample in ck_name:
            # condition 1: fill the missing values with mean RPFs of total genes
            if self.fill == 0:
                # gene_rpm[sample].loc[gene_rpm[sample] == 0] = self.mean_aa_rpm[sample]
                gene_rpm.loc[gene_rpm[sample] == 0, sample] = self.mean_aa_rpm[sample]

            # condition 2: fill the missing values with mean RPFs of currently gene cds
            elif self.fill == -1:
                # gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_rpm_mean[sample]
                gene_rpm.loc[gene_rpm[sample] == 0, sample] = gene_rpm_mean[sample]

            # condition 3: fill the missing values with mean RPFs of specific region
            else:
                # Get the background mean value of the specific region
                gene_background = gene_rpm[sample].loc[0:self.fill - 1]
                counter_aa = len(gene_background[gene_background == 0].index)
                gene_background_mean = gene_background.mean()
                # If two-thirds of the specific region is covered by RPF.
                if counter_aa >= self.fill * 2 / 3:
                    # gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_background_mean
                    gene_rpm.loc[gene_rpm[sample] == 0, sample] = gene_background_mean

                # If half of the specific region is covered by RPF, and gene_rpm_mean <= gene_background_mean
                elif counter_aa >= self.fill / 2 and gene_rpm_mean[sample] <= gene_background_mean:
                    # gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_background_mean
                    gene_rpm.loc[gene_rpm[sample] == 0, sample] = gene_background_mean

                # if the specific region is noisy, use gene_rpm_mean instead
                elif counter_aa >= self.fill / 2 and gene_rpm_mean[sample] > gene_background_mean:
                    # gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_rpm_mean[sample]
                    gene_rpm.loc[gene_rpm[sample] == 0, sample] = gene_rpm_mean[sample]

                # if the specific region is noisy, use gene_rpm_mean instead
                elif counter_aa < self.fill / 2:
                    # gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_rpm_mean[sample]
                    gene_rpm.loc[gene_rpm[sample] == 0, sample] = gene_rpm_mean[sample]

                else:
                    # gene_rpm[sample].loc[gene_rpm[sample] == 0] = self.mean_aa_rpm[sample]
                    gene_rpm.loc[gene_rpm[sample] == 0, sample] = self.mean_aa_rpm[sample]

        return gene_rpm

    @staticmethod
    def calculate_enrichment_ratio(gene_rpm, ck_name, ip_name):
        """Calculate all IP-to-control enrichment ratios."""
        gene_ratio = pd.DataFrame()
        # Ratio calculation is performed on any two samples in the two grouped data.
        for ip_sample in ip_name:
            for ck_sample in ck_name:
                gene_ratio[ip_sample + '_' + ck_sample] = gene_rpm[ip_sample].div(gene_rpm[ck_sample])

        return gene_ratio

    def smooth_enrichment(self, gene_ratio, gene_rpf):
        """Average and smooth per-position enrichment ratios."""
        # gene_ratio = gene_ratio.replace(np.inf, np.nan)
        gene_ratio_mean = gene_ratio.mean(axis=1)
        gene_ratio_mean.index = gene_rpf.from_tis

        if self.size == 0:
            gene_ratio_smooth = pd.DataFrame(gene_ratio_mean, index=gene_rpf.from_tis, columns=['enrich'])
            return gene_ratio_mean, gene_ratio_smooth
        else:
            gene_ratio_smooth = savgol_filter(gene_ratio_mean, self.size, self.polyorder, mode='nearest')
            gene_ratio_smooth = pd.DataFrame(gene_ratio_smooth, index=gene_rpf.from_tis, columns=['enrich'])
            return gene_ratio_mean, gene_ratio_smooth

    def calculate_consensus_enrichment(self, raw_gene_rpm, ck_name, ip_name):
        """Calculate matched-replicate local enrichment for consensus calling.

        Parameters
        ----------
        raw_gene_rpm : pandas.DataFrame
            Per-codon RPM values indexed by position relative to the start codon.
        ck_name : list[str]
            Ordered control sample names.
        ip_name : list[str]
            Ordered IP sample names. Each IP sample is paired with the control
            sample at the same list position.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Replicate-specific enrichment profiles and a consensus profile with
            ``enrich`` and pointwise ``support`` columns.

        Notes
        -----
        This wrapper delegates to the shared SeRP enrichment implementation so
        peak calling and metaplot analysis use the same matched-replicate local
        enrichment definition.
        """
        return calculate_matched_local_enrichment(
            raw_gene_rpm=raw_gene_rpm,
            control_samples=ck_name,
            ip_samples=ip_name,
            window=self.consensus_window,
            pseudocount=self.pseudocount,
            threshold=self.enrich,
        )

    def consensus_analysis_positions(self, cds_region):
        """Return CDS positions eligible for consensus QC and peak calling."""
        positions = sorted(int(position) for position in cds_region)
        positions = [position for position in positions if position >= self.background]

        stop_trim = max(int(self.stop_trim), 0)
        if stop_trim > 0 and len(positions) > stop_trim:
            positions = positions[:-stop_trim]
        elif stop_trim > 0:
            positions = []

        return positions

    def _segment_threshold_regions(self, profile, allowed_positions):
        """Segment one enrichment profile into deterministic threshold regions.

        Core positions must satisfy ``enrich``. Short gaps may be bridged only
        when every gap position is present in the analysis region and remains at
        or above the lower ``collision``/edge threshold.
        """
        allowed = set(int(position) for position in allowed_positions)
        core_positions = [
            int(position)
            for position in profile.index
            if int(position) in allowed and float(profile.loc[position]) >= self.enrich
        ]
        core_positions.sort()
        if not core_positions:
            return []

        regions = []
        current = [core_positions[0]]

        def finalize(candidate):
            if not candidate:
                return
            start = int(candidate[0])
            end = int(candidate[-1])
            span = end - start + 1
            core_count = len(candidate)
            gap_count = span - core_count
            if core_count < self.width:
                return
            if span <= 0:
                return
            if gap_count / float(span) > self.proportion:
                return
            regions.append(
                {
                    "start": start,
                    "end": end,
                    "span": span,
                    "core_count": core_count,
                    "gap_count": gap_count,
                    "core_positions": list(candidate),
                }
            )

        for position in core_positions[1:]:
            previous = current[-1]
            gap_length = position - previous - 1
            bridge_gap = gap_length == 0

            if 0 < gap_length <= self.gaps:
                gap_positions = list(range(previous + 1, position))
                gap_available = all(gap in allowed and gap in profile.index for gap in gap_positions)
                if gap_available:
                    gap_values = profile.loc[gap_positions].astype(float)
                    bridge_gap = bool((gap_values >= self.collision).all())

            if bridge_gap:
                current.append(position)
            else:
                finalize(current)
                current = [position]

        finalize(current)
        return regions

    @staticmethod
    def _interval_overlap(start1, end1, start2, end2):
        """Return inclusive overlap length between two codon intervals."""
        return max(0, min(end1, end2) - max(start1, start2) + 1)

    def call_consensus_peak_regions(self, replicate_ratio, consensus_profile, cds_region):
        """Call replicate-supported SeRP peaks without positional permutation tests.

        Returns
        -------
        tuple[collections.OrderedDict, dict[int, int]]
            Historical peak structure plus accurate internal-gap counts for each
            consensus peak.
        """
        allowed_positions = self.consensus_analysis_positions(cds_region)
        if len(allowed_positions) < self.width:
            return OrderedDict(), {}

        replicate_regions = {
            sample: self._segment_threshold_regions(
                replicate_ratio[sample],
                allowed_positions,
            )
            for sample in replicate_ratio.columns
        }
        consensus_regions = self._segment_threshold_regions(
            consensus_profile["enrich"],
            allowed_positions,
        )

        replicate_number = max(len(replicate_regions), 1)
        required_replicates = max(
            1,
            int(np.ceil(replicate_number * self.min_support - 1e-12)),
        )
        min_overlap = int(self.min_overlap)
        if min_overlap <= 0:
            min_overlap = max(1, int(np.ceil(self.width / 2.0)))

        peak = OrderedDict()
        gap_counts = {}
        peak_num = 0

        for candidate in consensus_regions:
            supported = 0
            for regions in replicate_regions.values():
                if any(
                    self._interval_overlap(
                        candidate["start"],
                        candidate["end"],
                        region["start"],
                        region["end"],
                    )
                    >= min_overlap
                    for region in regions
                ):
                    supported += 1

            if supported < required_replicates:
                continue

            peak[peak_num] = {candidate["span"]: candidate["core_positions"]}
            gap_counts[peak_num] = candidate["gap_count"]
            peak_num += 1

        return peak, gap_counts

    def consensus_peak_statistics(self, peak, raw_gene_rpm, ck_name, ip_name, mrna):
        """Summarize peak RPM and test paired biological replicates.

        The statistical unit is one matched biological replicate pair. When at
        least three pairs are available, a one-sided paired t-test is calculated
        on log2 integrated peak RPM. With fewer than three pairs, the p-value is
        deliberately left unavailable because a two-replicate t-test is too
        unstable for meaningful inferential annotation. Statistical significance
        is never used to define consensus peaks.
        """
        pvalues = []
        ck_rpm = []
        ip_rpm = []
        rpm_dict = OrderedDict()

        for peak_num, peak_dict in peak.items():
            peak_positions = list(peak_dict.values())[0]
            peak_start = int(peak_positions[0])
            peak_end = int(peak_positions[-1])
            peak_range = [
                position
                for position in range(peak_start, peak_end + 1)
                if position in raw_gene_rpm.index
            ]

            ck_values = raw_gene_rpm.loc[peak_range, ck_name].sum(axis=0).astype(float)
            ip_values = raw_gene_rpm.loc[peak_range, ip_name].sum(axis=0).astype(float)
            ck_rpm.append(float(ck_values.mean()))
            ip_rpm.append(float(ip_values.mean()))

            names = mrna + "_peak_" + str(peak_num)
            rpm_dict[names] = list(
                map(
                    str,
                    raw_gene_rpm.loc[peak_range, ck_name + ip_name]
                    .sum(axis=0)
                    .to_list(),
                )
            )

            if len(ck_values) < 3 or len(ck_values) != len(ip_values):
                pvalues.append("-")
                continue

            integrated_pseudocount = self.pseudocount * max(len(peak_range), 1)
            ck_log = np.log2(ck_values.to_numpy() + integrated_pseudocount)
            ip_log = np.log2(ip_values.to_numpy() + integrated_pseudocount)
            paired_difference = ip_log - ck_log
            if np.std(paired_difference, ddof=1) == 0:
                pvalues.append("-")
                continue

            test = ttest_rel(ip_log, ck_log, alternative="greater")
            if np.isfinite(test.pvalue):
                pvalues.append(str(float(test.pvalue)))
            else:
                pvalues.append("-")

        return pvalues, ck_rpm, ip_rpm, rpm_dict

    @staticmethod
    def annotate_peak_region(utr5_region, cds_region, utr3_region, peak_start, peak_end):
        """Annotate a peak as 5UTR, CDS, 3UTR, or a spanning region."""
        # check the peak start site
        peak_start_loci = ''
        if peak_start in utr5_region:
            peak_start_loci = 'utr5'
        elif peak_start in cds_region:
            peak_start_loci = 'cds'
        elif peak_start in utr3_region:
            peak_start_loci = 'utr3'
        # check the peak end site
        peak_end_loci = ''
        if peak_end in utr5_region:
            peak_end_loci = 'utr5'
        elif peak_end in cds_region:
            peak_end_loci = 'cds'
        elif peak_end in utr3_region:
            peak_end_loci = 'utr3'

        if peak_start_loci == peak_end_loci:
            return peak_start_loci
        else:
            return peak_start_loci + ',' + peak_end_loci

    def extend_peak_edges(
        self,
        enrichment_profile,
        peak,
        max_extension=10,
        allowed_positions=None,
    ):
        """Extend consensus peaks through adjacent lower-threshold edge positions.

        Parameters
        ----------
        enrichment_profile : pandas.DataFrame
            Consensus enrichment profile containing an ``enrich`` column.
        peak : collections.OrderedDict
            Called consensus peaks.
        max_extension : int, default=10
            Maximum number of codons added independently to either peak edge.
        allowed_positions : list[int] or None, optional
            Positions eligible for extension. When provided, edges cannot extend
            into excluded start/stop regions.

        Returns
        -------
        tuple[list[int], list[int]]
            Inclusive left and right edge coordinates for each peak.

        Notes
        -----
        This method is used only by the consensus caller. It deliberately treats
        the lower threshold as a peak shoulder/edge criterion rather than evidence
        of ribosome collision.
        """
        available_positions = set(int(position) for position in enrichment_profile.index)
        if allowed_positions is not None:
            available_positions.intersection_update(
                int(position) for position in allowed_positions
            )
        max_extension = max(int(max_extension), 0)
        left_edges = []
        right_edges = []

        for peak_info in peak.values():
            peak_positions = list(peak_info.values())[0]
            peak_start = int(peak_positions[0])
            peak_end = int(peak_positions[-1])
            left_edge = peak_start
            right_edge = peak_end

            for step in range(1, max_extension + 1):
                position = peak_start - step
                if position not in available_positions:
                    break
                if float(enrichment_profile.loc[position, "enrich"]) < self.collision:
                    break
                left_edge = position

            for step in range(1, max_extension + 1):
                position = peak_end + step
                if position not in available_positions:
                    break
                if float(enrichment_profile.loc[position, "enrich"]) < self.collision:
                    break
                right_edge = position

            left_edges.append(left_edge)
            right_edges.append(right_edge)

        return left_edges, right_edges

    def extend_collision_regions(self, gene_ratio_smooth, peak):
        """Extend called peaks into neighboring lower-threshold edge positions.

        Notes
        -----
        ``collision`` is retained as the historical parameter name for output
        compatibility. The extended region is an enrichment shoulder/edge and
        should not be interpreted as direct evidence of ribosome collision.
        """
        gene_index = gene_ratio_smooth.index.tolist()
        peak_left_start, peak_left_end, peak_right_start, peak_right_end = [], [], [], []

        for peak_nums, peak_info in peak.items():
            peak_left_num = 0
            peak_right_num = 0

            peak_regions = list(peak_info.values())[0]
            # check the start position and gene start site
            start_posi = peak_regions[0] - gene_index[0]
            if start_posi == 0:
                pass
            else:
                left_shift_aa = 0
                for aa_num in range(start_posi):
                    left_shift_aa += 1
                    shift_aa_enrich = gene_ratio_smooth.loc[peak_regions[0] - aa_num - 1].enrich
                    if left_shift_aa >= 10:
                        break
                    elif shift_aa_enrich < self.collision:
                        break
                    elif shift_aa_enrich >= self.collision:
                        peak_left_num += 1

            # check the stop position and gene stop site
            stop_posi = gene_index[-1] - peak_regions[-1]
            if stop_posi == 0:
                pass
            else:
                right_shift_aa = 0
                for aa_num in range(stop_posi):
                    right_shift_aa += 1
                    shift_aa_enrich = gene_ratio_smooth.loc[peak_regions[-1] + aa_num + 1].enrich
                    if right_shift_aa >= 10:
                        break
                    elif shift_aa_enrich < self.collision:
                        break
                    elif shift_aa_enrich >= self.collision:
                        peak_right_num += 1

            peak_left_start.append(peak_regions[0] - peak_left_num)
            # peak_left_end.append(peak_regions[0])
            # peak_right_start.append(peak_regions[-1])
            peak_right_end.append(peak_regions[-1] + peak_right_num)

        return peak_left_start, peak_right_end

    @staticmethod
    def summarize_peak_enrichment(peak, gene_ratio_smooth):
        """Summarize maximum site, maximum fold, and mean fold for each peak."""
        max_peak_enrich = []
        max_peak_site = []
        mean_peak_enrich = []

        for peak_num, peak_dict in peak.items():
            peak_list = list(peak_dict.values())[0]
            gene_ratio_list = gene_ratio_smooth.loc[peak_list]

            # max_enrich = round(gene_ratio_list.max()[0], 4)
            max_enrich = round(gene_ratio_list.max().iloc[0], 4)
            max_peak_enrich.append(max_enrich)

            # max_site = str(gene_ratio_list.idxmax()[0])
            max_site = str(gene_ratio_list.idxmax().iloc[0])

            max_peak_site.append(max_site)

            mean_enrich = str(round(gene_ratio_list.mean().enrich, 4))
            mean_peak_enrich.append(mean_enrich)

        return max_peak_enrich, max_peak_site, mean_peak_enrich

    @staticmethod
    def rank_sum_positions(posi_list, raw_gene_rpm, ck_name, ip_name):
        """Run the legacy rank-sum test for one candidate position list."""
        rank_sum_results = []
        peak_range = posi_list
        ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
        ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
        t, p = ranksums(ck_gene_rpf, ip_gene_rpf)
        rank_sum_results.append(str(p))

        return rank_sum_results

    @staticmethod
    def rank_sum_peaks(peak, raw_gene_rpm, ck_name, ip_name, mrna):
        """Run the legacy rank-sum test for all called peaks."""
        rank_sum_results = []
        ck_rpm = []
        ip_rpm = []
        rpm_dict = OrderedDict()
        for peak_num, peak_dict in peak.items():
            peak_range = list(peak_dict.values())[0]
            ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
            ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
            t, p = ranksums(ck_gene_rpf, ip_gene_rpf)
            rank_sum_results.append(str(p))
            ck_rpm.append(ck_gene_rpf.sum())
            ip_rpm.append(ip_gene_rpf.sum())

            names = mrna + '_peak_' + str(peak_num)
            rpm_dict[names] = list(map(str, raw_gene_rpm.loc[peak_range][ck_name + ip_name].sum(axis=0).to_list()))

        return rank_sum_results, ck_rpm, ip_rpm, rpm_dict

    @staticmethod
    def ttest_positions(posi_list, raw_gene_rpm, ck_name, ip_name):
        """Run the legacy paired t-test for one candidate position list."""
        ttest_results = []
        peak_range = posi_list
        ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
        ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
        t, p = ttest_rel(ck_gene_rpf, ip_gene_rpf)
        ttest_results.append(str(p))

        return ttest_results

    @staticmethod
    def ttest_peaks(peak, raw_gene_rpm, ck_name, ip_name, mrna):
        """Run the legacy paired t-test and RPM summary for called peaks."""
        ttest_results = []
        ck_rpm = []
        ip_rpm = []
        rpm_dict = OrderedDict()
        for peak_num, peak_dict in peak.items():
            peak_range = list(peak_dict.values())[0]
            ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
            ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
            t, p = ttest_rel(ck_gene_rpf, ip_gene_rpf)
            ttest_results.append(str(p))
            ck_rpm.append(ck_gene_rpf.sum())
            ip_rpm.append(ip_gene_rpf.sum())

            names = mrna + '_peak_' + str(peak_num)
            rpm_dict[names] = list(map(str, raw_gene_rpm.loc[peak_range][ck_name + ip_name].sum(axis=0).to_list()))

        return ttest_results, ck_rpm, ip_rpm, rpm_dict

    def call_peak_regions(self, eligible_index, gene_ratio_smooth, raw_gene_rpm, ck_name, ip_name):
        """Call enriched regions with the historical SeRP peak segmentation logic."""
        peak = OrderedDict()
        pre_peak = OrderedDict()
        temp_store = []
        pre_peak_num = 1
        now_length = 0

        # if none gaps allowed
        if self.gaps == 0:
            for posi in range(len(eligible_index)):
                if not temp_store:
                    temp_store.append(eligible_index[posi])
                    now_length = 1
                else:
                    now_gap_length = eligible_index[posi] - eligible_index[posi - 1] - 1
                    # Check for gaps inside the potential peak.
                    if now_gap_length == 0:
                        temp_store.append(eligible_index[posi])
                        now_length += 1
                    else:
                        if now_length >= self.width:
                            pre_peak[pre_peak_num] = temp_store
                            pre_peak_num += 1
                            temp_store = [eligible_index[posi]]
                            now_length = 1

                        elif now_length < self.width:
                            temp_store = [eligible_index[posi]]
                            now_length = 1
            # save the last peak
            if now_length >= self.width:
                pre_peak[pre_peak_num] = temp_store

            # Traverse each potential peak.
            now_peak_num = 0
            for peak_num, peak_posi in pre_peak.items():
                # [peak length, peak position]
                peak[now_peak_num] = [peak_posi[-1] - peak_posi[0] + 1, peak_posi]
                now_peak_num += 1

        else:
            # This step used to search the potential peak but do not combine them
            # Traverse all potential peak regions and use the gaps_threshold for segmentation.
            for posi in range(len(eligible_index)):
                if not temp_store:
                    temp_store.append(eligible_index[posi])
                    now_length = 1
                else:
                    now_gap_length = eligible_index[posi] - eligible_index[posi - 1] - 1

                    # Check for gaps inside the potential peak.
                    if now_gap_length == 0:
                        temp_store.append(eligible_index[posi])
                        now_length += 1

                    # The loop skips when the enrichment in the gap lower than Collision.
                    elif now_gap_length <= self.gaps:
                        # if the gap enrich high than specified collision ratio
                        gap_fold_ratio = gene_ratio_smooth.loc[eligible_index[posi - 1] + 1:eligible_index[posi] - 1] - self.collision
                        # The loop continues when the enrichment in the gap higher than Collision.
                        if gap_fold_ratio.values.min() >= 0:
                            temp_store.append(eligible_index[posi])
                            now_length = now_length + now_gap_length + 1

                        # Save the results and break the current loop,
                        # when the enrichment in the gap lower than Collision and the
                        # current peak width is wider than the peak_width.
                        elif gap_fold_ratio.values.min() < 0 and now_length >= self.width:
                            pre_peak[pre_peak_num] = temp_store
                            pre_peak_num += 1
                            temp_store = [eligible_index[posi]]
                            now_length = 1

                        else:
                            temp_store = [eligible_index[posi]]
                            now_length = 1

                    # If the length of the consecutive gap is greater than the gap_threshold
                    elif now_gap_length > self.gaps:
                        if now_length >= self.width:
                            pre_peak[pre_peak_num] = temp_store
                            pre_peak_num += 1
                            temp_store = [eligible_index[posi]]
                            now_length = 1

                        elif now_length < self.width:
                            temp_store = [eligible_index[posi]]
                            now_length = 1

            if now_length >= self.width:
                pre_peak[pre_peak_num] = temp_store

            # Traverse each potential peak.
            now_peak_num = 0
            for peak_num, peak_posi in pre_peak.items():
                # if the current peak meets the conditions.
                if len(peak_posi) / (peak_posi[-1] - peak_posi[0] + 1) >= 1 - self.proportion:
                    peak[now_peak_num] = [peak_posi[-1] - peak_posi[0] + 1, peak_posi]
                    now_peak_num += 1
                # if the current peak does not meets the conditions.
                else:

                    temp_split_list = []  # split the pre peak region by the gap
                    temp_list = []
                    now_region_num = 0

                    # Split the remaining candidate peaks again.
                    for posi in range(len(peak_posi)):
                        if not temp_split_list:
                            temp_split_list.append(peak_posi[posi])
                        # If the gaps does not exist
                        elif peak_posi[posi] - peak_posi[posi - 1] == 1:
                            temp_split_list.append(peak_posi[posi])
                        # If the gaps exist
                        elif peak_posi[posi] - peak_posi[posi - 1] > 1:
                            temp_list.append(temp_split_list)
                            temp_split_list = [peak_posi[posi]]
                            now_region_num += 1

                    if temp_split_list:
                        temp_list.append(temp_split_list)

                    # Regroup the divided candidate peaks in order.
                    # The loop continues when there is still candidate peaks in pre_peak Dict.
                    region_num = len(temp_list)
                    permutations = OrderedDict()
                    for start_posi in range(region_num - 1):
                        # Exclude peaks whose width is less than 5 AA.
                        if len(temp_list[start_posi]) >= self.width:
                            if len(temp_list[start_posi]) / (
                                    temp_list[start_posi][0] - temp_list[start_posi][-1] + 1) >= 1 - self.proportion:
                                ttest_results = self.ttest_positions(temp_list[start_posi], raw_gene_rpm, ck_name,
                                                                         ip_name)
                                permutations[str(start_posi)] = [ttest_results, len(temp_list[start_posi]), start_posi,
                                                                 temp_list[start_posi]]
                            else:
                                pass
                        for stop_posi in range(start_posi + 1, region_num):
                            merge_region = []
                            for i in temp_list[start_posi:stop_posi + 1]:
                                merge_region.extend(i)

                            merge_len = int(merge_region[-1]) - int(merge_region[0]) + 1
                            # If the current permutation meets the conditions.
                            if len(merge_region) / merge_len >= 1 - self.proportion:
                                permutations_num = '_'.join([str(num) for num in range(start_posi, stop_posi + 1)])
                                ttest_results = self.ttest_positions(temp_list[start_posi], raw_gene_rpm, ck_name,
                                                                         ip_name)
                                permutations[permutations_num] = [ttest_results, merge_len, permutations_num,
                                                                  merge_region]
                            else:
                                pass

                    # Filter out the best results from different permutations.
                    while len(permutations.keys()) >= 1:
                        # If there is only one element in the list.
                        if len(permutations.keys()) == 1:
                            peak_info = list(permutations.values())[0]
                            peak[now_peak_num] = [peak_info[1], peak_info[-1]]
                            break
                        else:
                            # If there are multiple elements in the list.
                            result = [i for i in permutations.keys() if '_' in i]

                            # filter the optimal permutation
                            if result and not(self.keep_all):
                                now_best_permutations = []
                                for permutations_num, permutations_info in permutations.items():
                                    if len(now_best_permutations) == 0:
                                        now_best_permutations = permutations_info

                                    # if the length of current permutation large than the optimal one.
                                    elif permutations_info[1] > now_best_permutations[1]:
                                        # compare the p_value and gap proportion of the current permutation and the optimal one.
                                        gap_per = (int(permutations_info[3][-1]) - int(permutations_info[3][0])) / len(
                                            permutations_info[3])
                                        now_best_gap_per = (int(permutations_info[3][-1]) - int(
                                            permutations_info[3][0])) / len(permutations_info[3])
                                        if float(permutations_info[0][0]) > float(now_best_permutations[0][0]) and gap_per > now_best_gap_per:
                                            continue
                                        else:
                                            now_best_permutations = permutations_info

                                    # if the length of current permutation equal to the optimal one.
                                    elif permutations_info[1] == now_best_permutations[1]:
                                        if float(permutations_info[0][0]) < float(now_best_permutations[0][0]):
                                            now_best_permutations = permutations_info

                                        # Compare the p_value of the current permutation and the optimal one.
                                        elif float(permutations_info[0][0]) == float(now_best_permutations[0][0]):
                                            # compare the gap proportion of the current permutation and the optimal one.
                                            gap_per = (int(permutations_info[3][-1]) - int(permutations_info[3][0]))/len(permutations_info[3])
                                            now_best_gap_per = (int(permutations_info[3][-1]) - int(permutations_info[3][0]))/len(permutations_info[3])
                                            if gap_per < now_best_gap_per:
                                                now_best_permutations = permutations_info
                                            elif gap_per > now_best_gap_per:
                                                continue

                                            # Compare the start position of the current permutation and the optimal one.
                                            elif int(permutations_num.split('_')[0]) < int(now_best_permutations[2].split('_')[0]):
                                                now_best_permutations = permutations_info

                                        elif float(permutations_info[0][0]) > float(now_best_permutations[0][0]):
                                            continue
                                        else:
                                            # now_best_permutations
                                            sys.stdout.write('unknown')

                                    # If the width of the current permutation is narrower than the optimal one.
                                    elif permutations_info[1] < now_best_permutations[1]:
                                        continue

                                # save the best permutations
                                peak_info = now_best_permutations[-1]
                                peak[now_peak_num] = [int(peak_info[-1]) - int(peak_info[0]) + 1, peak_info]
                                now_peak_num += 1
                                permutations.pop(now_best_permutations[2])

                                # delete the permutations and original elements
                                permutations_elements = now_best_permutations[2].split('_')
                                shift_posi = 0
                                for start_posi in permutations_elements:
                                    shift_posi += 1
                                    if permutations.get(start_posi):
                                        permutations.pop(start_posi)
                                    for stop_posi in permutations_elements[shift_posi:]:
                                        permutations_num = '_'.join([str(num) for num in range(int(start_posi), int(stop_posi) + 1)])
                                        if permutations.get(permutations_num):
                                            permutations.pop(permutations_num)
                            # If the elements in the list are all single
                            else:
                                for permutations_num, permutations_info in permutations.items():
                                    peak[now_peak_num] = [len(permutations_info[-1]), permutations_info[-1]]
                                    now_peak_num += 1
                                break

        # sort the peak by region
        peak_sort = OrderedDict()
        if not bool(peak):
            return peak_sort
        elif len(peak.keys()) == 1:
            peak_sort[0] = {peak[0][0]: peak[0][1]}
            return peak_sort
        else:
            peak_sort1 = sorted(peak.items(), key=lambda d: d[1][1][0])
            new_peak_num = 0
            for peak_info in peak_sort1:
                peak_sort[new_peak_num] = {peak_info[1][0]: peak_info[1][1]}
                new_peak_num += 1

            return peak_sort