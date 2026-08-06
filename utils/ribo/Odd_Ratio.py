#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-06
# Version: dev003
# Function: Detect local ribosome pauses and differential codon enrichment.
# Input: Frame-resolved RPF density data with control and treatment replicates.
# Output: Site-level pause evidence, differential statistics, summaries, and figures.

"""Local and differential ribosome-pause analysis."""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm

from . import RPFs


def _benjamini_hochberg(
    pvalue: np.ndarray | list[float],
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Benjamini-Hochberg rejection flags and adjusted p-values."""
    pvalue_array = np.asarray(pvalue, dtype=float)
    if pvalue_array.size == 0:
        return np.array([], dtype=bool), np.array([], dtype=float)

    order = np.argsort(pvalue_array, kind='mergesort')
    ranked = pvalue_array[order]
    rank = np.arange(1, ranked.size + 1, dtype=float)
    adjusted_ranked = ranked * ranked.size / rank
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0.0, 1.0)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted <= alpha, adjusted


class OddRatio:
    """Run local-pause and differential codon-enrichment analysis."""

    def __init__(self, args):
        # input and output files
        self.rpf_file = args.rpf
        self.output_prefix = args.output
        self.detail = args.detail

        # opts for odd ratio
        self.gene = args.list
        self.site = args.site
        self.frame = args.frame
        self.tis = args.tis
        self.tts = args.tts
        
        self.normal = args.normal
        self.flag = args.fdr
        self.pvalue = args.value
        self.fillna = args.zero
        self.test = args.test
        self.pseudocount = args.pseudocount
        self.min_log2_or = args.min_log2_or

        self.local_window = args.local_window
        self.pause_score = args.pause_score
        self.min_site_rpf = args.min_site_rpf
        self.min_local_coverage = args.min_local_coverage
        self.min_pause_replicates = args.min_pause_replicates

        self.thread = args.thread

        # samples
        self.control = [item.strip() for item in args.control.split(',') if item.strip()]
        self.treat = [item.strip() for item in args.treat.split(',') if item.strip()]
        self.rpf_num = args.min

        duplicated_samples = sorted(set(self.control).intersection(self.treat))
        if duplicated_samples:
            raise ValueError(
                "Samples cannot occur in both control and treatment groups: {samples}".format(
                    samples=", ".join(duplicated_samples)
                )
            )

        self.sample_name = self.control + self.treat
        self.sample_num = len(self.sample_name)

        # for the high expression gene filter
        self.merged_rpf = None
        self.count_rpf = None
        self.raw_rpf = None
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None
        self.dim2_df = None
        self.merge_odd_ratio = None

        # summary stalling codon
        self.scale = args.scale
        self.codon_odd_ratio = None
        self.codon_odd_ratio_num = None
        self.codon_odd_ratio_prop = None

        # codon list
        self.stop = args.stop
        self.codon_dict, self.codon_table = RPFs.codon_table()

        if self.stop:
            del self.codon_dict['TAA']
            del self.codon_dict['TAG']
            del self.codon_dict['TGA']
            
            self.codon_table = self.codon_table[self.codon_table['Abbr'] != 'Stop']

    def read_rpf(self):
        """Import RPF density while preserving the original count matrix.

        Returns
        -------
        None
            The function updates internal state or writes output files in place.

        Notes
        -----
        Legacy input information:
        - ``self.rpf_file``: RPF density file.
        - ``self.site``: the site of the codon.
        - ``self.frame``: the frame of the codon.
        - ``self.tis``: the start site of the CDS.
        - ``self.tts``: the end site of the CDS.
        - ``self.sample_num``: the number of the samples.
        - ``self.sample_name``: the name of the samples.
        - ``self.gene``: the gene list.
        - ``self.rpf_num``: the minimal number of the RPFs.

        Legacy return information:
        - merged_rpf, sample_name, sample_num, total_rpf_num, gene_rpf_sum, high_gene,
        high_rpf.

        Zero counts are retained because a zero in one condition can be the
        defining signal of a condition-specific pause. Positions with zero
        counts in every selected sample are removed only after transcript totals
        have been calculated. Statistical inference always uses raw counts.
        """

        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file,
                                      sites=self.site, frame=self.frame,
                                      tis=self.tis, tts=self.tts,
                                      sample_num=self.sample_num, sample_name=self.sample_name,
                                      gene=self.gene, rpf_num=self.rpf_num)
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3].copy()
        self.total_rpf_num = rpf_results[4]
        self.high_gene = rpf_results[6]
        del rpf_results

        self.merged_rpf = self.merged_rpf.loc[
            self.merged_rpf['codon'].isin(self.codon_dict.keys())
        ].copy()
        self.high_rpf = self.merged_rpf.loc[
            self.merged_rpf['name'].isin(self.high_gene)
        ].copy()
        self.high_rpf.sort_values(
            by=['name', 'from_tis', 'now_nt'],
            kind='mergesort',
            inplace=True,
        )
        self.high_rpf.reset_index(drop=True, inplace=True)

        self.count_rpf = self.high_rpf.copy()
        self.gene_rpf_sum = self.count_rpf.groupby(
            'name', sort=False
        )[self.sample_name].sum()

        covered = self.count_rpf[self.sample_name].sum(axis=1) > 0
        self.count_rpf = self.count_rpf.loc[covered].reset_index(drop=True)
        self.high_rpf = self.high_rpf.loc[covered].reset_index(drop=True)

        if self.normal:
            for sample in self.sample_name:
                total = float(self.total_rpf_num.get(sample, 0))
                rpm_column = sample + '_rpm'
                if total > 0:
                    self.high_rpf.loc[:, rpm_column] = (
                        self.high_rpf[sample].astype(float) * 1e6 / total
                    )
                else:
                    self.high_rpf.loc[:, rpm_column] = 0.0

        if self.fillna:
            print(
                "Warning: -z/--zero is deprecated; raw zero counts are retained "
                "without replacement.",
                flush=True,
            )

        print(
            "Retained transcripts={transcripts:,}, tested positions={positions:,}; "
            "zero counts preserved.".format(
                transcripts=self.count_rpf['name'].nunique(),
                positions=len(self.count_rpf),
            ),
            flush=True,
        )

    def make_two_dimensional_table(self):
        """Build raw site-versus-rest-of-transcript count tables.

        Returns
        -------
        None
            The function updates internal state or writes output files in place.

        Notes
        -----
        Legacy input information:
        - ``self.high_rpf``: high expression rpf table.
        - ``self.control``: control sample table.
        - ``self.treat``: treat sample table.

        Legacy return information:
        - ``self.gene_rpf_sum``: sum of high expression gene table.

        Workflow
        --------
        - 1 --> sum the rpf count of each gene with control and treatment group.
        - 2 --> calculate the gene delta.
        - delta = each gene sum - each codon of this gene.
        """

        self.high_rpf = self.high_rpf.reset_index(drop=True)
        self.count_rpf = self.count_rpf.reset_index(drop=True)

        self.high_rpf.loc[:, 'control'] = self.count_rpf[self.control].sum(axis=1)
        self.high_rpf.loc[:, 'treat'] = self.count_rpf[self.treat].sum(axis=1)

        gene_sums = self.gene_rpf_sum.copy()
        gene_sums.loc[:, 'control_sum'] = gene_sums[self.control].sum(axis=1)
        gene_sums.loc[:, 'treat_sum'] = gene_sums[self.treat].sum(axis=1)
        gene_sums = gene_sums.reset_index()

        self.high_rpf = pd.merge(
            self.high_rpf,
            gene_sums[['name', 'control_sum', 'treat_sum']],
            on='name',
            how='left',
            validate='many_to_one',
        )

        for sample in self.sample_name:
            self.high_rpf.loc[:, sample + '_gene_sum'] = (
                self.high_rpf['name'].map(self.gene_rpf_sum[sample]).to_numpy()
            )

        invalid_control = self.high_rpf['control_sum'] <= 0
        invalid_treat = self.high_rpf['treat_sum'] <= 0
        if invalid_control.any() or invalid_treat.any():
            keep = ~(invalid_control | invalid_treat)
            removed = int((~keep).sum())
            print(
                "Removed positions from transcripts without RPFs in one group: "
                "{removed:,}.".format(removed=removed),
                flush=True,
            )
            self.high_rpf = self.high_rpf.loc[keep].reset_index(drop=True)
            self.count_rpf = self.count_rpf.loc[keep].reset_index(drop=True)
        if self.high_rpf.empty:
            raise ValueError(
                "No positions remain after requiring transcript-level RPFs in "
                "both control and treatment groups."
            )

    @staticmethod
    def _safe_pause_score(site_count: pd.Series, local_mean: pd.Series) -> pd.Series:
        """Return site count divided by local-window mean with safe zeros."""
        score = np.divide(
            site_count.to_numpy(dtype=float),
            local_mean.to_numpy(dtype=float),
            out=np.zeros(len(site_count), dtype=float),
            where=local_mean.to_numpy(dtype=float) > 0,
        )
        return pd.Series(score, index=site_count.index)

    def calc_local_pause(self):
        """Calculate PausePred-style local scores and replicate support.

        The local score is the site RPF count divided by the centered-window
        mean count. Window coverage is the fraction of codon positions with at
        least one RPF. Group-level metrics use pooled counts, whereas replicate
        support is calculated independently for every sample.
        """
        metric_columns = ['control', 'treat'] + self.sample_name
        metric_values = pd.concat(
            [
                self.high_rpf[['control', 'treat']],
                self.count_rpf[self.sample_name],
            ],
            axis=1,
        )
        transcript_group = self.high_rpf['name']
        rolling_mean = (
            metric_values.groupby(transcript_group, sort=False)
            .rolling(
                window=self.local_window,
                center=True,
                min_periods=1,
            )
            .mean()
            .reset_index(level=0, drop=True)
            .sort_index()
        )
        rolling_coverage = (
            metric_values.gt(0).astype(float)
            .groupby(transcript_group, sort=False)
            .rolling(
                window=self.local_window,
                center=True,
                min_periods=1,
            )
            .mean()
            .reset_index(level=0, drop=True)
            .sort_index()
        )
        rolling_mean = rolling_mean.loc[:, metric_columns]
        rolling_coverage = rolling_coverage.loc[:, metric_columns]

        for group in ('control', 'treat'):
            local_mean = rolling_mean[group]
            local_coverage = rolling_coverage[group]
            self.high_rpf.loc[:, group + '_local_mean'] = local_mean
            self.high_rpf.loc[:, group + '_local_coverage'] = local_coverage
            self.high_rpf.loc[:, group + '_pause_score'] = self._safe_pause_score(
                self.high_rpf[group], local_mean
            )

        for group, samples in (('control', self.control), ('treat', self.treat)):
            support = np.zeros(len(self.high_rpf), dtype=np.int16)
            for sample in samples:
                local_mean = rolling_mean[sample]
                local_coverage = rolling_coverage[sample]
                sample_score = self._safe_pause_score(
                    self.count_rpf[sample], local_mean
                )
                support += (
                    (self.count_rpf[sample].to_numpy() >= self.min_site_rpf)
                    & (sample_score.to_numpy() >= self.pause_score)
                    & (local_coverage.to_numpy() >= self.min_local_coverage)
                ).astype(np.int16)
            self.high_rpf.loc[:, group + '_pause_replicates'] = support

            required_replicates = min(self.min_pause_replicates, len(samples))
            self.high_rpf.loc[:, group + '_local_pause'] = (
                (self.high_rpf[group] >= self.min_site_rpf)
                & (self.high_rpf[group + '_pause_score'] >= self.pause_score)
                & (
                    self.high_rpf[group + '_local_coverage']
                    >= self.min_local_coverage
                )
                & (
                    self.high_rpf[group + '_pause_replicates']
                    >= required_replicates
                )
            )

    def calc_odd_ratio(self):
        """Calculate Haldane-Anscombe-corrected odds ratios.

        Returns
        -------
        None
            The function updates internal state or writes output files in place.

        Notes
        -----
        Legacy input information:
        - ``self.high_rpf``: high expression rpf table.

        Legacy return information:
        - ``self.high_rpf``: contain the high rpf and odd-ratio.

        Workflow
        --------
        1. Calculate the delta of each gene.
        2. Calculate the odd ratio.
        """
        
        control = self.high_rpf['control'].to_numpy(dtype=float)
        treat = self.high_rpf['treat'].to_numpy(dtype=float)
        control_rest = (
            self.high_rpf['control_sum'].to_numpy(dtype=float) - control
        )
        treat_rest = self.high_rpf['treat_sum'].to_numpy(dtype=float) - treat
        pseudocount = float(self.pseudocount)

        log_odds_ratio = (
            np.log(treat + pseudocount)
            + np.log(control_rest + pseudocount)
            - np.log(control + pseudocount)
            - np.log(treat_rest + pseudocount)
        )
        self.high_rpf.loc[:, 'odd'] = np.exp(
            np.clip(log_odds_ratio, -700, 700)
        )
        self.high_rpf.loc[:, 'log2_odd'] = log_odds_ratio / np.log(2.0)
        self.high_rpf.loc[:, 'control_rest'] = control_rest
        self.high_rpf.loc[:, 'treat_rest'] = treat_rest

    def calc_differential_test(self):
        """Run a vectorized pooled or replicate-aware Wald test.

        Returns
        -------
        None
            Differential statistics are added to ``self.high_rpf``.

        Workflow
        --------
        1. Estimate group-specific site proportions within each transcript.
        2. Estimate beta-binomial overdispersion from replicate residuals.
        3. Test the group logit difference with a Wald statistic.
        4. Correct p-values with Benjamini-Hochberg.
        """
        use_beta_binomial = self.test == 'beta-binomial' or (
            self.test == 'auto'
            and min(len(self.control), len(self.treat)) >= 2
        )
        test_method = 'beta_binomial_wald' if use_beta_binomial else 'pooled_wald'

        site_count = self.count_rpf[self.sample_name].to_numpy(dtype=float)
        gene_count = self.high_rpf[
            [sample + '_gene_sum' for sample in self.sample_name]
        ].to_numpy(dtype=float)
        group_index = np.array(
            [0] * len(self.control) + [1] * len(self.treat),
            dtype=np.int8,
        )

        group_site = np.column_stack(
            (
                site_count[:, group_index == 0].sum(axis=1),
                site_count[:, group_index == 1].sum(axis=1),
            )
        )
        group_total = np.column_stack(
            (
                gene_count[:, group_index == 0].sum(axis=1),
                gene_count[:, group_index == 1].sum(axis=1),
            )
        )
        group_probability = (
            group_site + self.pseudocount
        ) / (group_total + 2.0 * self.pseudocount)

        if use_beta_binomial:
            expected_probability = group_probability[:, group_index]
            expected_count = gene_count * expected_probability
            binomial_variance = (
                gene_count
                * expected_probability
                * (1.0 - expected_probability)
            )
            valid = (gene_count > 0) & (binomial_variance > 0)
            pearson_component = np.divide(
                np.square(site_count - expected_count),
                binomial_variance,
                out=np.zeros_like(site_count, dtype=float),
                where=valid,
            )
            pearson_statistic = pearson_component.sum(axis=1)
            valid_sample_number = valid.sum(axis=1)
            valid_group_number = (
                valid[:, group_index == 0].any(axis=1).astype(int)
                + valid[:, group_index == 1].any(axis=1).astype(int)
            )
            residual_df = np.maximum(
                valid_sample_number - valid_group_number,
                0,
            )
            rho_denominator = (
                np.maximum(gene_count - 1.0, 0.0) * valid
            ).sum(axis=1)
            overdispersion = np.divide(
                pearson_statistic - residual_df,
                rho_denominator,
                out=np.zeros(len(self.high_rpf), dtype=float),
                where=rho_denominator > 0,
            )
            overdispersion = np.clip(overdispersion, 0.0, 0.99)
        else:
            overdispersion = np.zeros(len(self.high_rpf), dtype=float)

        probability_variance = np.zeros((len(self.high_rpf), 2), dtype=float)
        for group_id in (0, 1):
            selected = group_index == group_id
            sample_total = gene_count[:, selected]
            probability = group_probability[:, group_id]
            count_variance = (
                sample_total
                * probability[:, None]
                * (1.0 - probability[:, None])
                * (
                    1.0
                    + np.maximum(sample_total - 1.0, 0.0)
                    * overdispersion[:, None]
                )
            )
            total = group_total[:, group_id]
            probability_variance[:, group_id] = np.divide(
                count_variance.sum(axis=1),
                np.square(total),
                out=np.zeros(len(self.high_rpf), dtype=float),
                where=total > 0,
            )

        logit_variance = np.divide(
            probability_variance,
            np.square(group_probability * (1.0 - group_probability)),
            out=np.zeros_like(probability_variance),
            where=(group_probability > 0) & (group_probability < 1),
        )
        standard_error = np.sqrt(logit_variance.sum(axis=1))
        logit_difference = (
            np.log(group_probability[:, 1] / (1.0 - group_probability[:, 1]))
            - np.log(group_probability[:, 0] / (1.0 - group_probability[:, 0]))
        )
        statistic = np.divide(
            logit_difference,
            standard_error,
            out=np.zeros(len(self.high_rpf), dtype=float),
            where=standard_error > 0,
        )
        pvalue = 2.0 * norm.sf(np.abs(statistic))
        pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
        bhfdr = _benjamini_hochberg(pvalue, alpha=self.pvalue)

        self.high_rpf.loc[:, 'test_method'] = test_method
        self.high_rpf.loc[:, 'overdispersion'] = overdispersion
        self.high_rpf.loc[:, 'statistic'] = statistic
        self.high_rpf.loc[:, 'pvalue'] = pvalue
        self.high_rpf.loc[:, 'bhfdr'] = bhfdr[1]
        self.high_rpf.loc[:, 'flag'] = bhfdr[0]

        print(
            "Differential test={method}, positions={positions:,}, "
            "median_overdispersion={rho:.6f}.".format(
                method=test_method,
                positions=len(self.high_rpf),
                rho=float(np.median(overdispersion)),
            ),
            flush=True,
        )

    def calc_chi2_test2(self):
        """Run the differential test through the historical method name."""
        return self.calc_differential_test()

    def calc_chi2_test(self):
        """Run the differential test through the historical method name."""
        return self.calc_differential_test()

    def classify_pause(self):
        """Assign biologically interpretable local/differential pause classes."""
        significant = (
            (self.high_rpf[self.flag] < self.pvalue)
            & (self.high_rpf['log2_odd'].abs() >= self.min_log2_or)
        )
        control_pause = self.high_rpf['control_local_pause'].astype(bool)
        treat_pause = self.high_rpf['treat_local_pause'].astype(bool)
        treatment_enriched = significant & (self.high_rpf['log2_odd'] > 0)
        control_enriched = significant & (self.high_rpf['log2_odd'] < 0)

        pause_class = np.full(
            len(self.high_rpf),
            'not_significant',
            dtype=object,
        )
        pause_class[(control_pause | treat_pause) & ~significant] = (
            'local_pause_without_significant_shift'
        )
        pause_class[control_pause & treat_pause & ~significant] = (
            'constitutive_pause'
        )
        pause_class[significant] = 'relative_shift_without_strong_local_peak'
        pause_class[treatment_enriched & treat_pause] = (
            'treatment_enriched_pause'
        )
        pause_class[control_enriched & control_pause] = (
            'control_enriched_pause'
        )

        self.high_rpf.loc[:, 'differential'] = significant
        self.high_rpf.loc[:, 'pause_class'] = pause_class

    def output_odd_ratio(self):
        """Write differential positions, local pauses, and optional details."""
        all_result_file = self.output_prefix + '_codon_odd_ratio_all.txt'
        differential_file = self.output_prefix + '_codon_odd_ratio.txt'
        local_pause_file = self.output_prefix + '_codon_local_pause.txt'

        self.high_rpf.loc[:, 'control_group'] = ','.join(self.control)
        self.high_rpf.loc[:, 'treat_group'] = ','.join(self.treat)
        self.high_rpf.loc[:, 'analysis_site'] = self.site
        self.high_rpf.loc[:, 'analysis_frame'] = self.frame
        self.high_rpf.loc[:, 'local_window_codons'] = self.local_window
        self.high_rpf.loc[:, 'pause_score_threshold'] = self.pause_score
        self.high_rpf.loc[:, 'min_site_rpf_threshold'] = self.min_site_rpf
        self.high_rpf.loc[:, 'min_local_coverage_threshold'] = self.min_local_coverage

        if self.detail:
            self.high_rpf.to_csv(all_result_file, sep='\t', index=False)
        self.merge_odd_ratio = self.high_rpf.loc[
            self.high_rpf['differential']
        ].copy()
        self.merge_odd_ratio.to_csv(differential_file, sep='\t', index=False)

        local_pause = self.high_rpf.loc[
            self.high_rpf['control_local_pause']
            | self.high_rpf['treat_local_pause']
        ].copy()
        local_pause.to_csv(local_pause_file, sep='\t', index=False)

        print(
            "Output positions={positions:,}, differential={differential:,}, "
            "local_pauses={local_pauses:,}.".format(
                positions=len(self.high_rpf),
                differential=len(self.merge_odd_ratio),
                local_pauses=len(local_pause),
            ),
            flush=True,
        )

    @staticmethod
    def scale_method(scale, oddratio):
        """Scale paired summary columns for visualization."""
        if scale == 'minmax':
            denominator = oddratio.max().replace(0, np.nan)
            relative_oddratio = oddratio.div(denominator, axis=1).fillna(0)
        elif scale == 'zscore':
            denominator = oddratio.std().replace(0, np.nan)
            relative_oddratio = (
                (oddratio - oddratio.mean()).div(denominator, axis=1).fillna(0)
            )
        else:
            raise ValueError("Unknown scale method: {scale}".format(scale=scale))
        return relative_oddratio
    
    def summarize_odd_ratio(self):
        """Summarize differential calls for every sense codon."""
        codon_odd_file = self.output_prefix + '_sum_codon_odd_ratio.txt'

        summary = self.codon_table.reset_index().rename(columns={'codon': 'Codon'})
        grouped = self.high_rpf.groupby('codon', sort=False)
        summary_metrics = pd.DataFrame(
            {
                'codon_sum': grouped.size(),
                'control_codon_sum': self.high_rpf.loc[
                    self.high_rpf['log2_odd'] < 0
                ].groupby('codon').size(),
                'treat_codon_sum': self.high_rpf.loc[
                    self.high_rpf['log2_odd'] > 0
                ].groupby('codon').size(),
                'control_mean_odd_ratio': self.high_rpf.loc[
                    self.high_rpf['log2_odd'] < 0
                ].groupby('codon')['odd'].mean(),
                'treat_mean_odd_ratio': self.high_rpf.loc[
                    self.high_rpf['log2_odd'] > 0
                ].groupby('codon')['odd'].mean(),
                'control_number': self.high_rpf.loc[
                    self.high_rpf['differential']
                    & (self.high_rpf['log2_odd'] < 0)
                ].groupby('codon').size(),
                'treat_number': self.high_rpf.loc[
                    self.high_rpf['differential']
                    & (self.high_rpf['log2_odd'] > 0)
                ].groupby('codon').size(),
            }
        ).reset_index().rename(columns={'codon': 'Codon'})
        codon_odd_ratio = summary.merge(
            summary_metrics,
            on='Codon',
            how='left',
            validate='one_to_one',
        )
        count_columns = [
            'codon_sum',
            'control_codon_sum',
            'treat_codon_sum',
            'control_number',
            'treat_number',
        ]
        codon_odd_ratio.loc[:, count_columns] = (
            codon_odd_ratio[count_columns].fillna(0).astype(int)
        )

        codon_odd_ratio_prop = codon_odd_ratio.loc[
            :, ['control_number', 'treat_number']
        ].div(codon_odd_ratio['codon_sum'].replace(0, np.nan), axis=0) * 100
        codon_odd_ratio_prop = codon_odd_ratio_prop.fillna(0)
        codon_odd_ratio_prop.columns = ['control_proportion', 'treat_proportion']

        codon_odd_ratio_num_rel = self.scale_method(
            self.scale,
            codon_odd_ratio.loc[:, ['control_number', 'treat_number']],
        )
        codon_odd_ratio_prop_rel = self.scale_method(self.scale, codon_odd_ratio_prop)
        codon_odd_ratio_num_rel.columns = ['control_number_relative', 'treat_number_relative']
        codon_odd_ratio_prop_rel.columns = ['control_proportion_relative', 'treat_proportion_relative']

        self.codon_odd_ratio = pd.concat(
            [
                codon_odd_ratio,
                codon_odd_ratio_prop,
                codon_odd_ratio_num_rel,
                codon_odd_ratio_prop_rel,
            ],
            axis=1,
        )
        self.codon_odd_ratio.sort_values(by=['Abbr', 'Codon'], inplace=True)

        self.codon_odd_ratio['control_group'] = ",".join(self.control)
        self.codon_odd_ratio['treat_group'] = ",".join(self.treat)

        self.codon_odd_ratio_num = pd.melt(self.codon_odd_ratio, 
                                    id_vars=['Codon', 'AA', 'Abbr'], value_vars=['control_number', 'treat_number'],
                                    var_name='groups', value_name='Number')
        self.codon_odd_ratio_num['Codon'] = self.codon_odd_ratio_num['Codon'] + '(' + self.codon_odd_ratio_num['Abbr'] + ')'

        self.codon_odd_ratio_prop = pd.melt(self.codon_odd_ratio, 
                                    id_vars=['Codon', 'AA', 'Abbr'], value_vars=['control_proportion', 'treat_proportion'],
                                    var_name='groups', value_name='Proportion')
        self.codon_odd_ratio_prop['Codon'] = self.codon_odd_ratio_prop['Codon'] + '(' + self.codon_odd_ratio_prop['Abbr'] + ')'

        self.codon_odd_ratio.to_csv(codon_odd_file, sep='\t', index=False)


    def draw_odd_ratio_bar(self):
        """Draw the odd ratio barplot.

        Returns
        -------
        None
            The function updates internal state or writes output files in place.

        Notes
        -----
        Legacy input information:
        - ``self.codon_odd_ratio_num``: the count of codon odd ratio.
        - ``self.codon_odd_ratio_prop``: the proportion of codon odd ratio.

        Legacy return information:
        - odd_ratio_barplot.pdf.
        """

        out_pdf = self.output_prefix + "_odd_barplot.pdf"
        out_png = self.output_prefix + "_odd_barplot.png"
        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15, 12), dpi=300)
        sns.barplot(data=self.codon_odd_ratio_num, x='Codon', y='Number', hue='groups', ax=ax[0])
        sns.barplot(data=self.codon_odd_ratio_prop, x='Codon', y='Proportion', hue='groups', ax=ax[1])

        # add the gene annotation here, need to import the txt file
        ax[0].set_ylabel('Number')
        ax[1].set_ylabel('Proportion (%)')
        ax[1].set_xlabel('Codon')

        plt.sca(ax[0])
        plt.xticks(rotation=90)
        plt.sca(ax[1])
        plt.xticks(rotation=90)

        plt.suptitle('Significant differential pausing codon')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_odd_ratio_line(self):
        """Draw the odd ratio lineplot.

        Returns
        -------
        None
            The function updates internal state or writes output files in place.

        Notes
        -----
        Legacy input information:
        - ``self.codon_odd_ratio_num``: the count of codon odd ratio.
        - ``self.codon_odd_ratio_prop``: the proportion of codon odd ratio.

        Legacy return information:
        - odd_ratio_lineplot.pdf.
        """
        out_pdf = self.output_prefix + "_odd_lineplot.pdf"
        out_png = self.output_prefix + "_odd_lineplot.png"

        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 6), dpi=300)

        flag = 1
        for aa in self.codon_odd_ratio_num['AA'].unique().tolist():
            tmp = self.codon_odd_ratio_num.loc[self.codon_odd_ratio_num['AA']==aa, :]
            if flag == 1:
                sns.lineplot(data=tmp, x='Codon', y='Number', linewidth=1.5, hue='groups', style='groups',
                              markers=['d','o'], legend = True, ax=ax[0])
                flag += 1
            else:
                sns.lineplot(data=tmp, x='Codon', y='Number', linewidth=1.5, hue='groups', style='groups',
                              markers=['d', 'o'], legend = False, ax=ax[0])

        flag = 1
        for aa in self.codon_odd_ratio_prop['AA'].unique().tolist():
            tmp = self.codon_odd_ratio_prop.loc[self.codon_odd_ratio_prop['AA']==aa, :]
            if flag == 1:
                sns.lineplot(data=tmp, x='Codon', y='Proportion', linewidth=1.5, hue='groups', style='groups',
                             markers=['d', 'o'], legend = True, ax=ax[1])
                flag += 1
            else:
                sns.lineplot(data=tmp, x='Codon', y='Proportion', linewidth=1.5, hue='groups', style='groups',
                             markers=['d', 'o'], legend = False, ax=ax[1])

        # add the gene annotation here, need to import the txt file
        ax[0].legend(loc=0, ncol=2, fontsize=6)
        ax[1].legend(loc=0, ncol=2, fontsize=6)

        ax[0].set_ylabel('Number')
        ax[1].set_ylabel('Proportion (%)')
        ax[1].set_xlabel('Codon')

        plt.sca(ax[0])
        plt.xticks(rotation=90)
        plt.sca(ax[1])
        plt.xticks(rotation=90)

        plt.suptitle('Number of sig. diff. pausing codon')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_odd_ratio_scatter(self):
        """Draw the odd ratio scatterplot.

        Returns
        -------
        None
            The function updates internal state or writes output files in place.

        Notes
        -----
        Legacy input information:
        - ``self.codon_odd_ratio_num``: the count of codon odd ratio.
        - ``self.codon_odd_ratio_prop``: the proportion of codon odd ratio.

        Legacy return information:
        - odd_ratio_scatterplot.pdf.
        """

        out_pdf = self.output_prefix + "_odd_scatter.pdf"
        out_png = self.output_prefix + "_odd_scatter.png"

        # count_df = self.codon_odd_ratio_num.pivot_table(index=['Codon', 'AA', 'Abbr'], columns='groups', values='Number').reset_index()
        # per_df = self.codon_odd_ratio_prop.pivot_table(index=['Codon', 'AA', 'Abbr'], columns='groups', values='Proportion').reset_index()

        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 6.5), dpi=300)
        
        count_plot = sns.scatterplot(data=self.codon_odd_ratio, x='control_number', y='treat_number', 
                                     s = 80, hue='AA', legend = 'brief', ax=ax[0])
        count_plot.legend(loc=8, bbox_to_anchor=(0.5, -0.5), ncol=7, fontsize=6)

        per_plot = sns.scatterplot(data=self.codon_odd_ratio, x='control_proportion', y='treat_proportion', 
                                   s = 80, hue='AA', legend = 'brief', ax=ax[1])
        per_plot.legend(loc=8, bbox_to_anchor=(0.5, -0.5), ncol=7, fontsize=6)

        # add the gene annotation here, need to import the txt file
        ax[0].set_title('Number')
        ax[1].set_title('Proportion (%)')

        plt.suptitle('Number of sig. diff. pausing codon')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)