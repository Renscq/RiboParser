#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : SeRP.py


import os
import sys
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import polars as pl

from Bio.Seq import Seq

from scipy.signal import savgol_filter
from scipy.stats import ranksums
from scipy.stats import ttest_rel


class SeRP(object):
    def __init__(self, args):
        self.tmp = None

        # input files
        self.rpf_file = args.rpf
        self.total_rpf = args.norm
        self.anno = args.anno
        self.gene_list = []
        self.gene_dict = OrderedDict()

        # parameters for analysis
        self.ck = args.control
        self.ip = args.ip
        self.ck_name = self.ck.split(',')
        self.ip_name = self.ip.split(',')

        # parameters for mrna import
        self.gene_list = None
        self.title = ['transcripts', 'gene_name']

        # parameters for rpf import
        self.all_rpf = pd.DataFrame()
        self.rpf_seq = pd.DataFrame()
        self.all_rpf_dict = OrderedDict()
        self.rpf_gene_list = pd.DataFrame()
        self.total_rpf_num = pd.DataFrame()
        self.mean_aa_rpm = None
        self.peak_rpm = []

        # parameters for rpf normalise
        self.background = args.background
        self.scale = args.scale
        self.fill = args.fill
        self.correlation = args.corr
        self.threshold = args.min

        # parameters for enrichment filter
        self.size = args.size
        self.polyorder = args.k
        self.backFold = args.backFold
        self.enrich = args.enrich
        self.gaps = args.gaps
        self.proportion = args.proportion
        self.width = args.width
        self.collision = args.collision

        self.upstream = args.upstream
        self.downstream = args.downstream

        # output peak detection files
        self.all = args.all
        self.peak_out = args.output
        self.ratio_out = args.ratio
        
        self.peak_log_file = self.peak_out + '_peaks.log'
        self.peak_txt_file = self.peak_out + '_peaks.txt'
        self.peak_bed_file = self.peak_out + '_peaks.bed'
        self.peak_sig_bed_file = self.peak_out + '_sig_peaks.bed'
        self.peak_seq_file = self.peak_out + '_peaks_sequence.txt'
        # self.peak_rpm_file = self.peak_out + '_peaks_rpm.txt'
        self.peak_ratio_file = self.peak_out + '_peaks_ratio.txt'
        self.enrich_ratio_file = self.peak_out + '_enrich_ratio.txt'

        # self.peak_merge = pd.DataFrame()
        self.peak_merge = []
        self.bed_merge = pd.DataFrame()

        # ratio
        self.all_ratio = OrderedDict()
        self.ave_ratio = OrderedDict()

        # output scripts and figure
        self.figure_out = args.fig
        self.enrich_out = []
        self.matlab_merge = OrderedDict()
        self.fig_merge = OrderedDict()
        self.peak_scripts_file = self.peak_out + '_peaks_scripts.m'
        self.none_peak_scripts_file = self.peak_out + '_peaks_none_scripts.m'

    # check all the arguments
    def args_check(self):
        # check the input file
        if not os.path.exists(self.rpf_file):
            print('ERROR! The specified Input file does not exist!', flush=True)
            sys.exit()

        # check the total rpf file
        if self.total_rpf and not os.path.exists(self.total_rpf):
            print('ERROR! The specified Total RPFs file does not exist!', flush=True)
            sys.exit()

        # check the gene annotation file
        if self.anno and not os.path.exists(self.anno):
            print('ERROR! The specified gene annotation file does not exist!', flush=True)
            sys.exit()

        # check the output file
        if os.path.exists(self.peak_log_file):
            print('Warning. The output file Peak_log exists and will be overwritten.', flush=True)
            os.remove(self.peak_log_file)

        if os.path.exists(self.peak_txt_file):
            print('Warning. The output file Peak_Annotation exists and will be overwritten.', flush=True)
            os.remove(self.peak_txt_file)

        if os.path.exists(self.peak_scripts_file):
            print('Warning. The output file of Peak_Scripts exists and will be overwritten.', flush=True)
            os.remove(self.peak_scripts_file)

        if os.path.exists(self.none_peak_scripts_file):
            print('Warning. The output file of None_Peak_Scripts exists and will be overwritten.', flush=True)
            os.remove(self.none_peak_scripts_file)

        if os.path.exists(self.peak_bed_file):
            print('Warning. The output file of Bed_Annotation exists and will be overwritten.', flush=True)
            os.remove(self.peak_bed_file)

        if os.path.exists(self.peak_ratio_file):
            print('Warning. The output file of Peak_Ratio exists and will be overwritten.', flush=True)
            os.remove(self.peak_ratio_file)

        # check the figure output path
        if self.figure_out:
            figure_out_path = self.peak_out + '_figures'
            if not os.path.exists(figure_out_path):
                os.makedirs(figure_out_path)
            else:
                print('Warning. The output directory already exists, and the figure will be output to '
                      'the existing directory.', flush=True)

    # function for gene annotation txt import
    def gene_anno(self):

        # if the annotation file is provided
        if self.anno:
            with open(self.anno, 'r') as anno_in:
                for line in anno_in:
                    rec = line.strip().split('\t')
                    self.gene_dict[rec[2]] = rec[1]
                # add the mrna not in the annotation file
                for mrna in self.gene_list:
                    if self.gene_dict.get(mrna):
                        pass
                    else:
                        self.gene_dict[mrna] = '-'
        else:
            for mrna in self.gene_list:
                self.gene_dict[mrna] = '-'

    # function for rpf txt import
    def rpf_txt_read(self):
        # get the input the samples name
        # samples = ck.split(',') + ip.split(',')
        # sample_num = len(ck_name + ip_name)
        col_index = ['name', 'from_tis', 'region']

        # import the data with polars
        all_rpf = pl.read_csv(self.rpf_file, separator='\t')
        # all_rpf = pd.read_csv(self.rpf_file, sep='\t', header=0, names=None)
        all_rpf = all_rpf.to_pandas()

        self.rpf_seq = all_rpf.iloc[:, 0:6]
        self.all_rpf = all_rpf.loc[:, col_index]

        for now_sp in self.ck_name + self.ip_name:
            now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
            self.all_rpf[now_sp] = all_rpf.loc[:, now_sp_index].sum(axis=1)

        col_title = self.all_rpf.columns
        self.gene_list = self.all_rpf.name.unique().tolist()
        rpf_group_gene = self.all_rpf.groupby('name')

        # Create the gene library.
        print('Create the gene library.', flush=True)
        # sys.stdout.flush()
        for mrna, grouped in rpf_group_gene:
            self.all_rpf_dict[mrna] = grouped

        # for mrna in self.gene_list:
        #     self.all_rpf_dict[mrna] = {}

        # with open(self.rpf_file, 'r') as cover_in:
        #     for line in islice(cover_in, 2, None):
        #         rec = line.strip().split('\t')
        #         self.all_rpf_dict[rec[0]][int(rec[3])] = rec
        #
        #     for mrna, mess in self.all_rpf_dict.items():
        #         self.all_rpf_dict[mrna] = pd.DataFrame.from_dict(mess, orient='index', columns=col_title)

        # import the total_rpf txt
        # the total rpf num for the RPM normalization
        if self.total_rpf:
            self.total_rpf_num = pd.read_csv(self.total_rpf, sep='\t', header=None, names=None, index_col=0)
            self.total_rpf_num = pd.Series(self.total_rpf_num[1], index=self.ck_name + self.ip_name)
            self.total_rpf_num = self.total_rpf_num.astype(int)
        else:
            self.total_rpf_num = self.all_rpf[self.ck_name + self.ip_name].sum()

        print('Total RPFs number:', flush=True)
        print(self.total_rpf_num, flush=True)

        # Calculate the mean value of CDS rpm for all genes and use it for background.
        ck_cds_rpf = pd.DataFrame(self.all_rpf[self.all_rpf.region == 'cds'], columns=col_index + self.ck_name)
        ck_cds_rpm = ck_cds_rpf[self.ck_name].astype('float').div(self.total_rpf_num.loc[self.ck_name].values).div(self.gene_list.__len__()) * self.scale
        self.mean_aa_rpm = ck_cds_rpm.sum() / self.background

    # define the format of output files
    def format_title(self):
        self.title = ['transcripts', 'gene_name']
        rpf_names = self.ck_name + self.ip_name
        rpm_names = [i + '_rpm' for i in rpf_names]
        results_title = ['min_corr_ck', 'min_corr_ip', 'comments', 'peak_num', 'peak_length', 'gaps', 'collision_start',
                         'peak_start', 'peak_end', 'collision_end', 'max_site', 'max_fold', 'mean_fold', 'ck_peak_rpm',
                         'ip_peak_rpm', 'P_Value', 'BHFDR', 'Peak_loci']
        self.title.extend(rpf_names + rpm_names + results_title)

    @staticmethod
    def safe_to_numeric(series):
        try:
            return pd.to_numeric(series)
        except ValueError:
            # ignore the error, and return the original series
            return series
    
    # define the format of rpf output list,
    def make_rpf_list(self, gene_rpf):
        gene_rpf_sum = gene_rpf[self.ck_name + self.ip_name].sum().tolist()
        # gene_rpf_mean = gene_rpf[ck_name + ip_name].mean().tolist()
        raw_gene_rpm = gene_rpf[self.ck_name + self.ip_name].div(self.total_rpf_num.values) * self.scale
        raw_gene_rpm_sum = raw_gene_rpm[self.ck_name + self.ip_name].sum().tolist()

        return gene_rpf_sum, raw_gene_rpm, raw_gene_rpm_sum

    # define the function for different protein binding modes.
    def flt_background_max_value(self, ip_name, raw_gene_rpm, cds_region):
        # condition 1: bound to ribosome
        if self.background == 0:
            return True
        # condition 2: bound to nascent polypeptides
        elif self.background == 30:

            max_back_value = raw_gene_rpm.loc[0:self.background - 1].max()
            # Exclude outliers of stop codon, "cds_region[-1]-1"
            max_bound_value = raw_gene_rpm.loc[self.background:cds_region[-1] - 1].max()

            # filter the strong binding score
            if self.backFold:
                delta_value = max_bound_value - max_back_value * self.enrich
            else:
                delta_value = max_bound_value - max_back_value * self.collision

            # flag-ip samples after background need 2 fold change higher than mock-ip
            if any(delta_value[ip_name] < 0):
                return False
            else:
                return True

    # define the function for different protein binding modes.
    def flt_background_mean_value(self, ip_name, raw_gene_rpm, cds_region):
        # condition 1: bound to ribosome
        if self.background == 0:
            return True
        
        # condition 2: bound to nascent polypeptides
        elif self.background == 30:
            
            mean_back_value = raw_gene_rpm.loc[0:self.background - 1].mean()
            # Exclude outliers of stop codon, "cds_region[-1]-1"
            mean_bound_value = raw_gene_rpm.loc[self.background:cds_region[-1] - 1].mean()

            # filter the strong binding score
            if self.backFold:
                delta_value = mean_bound_value - mean_back_value * self.enrich
            else:
                delta_value = mean_bound_value - mean_back_value * self.collision

            # flag-ip samples after background need 2 fold change higher than mock-ip
            if any(delta_value[ip_name] < 0):
                return False
            else:
                return True
            
    # define the function to calculate the correlation between the replicates
    def corr_calc(self, sp_gene_rpf_corr, sp_name):
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

    # define the function to fill the missing values with mean RPF
    def fill_empty(self, ck_name, raw_gene_rpm, cds_region):
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

    # define the function to calculate the enrichment ratio for each sample
    @staticmethod
    def calc_ratio(gene_rpm, ck_name, ip_name):
        gene_ratio = pd.DataFrame()
        # Ratio calculation is performed on any two samples in the two grouped data.
        for ip_sample in ip_name:
            for ck_sample in ck_name:
                gene_ratio[ip_sample + '_' + ck_sample] = gene_rpm[ip_sample].div(gene_rpm[ck_sample])

        return gene_ratio

    # define the function to smooth the data
    def smooth_data(self, gene_ratio, gene_rpf):
        # gene_ratio = gene_ratio.replace(np.inf, np.nan)
        gene_ratio_mean = gene_ratio.mean(axis=1)
        gene_ratio_mean.index = gene_rpf.from_tis

        if self.size == 0:
            gene_ratio_smooth = pd.DataFrame(gene_ratio_mean, index=gene_rpf.from_tis, columns=['enrich'])
            return gene_ratio_mean, gene_ratio_smooth
        else:
            gene_ratio_smooth = savgol_filter(gene_ratio_mean, self.size, self.k, mode='nearest')
            gene_ratio_smooth = pd.DataFrame(gene_ratio_smooth, index=gene_rpf.from_tis, columns=['enrich'])
            return gene_ratio_mean, gene_ratio_smooth

    # define the function to get the peak annotation
    @staticmethod
    def peak_anno(utr5_region, cds_region, utr3_region, peak_start, peak_end):
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

    # define the function to retrieve the collision region
    def retrieve_collision(self, gene_ratio_smooth, peak):
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

    # define the function to retrieve the maximum fold and position
    @staticmethod
    def get_max_enrich(peak, gene_ratio_smooth):
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

    # define the function to calculate the PValue with Wilcoxon-test (rank sum test)
    @staticmethod
    def calc_rank_sum_v1(posi_list, raw_gene_rpm, ck_name, ip_name):
        rank_sum_results = []
        peak_range = posi_list
        ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
        ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
        t, p = ranksums(ck_gene_rpf, ip_gene_rpf)
        rank_sum_results.append(str(p))

        return rank_sum_results

    # define the function to calculate the PValue with Wilcoxon-test (rank sum test)
    @staticmethod
    def calc_rank_sum_v2(peak, raw_gene_rpm, ck_name, ip_name, mrna):
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

    # define the function to calculate the PValue with paired t-test
    @staticmethod
    def calc_ttest_v1(posi_list, raw_gene_rpm, ck_name, ip_name):
        ttest_results = []
        peak_range = posi_list
        ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
        ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
        t, p = ttest_rel(ck_gene_rpf, ip_gene_rpf)
        ttest_results.append(str(p))

        return ttest_results

    # define the function to calculate the PValue with Wilcoxon-test (rank sum test)
    @staticmethod
    def calc_ttest_v2(peak, raw_gene_rpm, ck_name, ip_name, mrna):
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

    # define the function to find the peak
    def find_peak(self, eligible_index, gene_ratio_smooth, raw_gene_rpm, ck_name, ip_name):
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
                                ttest_results = self.calc_ttest_v1(temp_list[start_posi], raw_gene_rpm, ck_name,
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
                                ttest_results = self.calc_ttest_v1(temp_list[start_posi], raw_gene_rpm, ck_name,
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
                            if result and not(self.all):
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

    # define the function to save the matlab scripts for none peak results
    def matlab_none_peak_scripts(self, scripts_out, mrna, gene_ratio_smooth, cds_region):
        mrna_name = mrna.replace('-', '_').replace('.', '_')
        title_name = mrna_name.replace('_', '-')
        # ylim_down = gene_ratio_smooth.min().enrich - gene_ratio_smooth.min().enrich * 0.1
        # ylim_up = gene_ratio_smooth.max().enrich + gene_ratio_smooth.max().enrich * 0.1

        scripts_out.writelines(''.join("clear\n"))
        scripts_out.writelines(' '.join([mrna_name, '= [\n']))
        scripts_out.writelines(' '.join([str(i) for i in gene_ratio_smooth.index.to_list()]) + '\n')
        scripts_out.writelines(' '.join([str(round(float(i[0]), 4)) for i in gene_ratio_smooth.values.tolist()]) + '\n];\n')
        scripts_out.writelines(''.join("figure('Visible','off');") + '\n')
        scripts_out.writelines(''.join("set(gcf, 'units', 'normalized', 'position', [0.1 0.1 0.4 0.3]);") + '\n')
        scripts_out.writelines(''.join("hold on") + '\n')
        scripts_out.writelines(
            ''.join("plot(%s(1,:),%s(2,:),'Color','#0072bd','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
        scripts_out.writelines(''.join("yline(%s,'Color','#999999','LineStyle','-.');" % self.enrich) + '\n')
        scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[0]) + '\n')
        scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[-1]) + '\n')
        scripts_out.writelines(''.join("xlim([%s(1,1) %s(1,end)]);" % (mrna_name, mrna_name)) + '\n')
        # scripts_out.writelines(''.join("%%ylim([%s %s]);" % (math.floor(ylim_down), math.ceil(ylim_up))) + '\n')
        scripts_out.writelines(''.join("box on;") + '\n')
        scripts_out.writelines(''.join("legend({'un-bound'}, 'location', 'EO');") + '\n')
        scripts_out.writelines(''.join("xlabel('AA position on mRNA');") + '\n')
        scripts_out.writelines(''.join("ylabel('Enrichment [a.u]');") + '\n')
        scripts_out.writelines(''.join("title('%s (%s)');" % (title_name, self.gene_dict[mrna])) + '\n')
        scripts_out.writelines(''.join("set(gca, 'FontName', 'Arial', 'FontSize', 12)") + '\n')
        scripts_out.writelines(''.join("hold off") + '\n')
        scripts_out.writelines(''.join("print('%s', '-r300', '-dpng');" % mrna) + '\n')
        scripts_out.writelines(''.join("saveas(gcf, 'pdf', '%s');" % mrna) + '\n\n')

    # define the function to save the matlab scripts for peak results
    def matlab_peak_scripts(self, scripts_out, mrna, gene_ratio_smooth, peak, peak_left_start, peak_right_end,
                            cds_region):
        mrna_name = mrna.replace('-', '_').replace('.', '_')
        title_name = mrna_name.replace('_', '-')
        collision_region = []
        peak_region = []
        collision_ratio = gene_ratio_smooth.copy()
        peak_ratio = gene_ratio_smooth.copy()
        # ylim_down = gene_ratio_smooth.min().enrich - abs(gene_ratio_smooth.min().enrich * 0.1)
        # ylim_up = gene_ratio_smooth.max().enrich + abs(gene_ratio_smooth.max().enrich * 0.1)

        for peak_num, peak_dict in peak.items():
            peak_len = list(peak_dict.keys())[0]
            peak_start = str(peak[peak_num][peak_len][0])
            peak_end = str(peak[peak_num][peak_len][-1])
            collision_start = peak_left_start[peak_num]
            collision_end = peak_right_end[peak_num]
            collision_region.extend([i for i in range(collision_start, int(peak_start) + 1)] +
                                    [j for j in range(int(peak_end), collision_end + 1)])
            peak_region.extend([i for i in range(int(peak_start), int(peak_end) + 1)])

        collision_ratio.loc[collision_ratio.index.difference(collision_region)] = np.nan
        peak_ratio.loc[collision_ratio.index.difference(peak_region)] = np.nan

        scripts_out.writelines(''.join("clear\n"))
        scripts_out.writelines(' '.join([mrna_name, '= [\n']))
        scripts_out.writelines(' '.join([str(i) for i in gene_ratio_smooth.index.to_list()]) + '\n')
        scripts_out.writelines(' '.join([str(round(float(i[0]), 4)) for i in gene_ratio_smooth.values.tolist()]) + '\n')
        scripts_out.writelines(' '.join([str(round(float(i[0]), 4)) for i in collision_ratio.values.tolist()]) + '\n')
        scripts_out.writelines(' '.join([str(round(float(i[0]), 4)) for i in peak_ratio.values.tolist()]) + '\n];\n')
        scripts_out.writelines(''.join("figure('Visible','off');") + '\n')
        scripts_out.writelines(''.join("set(gcf, 'units', 'normalized', 'position', [0.1 0.1 0.4 0.3]);") + '\n')
        scripts_out.writelines(''.join("hold on") + '\n')
        scripts_out.writelines(
            ''.join("plot(%s(1,:),%s(2,:),'Color','#0072bd','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
        scripts_out.writelines(
            ''.join("plot(%s(1,:),%s(3,:),'Color','#edb120','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
        scripts_out.writelines(
            ''.join("plot(%s(1,:),%s(4,:),'Color','#d95319','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
        scripts_out.writelines(''.join("yline(%s,'Color','#999999','LineStyle','-.');" % self.enrich) + '\n')
        scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[0]) + '\n')
        scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[-1]) + '\n')
        scripts_out.writelines(''.join("xlim([%s(1,1) %s(1,end)]);" % (mrna_name, mrna_name)) + '\n')
        # scripts_out.writelines(''.join("%%ylim([%f %f]);" % (math.floor(ylim_down), math.ceil(ylim_up))) + '\n')
        scripts_out.writelines(''.join("box on;") + '\n')
        scripts_out.writelines(''.join("legend({'un-bound', 'edge', 'bound'}, 'location', 'EO');") + '\n')
        scripts_out.writelines(''.join("xlabel('AA position on mRNA');") + '\n')
        scripts_out.writelines(''.join("ylabel('Enrichment [a.u]');") + '\n')
        scripts_out.writelines(''.join("title('%s (%s)');" % (title_name, self.gene_dict[mrna])) + '\n')
        scripts_out.writelines(''.join("set(gca, 'FontName', 'Arial', 'FontSize', 10)") + '\n')
        scripts_out.writelines(''.join("hold off") + '\n')
        scripts_out.writelines(''.join("print('%s', '-r300', '-dpng');" % mrna) + '\n')
        scripts_out.writelines(''.join("save_pdf('%s');" % mrna) + '\n\n')

        return collision_ratio, peak_ratio

    # define the function to draw the figure
    def draw_figure(self, mrna, gene_ratio_smooth, collision_ratio, peak_ratio, cds_region):
        figure_out_path = './' + self.peak_out + '_figures'

        # draw the line plot
        mpl.rcParams['axes.linewidth'] = 0.64
        mpl.rcParams['xtick.labelsize'] = 6
        mpl.rcParams['ytick.labelsize'] = 6
        mpl.rcParams['legend.fontsize'] = 8

        fig = plt.figure(figsize=(5, 3), dpi=300)
        plt.subplots_adjust(0.1, 0.15, 0.9, 0.85)
        ax = fig.add_subplot(111)

        # gene_ratio.plot()
        plt.tick_params(axis='x', which='both',
                        bottom='on', top='off',
                        labelbottom='on',
                        direction='out', width=0.5, length=2.0)
        plt.tick_params(axis='y', which='both',
                        left='on', right='off',
                        labelleft='on', labelright='off',
                        direction='out', width=0.5, length=2.0)
        x = gene_ratio_smooth.index.tolist()
        plt.plot(x, gene_ratio_smooth.values.tolist(), color='#0072bd', linewidth=0.8, label='un-bound')
        plt.plot(x, collision_ratio.values.tolist(), color='#edb120', linewidth=0.8, label='bound')
        plt.plot(x, peak_ratio.values.tolist(), color='#d95319', linewidth=0.8, label='strongly-bound')
        plt.axvline(x=cds_region[0], c="#999999", ls="-.", lw=0.5)
        plt.axvline(x=cds_region[-1], c="#999999", ls="-.", lw=0.5)
        plt.axhline(y=self.enrich, c="#999999", ls="-.", lw=0.5)
        plt.legend()
        # plt.fill_between(x, gene_ratio.iloc[:, :].max(axis=1),
        #                 gene_ratio.iloc[:, :].min(axis=1), color='blue', linewidth=0.0, alpha=0.1)
        # plt.yscale('log', basey=2)
        ax.set_xlabel("AA position on mRNA", fontsize=8)
        ax.set_ylabel("Enrichment [a.u]", fontsize=8)
        ax.set_xlim(x[0], x[-1])
        ax.set_title(mrna + ' (' + self.gene_dict[mrna] + ')', fontsize=10)
        plt.savefig(figure_out_path + '/' + mrna + '.png')
        plt.close()

    # define the function to merge the output results
    def result_output(self, mrna, gene_rpf_sum, raw_gene_rpm_sum, ck_min_corr, ip_min_corr, peak,
                      peak_left_start, peak_right_end, max_peak_enrich,
                      max_peak_site, mean_peak_enrich, ck_rpm, ip_rpm, ttest_results, peak_results, utr5_region,
                      cds_region, utr3_region):
        
        for peak_num, peak_dict in peak.items():
            peak_len = list(peak_dict.keys())[0]
            peak_start = peak[peak_num][peak_len][0]
            peak_end = peak[peak_num][peak_len][-1]
            gaps = (peak_end - peak_start + 1) - peak_len
            # annotate the peak location in utr/cds
            peal_loci = self.peak_anno(utr5_region, cds_region, utr3_region, peak_start, peak_end)

            # save the peak results
            peak_tmp = '\t'.join([mrna, self.gene_dict[mrna], '\t'.join([str(i) for i in gene_rpf_sum]),
                                  '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                  ck_min_corr, ip_min_corr, 'Qualified', str(peak_num + 1),
                                  str(peak_len), str(gaps),
                                  str(peak_left_start[peak_num]),
                                  str(peak_start), str(peak_end),
                                  str(peak_right_end[peak_num]),
                                  str(max_peak_site[peak_num]),
                                  str(max_peak_enrich[peak_num]),
                                  mean_peak_enrich[peak_num],
                                  str(ck_rpm[peak_num]), str(ip_rpm[peak_num]),
                                  ttest_results[peak_num], '-',
                                  peal_loci])
            # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
            self.peak_merge.append(peak_tmp.split('\t'))

            peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')

        peak_results.flush()

    # define the function to detect binding peaks
    def detect_binding_peaks(self):
        # check and rename the output files
        peak_results = open(self.peak_log_file, 'w')

        # output the scripts of figures

        scripts_peak = open(self.peak_scripts_file, 'w')
        scripts_none_peak = open(self.none_peak_scripts_file, 'w')

        # get the samples name
        ck_name = self.ck.split(',')
        ip_name = self.ip.split(',')
        # col_index = ['transcript', 'from_tis', 'region']
        sample_num = len(ck_name + ip_name)

        # save the results
        self.format_title()
        peak_results.writelines('\t'.join(self.title) + '\n')

        # calc the ratio and find the peak
        for mrna in self.gene_list:
            # now gene is :
            print(mrna, flush=True)

            # get the gene rpf
            # gene_rpf = pd.DataFrame(self.all_rpf_dict[mrna], columns=col_index + self.ck_name + self.ip_name)
            gene_rpf = self.all_rpf_dict[mrna]
            # gene_rpf = gene_rpf.apply(pd.to_numeric, errors='ignore')
            gene_rpf = gene_rpf.apply(self.safe_to_numeric)
            gene_rpf.index = gene_rpf.from_tis

            # normalise the data to rpm
            gene_rpf_sum, raw_gene_rpm, raw_gene_rpm_sum = self.make_rpf_list(gene_rpf)

            # filter the rpf count >= threshold
            sp_rpf_threshold = 0
            for rpf_sum in gene_rpf_sum:
                if rpf_sum > self.threshold:
                    sp_rpf_threshold += 1

            if sp_rpf_threshold < sample_num:
                peak_tmp = '\t'.join([mrna, self.gene_dict[mrna],
                                      '\t'.join([str(i) for i in gene_rpf_sum]),
                                      '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                      '-', '-', 'Too few RPFs', '-', '-', '-', '-', '-', '-', '-',
                                      '-', '-', '-', '-', '-', '-', '-', '-'])
                # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
                self.peak_merge.append(peak_tmp.split('\t'))
                peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')
                continue

            # get the utr and cds region，, for the peak region annotation
            utr5_region = gene_rpf[gene_rpf.region == '5utr'].index.tolist()
            cds_region = gene_rpf[gene_rpf.region == 'cds'].index.tolist()
            utr3_region = gene_rpf[gene_rpf.region == '3utr'].index.tolist()

            # set the threshold of the background region
            first_aa = self.flt_background_max_value(ip_name, raw_gene_rpm, cds_region)
            # first_aa = self.flt_background_mean_value(ip_name, raw_gene_rpm, cds_region)

            if first_aa:
                pass
            else:
                peak_tmp = '\t'.join([mrna, self.gene_dict[mrna],
                                      '\t'.join([str(i) for i in gene_rpf_sum]),
                                      '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                      '-', '-', 'Background value too large', '-', '-', '-', '-', '-', '-',
                                      '-', '-', '-', '-', '-', '-', '-', '-', '-'])
                # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
                self.peak_merge.append(peak_tmp.split('\t'))
                peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')


                # use the mean_aa_rpm to fill the empty rpf.
                gene_rpm = self.fill_empty(ck_name, raw_gene_rpm, cds_region)

                # calculate the enrichment ratio for each sample
                gene_ratio = self.calc_ratio(gene_rpm, ck_name, ip_name)

                # smooth the data
                gene_ratio_mean, gene_ratio_smooth = self.smooth_data(gene_ratio, gene_rpf)

                # add the collision and peak ratio column with nan value
                enrich_out = pd.concat([gene_rpf, gene_ratio_smooth], axis=1).reset_index(drop=True)
                enrich_out = pd.concat([enrich_out, pd.DataFrame(np.nan, index=enrich_out.index, columns=['edge', 'bound'])], axis=1)

                self.enrich_out.append(enrich_out.reset_index(drop=True))

                continue

            # use the mean_aa_rpm to fill the empty rpf.
            gene_rpm = self.fill_empty(ck_name, raw_gene_rpm, cds_region)

            # calculate the enrichment ratio for each sample
            gene_ratio = self.calc_ratio(gene_rpm, ck_name, ip_name)

            # smooth the data
            gene_ratio_mean, gene_ratio_smooth = self.smooth_data(gene_ratio, gene_rpf)

            # combine the results of each gene
            self.all_ratio[mrna] = gene_ratio
            # ave_ratio[mrna] = gene_ratio_mean
            # cds_ratio[mrna] = gene_ratio.loc[cds_region]

            # calculate the correlation of samples
            ck_gene_rpf_corr = raw_gene_rpm[ck_name].corr()
            ck_min_corr = str(round(ck_gene_rpf_corr.values.min(), 4))
            
            ip_gene_rpf_corr = raw_gene_rpm[ip_name].corr()
            ip_min_corr = str(round(ip_gene_rpf_corr.values.min(), 4))

            # filter the correlation of the ck replicates： default 0.5
            ck_corr_pairs, total_pairs = self.corr_calc(ck_gene_rpf_corr, ck_name)

            if ck_corr_pairs < total_pairs:
                self.matlab_none_peak_scripts(scripts_none_peak, mrna, gene_ratio_smooth, cds_region)
                peak_tmp = '\t'.join([mrna, self.gene_dict[mrna],
                                      '\t'.join([str(i) for i in gene_rpf_sum]),
                                      '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                      ck_min_corr, ip_min_corr, 'Poor repeatability of ck samples',
                                      '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
                                      '-', '-'])
                
                # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
                self.peak_merge.append(peak_tmp.split('\t'))
                peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')
                continue

            # filter the correlation of the ip replicates： default 0.5
            ip_corr_pairs, total_pairs = self.corr_calc(ip_gene_rpf_corr, ip_name)

            if ip_corr_pairs < total_pairs:
                self.matlab_none_peak_scripts(scripts_none_peak, mrna, gene_ratio_smooth, cds_region)
                peak_tmp = '\t'.join([mrna, self.gene_dict[mrna],
                                      '\t'.join([str(i) for i in gene_rpf_sum]),
                                      '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                      ck_min_corr, ip_min_corr, 'Poor repeatability of ip samples',
                                      '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
                                      '-', '-'])
                
                # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
                self.peak_merge.append(peak_tmp.split('\t'))
                peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')
                continue

            # get the index of the rpm >= fold change,
            # can specify 1.2 for the collision and the 2 for the binding region
            eligible_index = gene_ratio_smooth[gene_ratio_smooth.enrich >= self.enrich].index.tolist()
            # index_12 = gene_ratio_smooth[gene_ratio_smooth.enrich >= 1.2].index.tolist()

            # filter the total aa number >= 5
            if len(eligible_index) < self.width:
                self.matlab_none_peak_scripts(scripts_none_peak, mrna, gene_ratio_smooth, cds_region)
                peak_tmp = '\t'.join([mrna, self.gene_dict[mrna],
                                      '\t'.join([str(i) for i in gene_rpf_sum]),
                                      '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                      ck_min_corr, ip_min_corr, 'Qualified',
                                      '0', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
                                      '-', '-'])
                # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
                self.peak_merge.append(peak_tmp.split('\t'))
                peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')
                continue

            # find_peak main programme is here !!!
            # get the peaks from the smooth data
            peak = self.find_peak(eligible_index, gene_ratio_smooth, raw_gene_rpm, ck_name, ip_name)

            # if no peak in current gene
            if not bool(peak):
                self.matlab_none_peak_scripts(scripts_none_peak, mrna, gene_ratio_smooth, cds_region)
                peak_tmp = '\t'.join([mrna, self.gene_dict[mrna],
                                      '\t'.join([str(i) for i in gene_rpf_sum]),
                                      '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                      ck_min_corr, ip_min_corr, 'Qualified',
                                      '0', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
                                      '-', '-'])
                # self.peak_merge = self.peak_merge.append(pd.Series(peak_tmp.split('\t')), ignore_index=True)
                self.peak_merge.append(peak_tmp.split('\t'))
                peak_results.writelines('\t'.join(peak_tmp.split('\t')) + '\n')

                enrich_out = pd.concat([gene_rpf, gene_ratio_smooth], axis=1).reset_index(drop=True)
                # add the collision and peak ratio column with nan value
                enrich_out = pd.concat([enrich_out, pd.DataFrame(np.nan, index=enrich_out.index, columns=['edge', 'bound'])], axis=1)

                self.enrich_out.append(enrich_out.reset_index(drop=True))

            else:
                # get the collision region from the smooth data
                peak_left_start, peak_right_end = self.retrieve_collision(gene_ratio_smooth, peak)

                # get the maximum site and fold change
                max_peak_enrich, max_peak_site, mean_peak_enrich = self.get_max_enrich(peak, gene_ratio_smooth)

                # calculate the p value of peaks by rank sum test
                ttest_results, ck_rpm, ip_rpm, rpm_dict = self.calc_ttest_v2(peak, raw_gene_rpm, ck_name, ip_name,
                                                                                   mrna)

                for peak_names, now_peak_rpm in rpm_dict.items():
                    mess = peak_names + '\t' + '\t'.join(now_peak_rpm) + '\n'
                    self.peak_rpm.append(mess)

                # output the peaks in txt format
                self.result_output(mrna, gene_rpf_sum, raw_gene_rpm_sum, ck_min_corr, ip_min_corr, peak,
                                   peak_left_start, peak_right_end,
                                   max_peak_enrich, max_peak_site, mean_peak_enrich, ck_rpm, ip_rpm, ttest_results,
                                   peak_results, utr5_region, cds_region, utr3_region)

                # output the figure scripts
                collision_ratio, peak_ratio = self.matlab_peak_scripts(scripts_peak, mrna, gene_ratio_smooth, peak,
                                                                       peak_left_start, peak_right_end, cds_region)

                # Output the demo figure
                if self.figure_out:
                    self.draw_figure(mrna, gene_ratio_smooth, collision_ratio, peak_ratio, cds_region)
                
                collision_ratio.rename(columns={'enrich': 'edge'}, inplace=True)
                peak_ratio.rename(columns={'enrich': 'bound'}, inplace=True)
                self.enrich_out.append(pd.concat([gene_rpf, gene_ratio_smooth, collision_ratio, peak_ratio], axis=1).reset_index(drop=True))

        self.peak_merge = pd.DataFrame(self.peak_merge)
        self.peak_merge.columns = self.title

        peak_results.close()
        scripts_none_peak.close()
        scripts_peak.close()

        print("transcripts number: %s" % len(self.gene_list), flush=True)

    # define the function to output the RPM of each peak
    def output_peak_rpm(self):
        # save the peak rpm
        split_line = pd.DataFrame(["##############################################"])
        peak_rpm = open(self.peak_rpm_file, 'w')
        # for mrna, ratio in self.all_ratio.items():
        #     ratio.index.name = mrna
        #     split_line.to_csv(self.peak_ratio_file, mode='a+', header=False, index=False)
        #     ratio.to_csv(self.peak_ratio_file, mode='a+', header=True, index=True)

    # define the function to output the original ratio data
    def output_ratio(self):
        split_line = pd.DataFrame(["##############################################"])
        for mrna, ratio in self.all_ratio.items():
            ratio.index.name = mrna
            split_line.to_csv(self.enrich_ratio_file, mode='a+', header=False, index=False)
            ratio.to_csv(self.enrich_ratio_file, mode='a+', sep = '\t', header=True, index=True)

    # define the function to output the peak message
    def output_peak(self):
        peak_all = self.peak_merge.loc[self.peak_merge.iloc[:, -3] != '-', ].copy()

        # peak_all.columns = self.title
        peak_all['P_Value'] = peak_all['P_Value'].astype('float')
        peak_all['rank'] = peak_all['P_Value'].rank(method='first')
        peak_num = peak_all['rank'].max()

        # calculate the FDR
        for idx, rows in peak_all.iterrows():
            peak_all.loc[idx, 'fdr'] = rows['P_Value'] * peak_num / rows['rank']

        # calculate the BHFDR
        peak_all = peak_all.sort_values('rank')
        peak_rank = peak_all['fdr'].to_list()
        for now_peak_num in range(int(peak_num)):
            now_fdr = peak_rank[now_peak_num]
            if now_peak_num + 1 < int(peak_num):
                now_min_fdr = min(peak_rank[now_peak_num + 1::])
                if now_fdr > now_min_fdr:
                    peak_rank[now_peak_num] = now_min_fdr
                else:
                    pass
            else:
                pass

        peak_all.loc[:, 'BHFDR'] = peak_rank
        del peak_all['fdr']
        peak_all = peak_all.sort_index()

        peak_all.to_csv(self.peak_txt_file, sep='\t', header=True, index=False)

        # save the results in bed format, start with 0
        # bed12 format: chrom / chromStart / chromEn / name / score / strand /
        # thickStart / thickEnd / itemRGB / blockCount / blockSize / blockStarts-
        peak_bed = peak_all.loc[:, ['transcripts', 'peak_start', 'peak_end', 'peak_num', 'mean_fold', 'max_site', 'P_Value']]
        peak_bed.loc[:, 'name'] = peak_bed['transcripts'] + '_peak_' + peak_bed['peak_num']
        peak_bed.loc[:, 'strand'] = '+'

        # peak_bed = peak_bed.apply(pd.to_numeric, errors='ignore')
        peak_bed= peak_bed.apply(self.safe_to_numeric)

        peak_bed.loc[:, 'peak_start'] = peak_bed.loc[:, 'peak_start'] - 1
        peak_bed = peak_bed.loc[peak_bed.loc[:, 'peak_start'] >= 0, :]
        peak_bed.loc[:, 'start'] = peak_bed.loc[:, 'peak_start'] * 3
        peak_bed.loc[:, 'end'] = peak_bed.loc[:, 'peak_end'] * 3 + 2

        bed_merge = peak_bed.loc[:, ['transcripts', 'start', 'end', 'name', 'mean_fold', 'strand']]
        bed_merge.columns = ['#transcripts', 'start', 'end', 'name', 'mean_fold', 'strand']
        bed_merge.to_csv(self.peak_bed_file, sep='\t', header=True, index=False)
        
        peak_sig_bed = peak_bed.loc[peak_bed.loc[:, 'P_Value'] < 0.05, :]
        bed_sig_merge = peak_sig_bed.loc[:, ['transcripts', 'start', 'end', 'name', 'mean_fold', 'strand']]
        bed_sig_merge.columns = ['#transcripts', 'start', 'end', 'name', 'mean_fold', 'strand']
        bed_sig_merge.to_csv(self.peak_sig_bed_file, sep='\t', header=True, index=False)

        # save the results in txt format contains the peak sequences
        # txt format: transcripts / gene_name / max_size / upstream sequence / peak sequence / downstream sequence

        peak_gene_name = peak_all.loc[:, 'transcripts'].unique().tolist()
        self.rpf_seq = self.rpf_seq.loc[self.rpf_seq.loc[:, 'name'].isin(peak_gene_name), :]
        rpf_seq_group = self.rpf_seq.groupby('name')

        peak_seq_list = []

        for idx, rows in peak_bed.iterrows():
            transcripts = rows['transcripts']
            gene_name = self.gene_dict[transcripts]
            peak_name = rows['name']

            up_start = rows['peak_start'] - self.upstream
            peak_start = rows['peak_start']
            peak_end = rows['peak_end']
            down_end = rows['peak_end'] + self.downstream

            gene_rpf = rpf_seq_group.get_group(transcripts)

            upstream_seq = gene_rpf.loc[(gene_rpf.from_tis >= up_start) & (gene_rpf.from_tis < peak_start), 'codon'].tolist()
            peak_seq = gene_rpf.loc[(gene_rpf.from_tis >= peak_start) & (gene_rpf.from_tis <= peak_end), 'codon'].tolist()
            downstream_seq = gene_rpf.loc[(gene_rpf.from_tis > peak_end) & (gene_rpf.from_tis <= down_end), 'codon'].tolist()
            max_site = rows['max_site']

            upstream_aa = Seq(''.join(upstream_seq)).translate()
            peak_aa = Seq(''.join(peak_seq)).translate()
            downstream_aa = Seq(''.join(downstream_seq)).translate()

            temp_peak_seq = [transcripts, gene_name, peak_name,
                             peak_start, peak_end, max_site,
                             ''.join(upstream_seq), ''.join(peak_seq), ''.join(downstream_seq),
                             str(upstream_aa), str(peak_aa), str(downstream_aa)]

            peak_seq_list.append(temp_peak_seq)
        
        self.peak_seq = pd.DataFrame(peak_seq_list, 
                                     columns=['transcripts', 'gene_name', 'peak_name',
                                              'peak_start', 'peak_end', 'max_site', 
                                              'upstream_nt', 'peak_nt', 'downstream_nt',
                                              'upstream_aa', 'peak_aa', 'downstream_aa'])
        
        self.peak_seq.to_csv(self.peak_seq_file, sep='\t', header=True, index=False)

        # # save the rpm of each peak
        # if self.rpm:
        #     self.output_peak_rpm()
        # else:
        #     pass

        # output the significant peaks enrichment ratio to self.peak_ratio_file
        peak_ratio_df = pd.concat(self.enrich_out, axis=0)
        peak_ratio_df.to_csv(self.peak_ratio_file, sep='\t', header=True, index=False)

        # save all the enrichment ratio of each gene
        if self.ratio_out:
            self.output_ratio()
        else:
            pass

    # define the function to output the matlab scripts
    def output_matlab(self):
        split_line = pd.DataFrame(["##############################################"])
        for mrna, ratio in self.all_ratio.items():
            ratio.index.name = mrna
            split_line.to_csv(self.peak_ratio_file, mode='a+', index=False, header=False)
            ratio.to_csv(self.peak_ratio_file, mode='a+', index=True, header=True)
