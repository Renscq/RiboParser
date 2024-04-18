#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/2/21 1:28
# @Project : ribo-pipe
# @Script  : find_peak.py
# @Version : python 3.8.5
# @product : PyCharm
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
import gc
import os
import re
import sys
import time
from collections import OrderedDict
from itertools import islice

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import ranksums


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(description="This script is used to summary the rpf and codon density.")
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.5"
    )
    # arguments for the input files
    input_group = parser.add_argument_group('Input files arguments')
    input_group.add_argument(
        "-i", dest="input", required=True, type=str, help="input the RPF coverage file"
    )
    input_group.add_argument(
        "-r", dest="rpf", required=False, type=str,
        help="specify the total RPFs number for normalise. (default: The sum of input samples)"
    )
    input_group.add_argument(
        "-a", dest="anno", required=False, type=str,
        help="specify the gene annotation file."
    )
    input_group.add_argument(
        "--ck", dest="control", required=True, type=str,
        help="specify the name of control samples, separated by commas."
    )
    input_group.add_argument(
        "--ip", dest="ip", required=True, type=str,
        help="specify the name of immunoprecipitation samples, separated by commas."
    )

    # arguments for the data filtering
    data_group = parser.add_argument_group('Data filtering arguments')
    data_group.add_argument(
        "-t", dest="threshold", required=False, type=int, default=10,
        help="specify the minimum number of RPFs coverage for genes to be retained (default %(default)s)"
    )
    data_group.add_argument(
        "--corr", dest="corr", required=False, type=float, default=0.5,
        help="specify the minimum correlation of replicates (default %(default)s)"
    )

    # arguments for the peaks scanning, gaps and ratio filtering
    peak_group = parser.add_argument_group('Peak scanning arguments')
    peak_group.add_argument(
        "-f", dest="fill", required=False, type=int, default=30,
        help="Use the mean value of gene RPF to fill the missing values. (Commonly used first "
             "30 AA, default %(default)s codon). Because the ribo-seq data is noisy in single"
             " gene level and many locations are not covered by RPF, it is necessary to fill "
             "the blank position use the mean value of RPF as the background. Three modes are "
             "provided: (1): [30] the first 30 AA (or other number); (2): [-1] the average RPFs"
             " value of current gene; (3): [0] the average RPF of all genes is used."
    )
    peak_group.add_argument(
        "-w", dest="width", required=False, type=int, default=5,
        help="specify the width of binding peaks (default %(default)s AA)."
             "By default, the width of each peak is not less than 5 amino acids."
    )
    peak_group.add_argument(
        "-e", dest="enrich", required=False, type=float, default=2.0,
        help="specify the enrichment fold of each peak height (default %(default)s)"
    )
    peak_group.add_argument(
        "-c", dest="collision", required=False, type=float, default=1.5,
        help="specify the fold of collision before/after the real binding peak (default %(default)s)"
    )
    peak_group.add_argument(
        "-g", dest="gaps", required=False, type=int, default=2,
        help="specify the gap width inside each peak. (default %(default)s AA)."
             "The ribo-seq profiling is noisy at the single gene level (for example,"
             "the ribosomal density along a gene may change from a positive number"
             "to zero and, again, to a positive number). Therefore, gaps in the peak"
             "are allowed. By default, the gaps there can be no more than 2 consecutive"
             "amino acids in a peak."
    )
    peak_group.add_argument(
        "-p", dest="proportion", required=False, type=float, default=0.15,
        help="specify the proportion of gap in each peak (default %(default)s)."
             "Data with too much noise has low credibility, so it is necessary to "
             "add a limit to the gap proportion in each peak. By default, the total gaps"
             " length within a peak cannot exceed 20%%."
    )
    peak_group.add_argument(
        "--back", dest="background", required=False, type=int, default=0,
        help="use the 5 prime AA to filter the enrichment fold (default %(default)s codon, "
             "not applicable). Two modes are provided: (1): [30](or other number) Some "
             "proteins that bind to the nascent polypeptide chain cannot affect the translation"
             " process of the first 30 AA. Theoretically, the RPFs of the unbound area should"
             " be smaller than the RPFs of the bound area. (2): [0] However, some proteins "
             "only affect the elongation of the ribosome, so the peak may exist in any location."
    )
    peak_group.add_argument(
        "--scale", dest="scale", required=False, type=float, default=1e6,
        help="specify the scale level for reads normalise (default %(default)s RPM)"
    )

    # arguments for output figures, ratio, peak annotation
    output_group = parser.add_argument_group('Output files arguments')
    output_group.add_argument(
        "-o", dest="output", required=False, type=str, default='results',
        help="prefix of output file name (default %(default)s_ratio.txt)."
             "This prefix will be used as the output of Peak/Ratio/Matlab Script/Bed files. "
    )
    output_group.add_argument(
        "-R", dest="ratio", required=False, type=str, default='F',
        help="Whether to output the original ratio, [T / F], (default %(default)s)"
    )
    output_group.add_argument(
        "--fig", dest="fig", required=False, type=str, default='F',
        help="whether to draw the demo graph of peak scan results, [T / F], (default %(default)s)")

    args = parser.parse_args()

    args_dict = vars(args)
    for k, v in args_dict.items():
        print("{:<12}: {:<}".format(k, str(v)))

    sys.stdout.flush()
    return args


# define all the global variables
# define the dictionary to save the ratio
all_ratio = OrderedDict()
ave_ratio = OrderedDict()


# check all the arguments
def args_check(rpf_file, total_rpf, gene_anno_file, peak_out, figure_out):
    # check the input file
    if not os.path.exists(rpf_file):
        print('ERROR! The specified Input file does not exist!')
        sys.exit()

    # check the total rpf file
    if total_rpf and not os.path.exists(total_rpf):
        print('ERROR! The specified Total RPFs file does not exist!')
        sys.exit()

    # check the gene annotation file
    if gene_anno_file and not os.path.exists(gene_anno_file):
        print('ERROR! The specified gene annotation file does not exist!')
        sys.exit()

    # check the output file
    peak_txt_file = peak_out + '_peaks.txt'
    if os.path.exists(peak_txt_file):
        print('Warning. The output file Peak_Annotation exists and will be overwritten.')

    peak_scripts_file = peak_out + 'peak_scripts.m'
    if os.path.exists(peak_scripts_file):
        print('Warning. The output file of Peak_Scripts exists and will be overwritten.')

    none_peak_scripts_file = peak_out + '_none_peak_scripts.m'
    if os.path.exists(none_peak_scripts_file):
        print('Warning. The output file of None_Peak_Scripts exists and will be overwritten.')

    peak_bed_file = peak_out + '_peaks.bed'
    if os.path.exists(peak_bed_file):
        print('Warning. The output file of Bed_Annotation exists and will be overwritten.')

    # check the figure output path
    if figure_out == "T":
        figure_out_path = peak_out + '_figures'
        if not os.path.exists(figure_out_path):
            os.makedirs(figure_out_path)
        else:
            print('Warning. The output directory already exists, and the figure will be output to the existing directory.')


# function for rpf txt import
def gene_anno(anno_file, gene_list):
    gene_dict = {}
    # if the annotation file is provided
    if anno_file:
        with open(anno_file, 'r') as anno_in:
            for line in anno_in:
                rec = line.strip().split('\t')
                gene_dict[rec[0]] = rec[1::]
            # add the mrna not in the annotation file
            for mrna in gene_list:
                if gene_dict.get(mrna):
                    pass
                else:
                    gene_dict[mrna] = '-'

            return gene_dict
    else:
        for mrna in gene_list:
            gene_dict[mrna] = '-'

        return gene_dict


# function for rpf txt import
def rpf_txt_read(rpf_file, total_rpf, ck, ip, background, scaled):
    # get the input the samples name
    ck_name = ck.split(',')
    ip_name = ip.split(',')
    # samples = ck.split(',') + ip.split(',')
    # sample_num = len(ck_name + ip_name)
    col_index = ['transcript', 'from_cds_start', 'region']

    # import the data
    all_rpf = pd.read_csv(rpf_file, sep='\t', header=0, names=None)
    col_title = all_rpf.columns
    gene_list = all_rpf.transcript.unique().tolist()

    # Create the gene library.
    print('Create the gene library.')
    sys.stdout.flush()

    all_rpf_dict = OrderedDict()
    for mrna in gene_list:
        all_rpf_dict[mrna] = {}

    with open(rpf_file, 'r') as cover_in:
        for line in islice(cover_in, 2, None):
            rec = line.strip().split('\t')
            all_rpf_dict[rec[0]][int(rec[3])] = rec

        for mrna, mess in all_rpf_dict.items():
            all_rpf_dict[mrna] = pd.DataFrame.from_dict(mess, orient='index', columns=col_title)

    # import the total_rpf txt
    if total_rpf:
        total_rpf_num = pd.read_csv(total_rpf, sep='\t', header=None, names=None, index_col=0)
        total_rpf_num = pd.Series(total_rpf_num[1], index=ck_name + ip_name)
    else:
        total_rpf_num = all_rpf[ck_name + ip_name].sum()

    print('Total RPFs number:')
    print(total_rpf_num)

    # Calculate the mean value of rpm for all genes and use it for background.
    ck_cds_rpf = pd.DataFrame(all_rpf[all_rpf.region == 'cds'], columns=col_index + ck_name)
    ck_cds_rpm = ck_cds_rpf[ck_name].astype('float').div(total_rpf_num.loc[ck_name].values).div(gene_list.__len__()) * scaled
    mean_aa_rpm = ck_cds_rpm.sum() / background

    return all_rpf_dict, gene_list, total_rpf_num, mean_aa_rpm


# define the format of output files
def format_title(ck_name, ip_name):
    title = ['transcripts', 'gene_name']
    rpf_names = ck_name + ip_name
    rpm_names = [i + '_rpm' for i in rpf_names]
    results_title = ['min_corr_ck', 'min_corr_ip', 'comments', 'peak_num', 'peak_length', 'gaps', 'collision_start',
                     'peak_start', 'peak_end', 'collision_end', 'max_site', 'max_fold', 'mean_fold', 'ck_peak_rpm',
                     'ip_peak_rpm', 'P_Value', 'Peak_loci']
    title.extend(rpf_names + rpm_names + results_title)
    return title


# define the format of rpf output list
def make_rpf_list(gene_rpf, ck_name, ip_name, total_rpf_num, scaled):
    gene_rpf_sum = gene_rpf[ck_name + ip_name].sum().tolist()
    # gene_rpf_mean = gene_rpf[ck_name + ip_name].mean().tolist()
    raw_gene_rpm = gene_rpf[ck_name + ip_name].div(total_rpf_num.values) * scaled
    raw_gene_rpm_sum = raw_gene_rpm[ck_name + ip_name].sum().tolist()

    return gene_rpf_sum, raw_gene_rpm, raw_gene_rpm_sum


# define the function for different protein binding modes.
def flt_background_value(ip_name, background, raw_gene_rpm, cds_region, enrich):
    # condition 1: bound to ribosome
    if background == 0:
        return True
    # condition 2: bound to nascent polypeptides
    elif background == 30:
        max_back_value = raw_gene_rpm.loc[0:background - 1].max()
        # Exclude outliers of stop codon, "cds_region[-1]-1"
        max_bound_value = raw_gene_rpm.loc[background:cds_region[-1] - 1].max()
        delta_value = max_bound_value - max_back_value * enrich

        if any(delta_value[ip_name] < 0):
            return False
        else:
            return True


# define the function to calculate the correlation between the replicates
def corr_calc(sp_gene_rpf_corr, sp_name, correlation):
    sp_corr_pairs = 0
    # more than one sample
    if len(sp_name) > 1:
        for rows in range(len(sp_name) - 1):
            for cols in range(rows, len(sp_name) - 1):
                if sp_gene_rpf_corr.iloc[rows, cols + 1] >= correlation:
                    sp_corr_pairs += 1
    # only one sample
    else:
        sp_corr_pairs = 1

    total_pairs = 0
    for counter in range(len(sp_name)):
        total_pairs += counter

    return sp_corr_pairs, total_pairs


# define the function to fill the missing values with mean RPF
def fill_empty(ck_name, raw_gene_rpm, mean_aa_rpm, fill, cds_region):
    gene_rpm = raw_gene_rpm.copy()
    gene_rpm_mean = gene_rpm.loc[cds_region].mean()

    for sample in ck_name:
        # condition 1: fill the missing values with mean RPFs of total genes
        if fill == 0:
            gene_rpm[sample].loc[gene_rpm[sample] == 0] = mean_aa_rpm[sample]
        # condition 2: fill the missing values with mean RPFs of currently gene cds
        elif fill == -1:
            gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_rpm_mean[sample]
        # condition 3: fill the missing values with mean RPFs of specific region
        else:
            # Get the background mean value of the specific region
            gene_background = gene_rpm[sample].loc[0:fill - 1]
            counter_aa = len(gene_background[gene_background == 0].index)
            gene_background_mean = gene_background.mean()
            # If two-thirds of the specific region is covered by RPF.
            if counter_aa >= fill * 2 / 3:
                gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_background_mean
            # If half of the specific region is covered by RPF, and gene_rpm_mean <= gene_background_mean
            elif counter_aa >= fill / 2 and gene_rpm_mean[sample] <= gene_background_mean:
                gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_background_mean
            # if the specific region is noisy, use gene_rpm_mean instead
            elif counter_aa >= fill / 2 and gene_rpm_mean[sample] > gene_background_mean:
                gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_rpm_mean[sample]
            # if the specific region is noisy, use gene_rpm_mean instead
            elif counter_aa < fill / 2:
                gene_rpm[sample].loc[gene_rpm[sample] == 0] = gene_rpm_mean[sample]
            else:
                gene_rpm[sample].loc[gene_rpm[sample] == 0] = mean_aa_rpm[sample]

    return gene_rpm


# define the function to calculate the enrichment ratio for each sample
def calc_ratio(gene_rpm, ck_name, ip_name):
    gene_ratio = pd.DataFrame()
    # Ratio calculation is performed on any two samples in the two grouped data.
    for ip_sample in ip_name:
        for ck_sample in ck_name:
            gene_ratio[ip_sample + '_' + ck_sample] = gene_rpm[ip_sample].div(gene_rpm[ck_sample])

    return gene_ratio


# define the function to smooth the data
def smooth_data(gene_ratio, gene_rpf):
    # gene_ratio = gene_ratio.replace(np.inf, np.nan)
    gene_ratio_mean = gene_ratio.mean(axis=1)
    gene_ratio_smooth = savgol_filter(gene_ratio_mean, 3, 1, mode='nearest')
    # gene_ratio_smooth = savgol_filter(gene_ratio_mean, 5, 3, mode='nearest')
    gene_ratio_mean.index = gene_rpf.from_cds_start
    gene_ratio_smooth = pd.DataFrame(gene_ratio_smooth, index=gene_rpf.from_cds_start, columns=['enrich'])

    return gene_ratio_mean, gene_ratio_smooth


# define the function to get the peak annotation
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
def retrieve_collision(gene_ratio_smooth, peak, collision):
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
                elif shift_aa_enrich < collision:
                    break
                elif shift_aa_enrich >= collision:
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
                elif shift_aa_enrich < collision:
                    break
                elif shift_aa_enrich >= collision:
                    peak_right_num += 1

        peak_left_start.append(peak_regions[0] - peak_left_num)
        # peak_left_end.append(peak_regions[0])
        # peak_right_start.append(peak_regions[-1])
        peak_right_end.append(peak_regions[-1] + peak_right_num)

    return peak_left_start, peak_right_end


# define the function to retrieve the maximum fold and position
def get_max_enrich(peak, gene_ratio_smooth):
    max_peak_enrich = []
    max_peak_site = []
    mean_peak_enrich = []

    for peak_num, peak_dict in peak.items():
        peak_list = list(peak_dict.values())[0]
        gene_ratio_list = gene_ratio_smooth.loc[peak_list]

        max_enrich = round(gene_ratio_list.max()[0], 4)
        max_peak_enrich.append(max_enrich)

        max_site = str(gene_ratio_list.idxmax()[0])
        max_peak_site.append(max_site)

        mean_enrich = str(round(gene_ratio_list.mean().enrich, 4))
        mean_peak_enrich.append(mean_enrich)

    return max_peak_enrich, max_peak_site, mean_peak_enrich


# define the function to calculate the PValue with Wilcoxon-test (rank sum test)
def calc_rank_sum_v1(posi_list, raw_gene_rpm, ck_name, ip_name):
    rank_sum_results = []
    peak_range = posi_list
    ck_gene_rpf = raw_gene_rpm.loc[peak_range][ck_name].mean(axis=1)
    ip_gene_rpf = raw_gene_rpm.loc[peak_range][ip_name].mean(axis=1)
    t, p = ranksums(ck_gene_rpf, ip_gene_rpf)
    rank_sum_results.append(str(p))

    return rank_sum_results


# define the function to calculate the PValue with Wilcoxon-test (rank sum test)
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
        rpm_dict[names] = list(map(str, raw_gene_rpm.loc[peak_range][ck_name+ip_name].sum(axis=0).to_list()))

    return rank_sum_results, ck_rpm, ip_rpm, rpm_dict


# define the function to find the peak
def find_peak(peak_width, collision, eligible_index, gaps_threshold, proportion, gene_ratio_smooth, raw_gene_rpm, ck_name, ip_name):
    peak = OrderedDict()
    pre_peak = OrderedDict()
    temp_store = []
    pre_peak_num = 1
    now_length = 0

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
            elif now_gap_length <= gaps_threshold:
                gap_fold_ratio = gene_ratio_smooth.loc[eligible_index[posi - 1] + 1:eligible_index[posi] - 1] >= collision
                # The loop continues when the enrichment in the gap higher than Collision.
                if False not in gap_fold_ratio.values:
                    temp_store.append(eligible_index[posi])
                    now_length = now_length + now_gap_length + 1

                # Save the results and break the current loop,
                # when the enrichment in the gap lower than Collision and the current peak width is wider than the peak_width.
                elif False in gap_fold_ratio.values and now_length >= peak_width:
                    pre_peak[pre_peak_num] = temp_store
                    pre_peak_num += 1
                    temp_store = [eligible_index[posi]]
                    now_length = 1

                else:
                    temp_store = [eligible_index[posi]]
                    now_length = 1

            # If the length of the consecutive gap is greater than the gap_threshold
            elif now_gap_length > gaps_threshold:
                if now_length >= peak_width:
                    pre_peak[pre_peak_num] = temp_store
                    pre_peak_num += 1
                    temp_store = [eligible_index[posi]]
                    now_length = 1

                elif now_length < peak_width:
                    temp_store = [eligible_index[posi]]
                    now_length = 1

    if now_length >= peak_width:
        pre_peak[pre_peak_num] = temp_store

    # Traverse each potential peak.
    now_peak_num = 0
    for peak_num, peak_posi in pre_peak.items():
        # if the current peak meets the conditions.
        if len(peak_posi) / (peak_posi[-1] - peak_posi[0] + 1) >= 1 - proportion:
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
                if len(temp_list[start_posi]) >= peak_width:
                    if len(temp_list[start_posi]) / (temp_list[start_posi][0] - temp_list[start_posi][-1] + 1) >= 1 - proportion:
                        rank_sum_results = calc_rank_sum_v1(temp_list[start_posi], raw_gene_rpm, ck_name, ip_name)
                        permutations[str(start_posi)] = [rank_sum_results, len(temp_list[start_posi]), start_posi, temp_list[start_posi]]
                    else:
                        pass
                for stop_posi in range(start_posi + 1, region_num):
                    merge_region = []
                    for i in temp_list[start_posi:stop_posi + 1]:
                        merge_region.extend(i)

                    merge_len = int(merge_region[-1]) - int(merge_region[0]) + 1
                    # If the current permutation meets the conditions.
                    if len(merge_region) / merge_len >= 1 - proportion:
                        permutations_num = '_'.join([str(num) for num in range(start_posi, stop_posi + 1)])
                        rank_sum_results = calc_rank_sum_v1(temp_list[start_posi], raw_gene_rpm, ck_name, ip_name)
                        permutations[permutations_num] = [rank_sum_results, merge_len, permutations_num, merge_region]
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
                    str_regex = re.compile(r'_')
                    result = [i for i in permutations.keys() if str_regex.match(i)]
                    # If there are permutations elements in the list
                    if result:
                        now_best_permutations = []
                        for permutations_num, permutations_info in permutations.items():
                            if len(now_best_permutations) == 0:
                                now_best_permutations = permutations_info
                            # if the length of current permutation large than the optimal one.
                            elif permutations_info[1] > now_best_permutations:
                                now_best_permutations = permutations_info
                            # if the length of current permutation equal to the optimal one.
                            elif permutations_info[1] == now_best_permutations[1]:
                                if permutations_info[0] < now_best_permutations[0]:
                                    now_best_permutations = permutations_info
                                # Compare the p_value of the current permutation and the optimal one.
                                elif permutations_info[0] == now_best_permutations[0]:
                                    # Compare the start position of the current permutation and the optimal one.
                                    if int(permutations_num.split('_')[0]) < int(now_best_permutations[2].split('_')[0]):
                                        now_best_permutations = permutations_info
                                elif permutations_info[0] > now_best_permutations[0]:
                                    pass
                            # If the width of the current permutation is narrower than the optimal one.
                            elif permutations_info[1] < now_best_permutations:
                                pass

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
                                permutations_num = '_'.join([str(num) for num in range(start_posi, stop_posi + 1)])
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
def matlab_none_peak_scripts(scripts_out, mrna, gene_dict, gene_ratio_smooth, enrich, cds_region):
    mrna_name = mrna.replace('-', '_').replace('.', '_')
    title_name = mrna_name.replace('_', '-')
    ylim_down = gene_ratio_smooth.min().enrich - gene_ratio_smooth.min().enrich * 0.1
    ylim_up = gene_ratio_smooth.max().enrich + gene_ratio_smooth.max().enrich * 0.1

    scripts_out.writelines(''.join("clear\n"))
    scripts_out.writelines(' '.join([mrna_name, '= [\n']))
    scripts_out.writelines(' '.join([str(i) for i in gene_ratio_smooth.index.to_list()]) + '\n')
    scripts_out.writelines(' '.join([str(i[0]) for i in gene_ratio_smooth.values.tolist()]) + '\n];\n')
    scripts_out.writelines(''.join("figure('Visible','off');") + '\n')
    scripts_out.writelines(''.join("set(gcf, 'units', 'normalized', 'position', [0.1 0.1 0.4 0.3]);") + '\n')
    scripts_out.writelines(''.join("hold on") + '\n')
    scripts_out.writelines(''.join("plot(%s(1,:),%s(2,:),'Color','#0072bd','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
    scripts_out.writelines(''.join("yline(%s,'Color','#999999','LineStyle','-.');" % enrich) + '\n')
    scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[0]) + '\n')
    scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[-1]) + '\n')
    scripts_out.writelines(''.join("xlim([%s(1,1) %s(1,end)]);" % (mrna_name, mrna_name)) + '\n')
    # scripts_out.writelines(''.join("%%ylim([%s %s]);" % (math.floor(ylim_down), math.ceil(ylim_up))) + '\n')
    scripts_out.writelines(''.join("box on;") + '\n')
    scripts_out.writelines(''.join("legend({'unbound'}, 'location', 'EO');") + '\n')
    scripts_out.writelines(''.join("xlabel('AA position on mRNA');") + '\n')
    scripts_out.writelines(''.join("ylabel('Enrichment [a.u]');") + '\n')
    scripts_out.writelines(''.join("title('%s (%s)');" % (title_name, gene_dict[mrna][1])) + '\n')
    scripts_out.writelines(''.join("set(gca, 'FontName', 'Arial', 'FontSize', 12)") + '\n')
    scripts_out.writelines(''.join("hold off") + '\n')
    scripts_out.writelines(''.join("print('%s', '-r300', '-dpng');" % mrna) + '\n')
    scripts_out.writelines(''.join("save_pdf('%s');" % mrna) + '\n\n')


# define the function to save the matlab scripts for peak results
def matlab_peak_scripts(scripts_out, mrna, gene_dict, gene_ratio_smooth, peak, peak_left_start, peak_right_end, enrich, cds_region):
    mrna_name = mrna.replace('-', '_').replace('.', '_')
    title_name = mrna_name.replace('_', '-')
    collision_region = []
    peak_region = []
    collision_ratio = gene_ratio_smooth.copy()
    peak_ratio = gene_ratio_smooth.copy()
    ylim_down = gene_ratio_smooth.min().enrich - abs(gene_ratio_smooth.min().enrich * 0.1)
    ylim_up = gene_ratio_smooth.max().enrich + abs(gene_ratio_smooth.max().enrich * 0.1)

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
    scripts_out.writelines(' '.join([str(i[0]) for i in gene_ratio_smooth.values.tolist()]) + '\n')
    scripts_out.writelines(' '.join([str(i[0]) for i in collision_ratio.values.tolist()]) + '\n')
    scripts_out.writelines(' '.join([str(i[0]) for i in peak_ratio.values.tolist()]) + '\n];\n')
    scripts_out.writelines(''.join("figure('Visible','off');") + '\n')
    scripts_out.writelines(''.join("set(gcf, 'units', 'normalized', 'position', [0.1 0.1 0.4 0.3]);") + '\n')
    scripts_out.writelines(''.join("hold on") + '\n')
    scripts_out.writelines(''.join("plot(%s(1,:),%s(2,:),'Color','#0072bd','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
    scripts_out.writelines(''.join("plot(%s(1,:),%s(3,:),'Color','#edb120','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
    scripts_out.writelines(''.join("plot(%s(1,:),%s(4,:),'Color','#d95319','LineWidth', 0.8);" % (mrna_name, mrna_name)) + '\n')
    scripts_out.writelines(''.join("yline(%s,'Color','#999999','LineStyle','-.');" % enrich) + '\n')
    scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[0]) + '\n')
    scripts_out.writelines(''.join("xline(%s,'Color','#999999','LineStyle','-.');" % cds_region[-1]) + '\n')
    scripts_out.writelines(''.join("xlim([%s(1,1) %s(1,end)]);" % (mrna_name, mrna_name)) + '\n')
    # scripts_out.writelines(''.join("%%ylim([%f %f]);" % (math.floor(ylim_down), math.ceil(ylim_up))) + '\n')
    scripts_out.writelines(''.join("box on;") + '\n')
    scripts_out.writelines(''.join("legend({'unbound', 'collision', 'bound'}, 'location', 'EO');") + '\n')
    scripts_out.writelines(''.join("xlabel('AA position on mRNA');") + '\n')
    scripts_out.writelines(''.join("ylabel('Enrichment [a.u]');") + '\n')
    scripts_out.writelines(''.join("title('%s (%s)');" % (title_name, gene_dict[mrna][1])) + '\n')
    scripts_out.writelines(''.join("set(gca, 'FontName', 'Arial', 'FontSize', 10)") + '\n')
    scripts_out.writelines(''.join("hold off") + '\n')
    scripts_out.writelines(''.join("print('%s', '-r300', '-dpng');" % mrna) + '\n')
    scripts_out.writelines(''.join("save_pdf('%s');" % mrna) + '\n\n')

    return collision_ratio, peak_ratio


# define the function to draw the figure
def draw_figure(peak_out, mrna, gene_dict, gene_ratio_smooth, collision_ratio, peak_ratio, enrich, cds_region):
    figure_out_path = './' + peak_out + '_figures'

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
    plt.plot(x, gene_ratio_smooth.values.tolist(), color='#0072bd', linewidth=0.8, label='unbound')
    plt.plot(x, collision_ratio.values.tolist(), color='#edb120', linewidth=0.8, label='collision')
    plt.plot(x, peak_ratio.values.tolist(), color='#d95319', linewidth=0.8, label='bound')
    plt.axvline(x=cds_region[0], c="#999999", ls="-.", lw=0.5)
    plt.axvline(x=cds_region[-1], c="#999999", ls="-.", lw=0.5)
    plt.axhline(y=enrich, c="#999999", ls="-.", lw=0.5)
    plt.legend()
    # plt.fill_between(x, gene_ratio.iloc[:, :].max(axis=1), gene_ratio.iloc[:, :].min(axis=1), color='blue', linewidth=0.0, alpha=0.1)
    # plt.yscale('log', basey=2)
    ax.set_xlabel("AA position on mRNA", fontsize=8)
    ax.set_ylabel("Enrichment [a.u]", fontsize=8)
    ax.set_xlim(x[0], x[-1])
    ax.set_title(mrna + ' (' + gene_dict[mrna][1] + ')', fontsize=10)
    plt.savefig(figure_out_path + '/' + mrna + '.png')
    plt.close()


# define the function to output the original ratio data
def output_ratio(ratio_out, peak_out):
    # whether to output the ratio
    if ratio_out == "T":
        print("Step5: Output the original ratio.")
        sys.stdout.flush()

        ratio_out_file = peak_out + '_ratio.txt'
        # Determine whether the file already exists.
        if os.path.exists(ratio_out_file):
            os.remove(ratio_out_file)

        split_line = pd.DataFrame(["##############################################"])
        for mrna, ratio in all_ratio.items():
            ratio.index.name = mrna
            split_line.to_csv(ratio_out_file, mode='a+', index=False, header=False)
            ratio.to_csv(ratio_out_file, mode='a+', index=True, header=True)
    else:
        pass


# define the function to output the results in txt format
def result_output(gene_dict, mrna, gene_rpf_sum, raw_gene_rpm_sum, ck_min_corr, ip_min_corr, peak, peak_left_start, peak_right_end, max_peak_enrich,
                  max_peak_site, mean_peak_enrich, ck_rpm, ip_rpm, rank_sum_results, peak_results, utr5_region, cds_region, utr3_region, mrna_peak_bed):
    for peak_num, peak_dict in peak.items():
        peak_len = list(peak_dict.keys())[0]
        peak_start = peak[peak_num][peak_len][0]
        peak_end = peak[peak_num][peak_len][-1]
        gaps = (peak_end - peak_start + 1) - peak_len
        # annotate the peak location in utr/cds
        peal_loci = peak_anno(utr5_region, cds_region, utr3_region, peak_start, peak_end)

        # save the peak results
        peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                           '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                           ck_min_corr, ip_min_corr, 'Qualified', str(peak_num + 1),
                                           str(peak_len), str(gaps),
                                           str(peak_left_start[peak_num]),
                                           str(peak_start), str(peak_end),
                                           str(peak_right_end[peak_num]),
                                           str(max_peak_enrich[peak_num]), max_peak_site[peak_num],
                                           mean_peak_enrich[peak_num],
                                           str(ck_rpm[peak_num]), str(ip_rpm[peak_num]),
                                           rank_sum_results[peak_num],
                                           peal_loci]) + '\n')

        # save the results in bed format
        # bed12 format: chrom/chromStart /chromEn /name /score /strand /thickStart /thickEnd /itemRGB /blockCount /blockSize /blockStarts-
        mrna_peak_bed.writelines(
            '\t'.join([mrna, str(peak_start), str(peak_end), mrna + '_peak_' + str(peak_num), mean_peak_enrich[peak_num]]) + '\n')

    peak_results.flush()
    mrna_peak_bed.flush()


# define the function to detect binding peaks
def detect_binding_peaks(gene_dict, ck, ip, peak_out, ratio, all_rpf, gene_list, total_rpf_num, background, fill, mean_aa_rpm, scaled, correlation,
                         threshold, enrich, gaps_threshold, proportion, peak_width, collision, figure_out):
    # check and rename the output files
    peak_out_file = peak_out + '_peaks.txt'
    peak_results = open(peak_out_file, 'w')

    script_peak_file = peak_out + '_peak_scripts.m'
    scripts_peak = open(script_peak_file, 'w')

    script_none_peak_file = peak_out + '_none_peak_scripts.m'
    scripts_none_peak = open(script_none_peak_file, 'w')

    mrna_peak_bed_file = peak_out + '_peaks.bed'
    mrna_peak_bed = open(mrna_peak_bed_file, 'w')

    # get the samples name
    ck_name = ck.split(',')
    ip_name = ip.split(',')
    col_index = ['transcript', 'from_cds_start', 'region']
    sample_num = len(ck_name + ip_name)

    # save the results
    title = format_title(ck_name, ip_name)
    peak_results.writelines('\t'.join(title) + '\n')

    # save the peak rpm
    peak_reads_file = peak_out + '_reads.txt'
    peak_reads = open(peak_reads_file, 'w')

    # calc the ratio and find the peak
    for mrna in gene_list:
        # now gene is :
        print(mrna)
        sys.stdout.flush()

        # get the gene rpf
        gene_rpf = pd.DataFrame(all_rpf[mrna], columns=col_index + ck_name + ip_name)
        gene_rpf = gene_rpf.apply(pd.to_numeric, errors='ignore')
        gene_rpf.index = gene_rpf.from_cds_start

        # normalise the data to rpm
        gene_rpf_sum, raw_gene_rpm, raw_gene_rpm_sum = make_rpf_list(gene_rpf, ck_name, ip_name, total_rpf_num, scaled)

        # filter the rpf count >= threshold
        sp_rpf_threshold = 0
        for rpf_sum in gene_rpf_sum:
            if rpf_sum > threshold:
                sp_rpf_threshold += 1
        if sp_rpf_threshold < sample_num:
            peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                               '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                               '-', '-', 'Too few RPFs', '-', '-', '-', '-', '-', '-', '-',
                                               '-', '-', '-', '-', '-', '-', '-']) + '\n')
            continue

        # get the utr and cds region，, for the peak region annotation
        utr5_region = gene_rpf[gene_rpf.region == '5utr'].index.tolist()
        cds_region = gene_rpf[gene_rpf.region == 'cds'].index.tolist()
        utr3_region = gene_rpf[gene_rpf.region == '3utr'].index.tolist()

        # set the threshold of the background region
        first_30aa = flt_background_value(ip_name, background, raw_gene_rpm, cds_region, enrich)
        if first_30aa:
            pass
        else:
            peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                               '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                               '-', '-', 'Background value too large', '-', '-', '-', '-', '-', '-', '-',
                                               '-', '-', '-', '-', '-', '-', '-']) + '\n')
            continue

        # use the mean_aa_rpm to fill the empty rpf.
        gene_rpm = fill_empty(ck_name, raw_gene_rpm, mean_aa_rpm, fill, cds_region)

        # calculate the enrichment ratio for each sample
        gene_ratio = calc_ratio(gene_rpm, ck_name, ip_name)

        # smooth the data
        gene_ratio_mean, gene_ratio_smooth = smooth_data(gene_ratio, gene_rpf)

        # combine the results of each gene
        all_ratio[mrna] = gene_ratio
        # ave_ratio[mrna] = gene_ratio_mean
        # cds_ratio[mrna] = gene_ratio.loc[cds_region]

        # calculate the correlation of samples
        ck_gene_rpf_corr = raw_gene_rpm[ck_name].corr()
        ck_min_corr = str(round(ck_gene_rpf_corr.min()[0], 4))
        ip_gene_rpf_corr = raw_gene_rpm[ip_name].corr()
        ip_min_corr = str(round(ip_gene_rpf_corr.min()[0], 4))

        # filter the correlation of the ck replicates： default 0.5
        ck_corr_pairs, total_pairs = corr_calc(ck_gene_rpf_corr, ck_name, correlation)
        if ck_corr_pairs < total_pairs:
            matlab_none_peak_scripts(scripts_none_peak, mrna, gene_dict, gene_ratio_smooth, enrich, cds_region)
            peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                               '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                               ck_min_corr, ip_min_corr, 'Poor repeatability of ck samples',
                                               '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']) + '\n')
            continue

        # filter the correlation of the ip replicates： default 0.5
        ip_corr_pairs, total_pairs = corr_calc(ip_gene_rpf_corr, ip_name, correlation)
        if ip_corr_pairs < total_pairs:
            matlab_none_peak_scripts(scripts_none_peak, mrna, gene_dict, gene_ratio_smooth, enrich, cds_region)
            peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                               '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                               ck_min_corr, ip_min_corr, 'Poor repeatability of ip samples',
                                               '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']) + '\n')
            continue

        # get the index of the rpm >= fold change,
        # 1.2 for the collision and the 2 for the stronger binding region
        eligible_index = gene_ratio_smooth[gene_ratio_smooth.enrich >= enrich].index.tolist()
        # index_12 = gene_ratio_smooth[gene_ratio_smooth.enrich >= 1.2].index.tolist()

        # filter the total aa number >= 5
        if len(eligible_index) < peak_width:
            matlab_none_peak_scripts(scripts_none_peak, mrna, gene_dict, gene_ratio_smooth, enrich, cds_region)
            peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                               '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                               ck_min_corr, ip_min_corr, 'Qualified',
                                               '0', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']) + '\n')
            continue

        # find_peak main programme is here !!!
        # get the peaks from the smooth data
        peak = find_peak(peak_width, collision, eligible_index, gaps_threshold, proportion, gene_ratio_smooth, raw_gene_rpm, ck_name, ip_name)

        # if no peak in current gene
        if not bool(peak):
            matlab_none_peak_scripts(scripts_none_peak, mrna, gene_dict, gene_ratio_smooth, enrich, cds_region)
            peak_results.writelines('\t'.join([mrna, gene_dict[mrna][1], '\t'.join([str(i) for i in gene_rpf_sum]),
                                               '\t'.join([str(round(i, 4)) for i in raw_gene_rpm_sum]),
                                               ck_min_corr, ip_min_corr, 'Qualified',
                                               '0', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']) + '\n')
        else:
            # get the collision region from the smooth data
            peak_left_start, peak_right_end = retrieve_collision(gene_ratio_smooth, peak, collision)

            # get the maximum site and fold change
            max_peak_enrich, max_peak_site, mean_peak_enrich = get_max_enrich(peak, gene_ratio_smooth)

            # calculate the p value of peaks by rank sum test
            rank_sum_results, ck_rpm, ip_rpm, rpm_dict = calc_rank_sum_v2(peak, raw_gene_rpm, ck_name, ip_name, mrna)

            for peak_names, peak_rpm in rpm_dict.items():
                peak_reads.writelines(peak_names + '\t' + '\t'.join(peak_rpm) + '\n')

            # output the peaks in txt format
            result_output(gene_dict, mrna, gene_rpf_sum, raw_gene_rpm_sum, ck_min_corr, ip_min_corr, peak, peak_left_start, peak_right_end,
                          max_peak_enrich, max_peak_site, mean_peak_enrich, ck_rpm, ip_rpm, rank_sum_results, peak_results, utr5_region, cds_region, utr3_region,
                          mrna_peak_bed)

            # output the figure scripts
            collision_ratio, peak_ratio = matlab_peak_scripts(scripts_peak, mrna, gene_dict, gene_ratio_smooth, peak, peak_left_start, peak_right_end,
                                                              enrich, cds_region)

            # Output the demo figure
            if figure_out == "T":
                draw_figure(peak_out, mrna, gene_dict, gene_ratio_smooth, collision_ratio, peak_ratio, enrich, cds_region)

    peak_results.close()
    scripts_none_peak.close()
    scripts_peak.close()
    peak_reads.close()

    # output the ratio
    output_ratio(ratio, peak_out)

    print("transcripts number: %s" % len(gene_list))


# main programme is here
def main():
    print("####################################################")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sys.stdout.flush()
    args = parse_args()

    print('\nStep1: Check all the files and directories.')
    args_check(args.input, args.rpf, args.anno, args.output, args.fig)

    print('\nStep2: Import the RPF data.')
    sys.stdout.flush()
    all_rpf, gene_list, total_rpf_num, mean_aa_rpm = rpf_txt_read(args.input, args.rpf, args.control, args.ip,
                                                                  args.background, args.scale)

    print('\nStep3: Import the gene annotation.')
    sys.stdout.flush()
    gene_dict = gene_anno(args.anno, gene_list)

    print('\nStep4: Scan peaks from the data.')
    sys.stdout.flush()
    detect_binding_peaks(gene_dict, args.control, args.ip, args.output, args.ratio, all_rpf, gene_list, total_rpf_num, args.background, args.fill,
                         mean_aa_rpm, args.scale, args.corr, args.threshold, args.enrich, args.gaps, args.proportion, args.width, args.collision,
                         args.fig)

    gc.collect()
    print("\nAll done!")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print("####################################################")


if __name__ == '__main__':
    main()

