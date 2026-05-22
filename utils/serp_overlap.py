#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : serp_overlap.py


import matplotlib
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
# from interval import Interval

from utils.ribo import ArgsParser


def import_peak_region(mock, flag):
    mock_peak = pd.read_csv(mock, sep='\t')
    sig_mock_peak = mock_peak[mock_peak['BHFDR'].replace('-', 'nan').astype(float) < 0.05].copy()

    flag_peak = pd.read_csv(flag, sep='\t')
    sig_flag_peak = flag_peak[flag_peak['BHFDR'].replace('-', 'nan').astype(float) < 0.05].copy()

    return sig_mock_peak, sig_flag_peak


# def get_overlap(left_peak, right_peak, marker):
#     for idx_mock, rows_mock in left_peak.iterrows():
#         if rows_mock.overlap != 'NaN':
#             continue
#         # get the same genes in another peak table
#         tmp_flag = right_peak[right_peak.transcripts == rows_mock.transcripts]
#         if tmp_flag.empty:
#             if left_peak.loc[idx_mock, 'overlap'] == 'NaN':
#                 left_peak.loc[idx_mock, 'overlap'] = marker
#         else:
#             # make the interval for test
#             mock_start = int(rows_mock.peak_start)
#             mock_end = int(rows_mock.peak_end)
#             mock_interval = Interval(mock_start, mock_end, lower_closed=True, upper_closed=True)
#             for idx_flag, rows_flag in tmp_flag.iterrows():
#                 # make the interval for test
#                 flag_start = int(rows_flag.peak_start)
#                 flag_end = int(rows_flag.peak_end)
#                 flag_interval = Interval(flag_start, flag_end, lower_closed=True, upper_closed=True)
#                 # annotate the overlap condition
#                 if flag_interval.overlaps(mock_interval):
#                     if left_peak.loc[idx_mock, 'overlap'] == 'NaN':
#                         left_peak.loc[idx_mock, 'overlap'] = 'overlap'
#                 else:
#                     if left_peak.loc[idx_mock, 'overlap'] == 'NaN':
#                         left_peak.loc[idx_mock, 'overlap'] = marker
#     return left_peak


def get_overlap(left_peak, right_peak, marker):
    for idx_mock, rows_mock in left_peak.iterrows():
        if rows_mock.overlap != 'NaN':
            continue
        # get the same genes in another peak table
        tmp_flag = right_peak[right_peak.transcripts == rows_mock.transcripts]
        if tmp_flag.empty:
            if left_peak.loc[idx_mock, 'overlap'] == 'NaN':
                left_peak.loc[idx_mock, 'overlap'] = marker
        else:
            mock_start = int(rows_mock.peak_start)
            mock_end = int(rows_mock.peak_end)
            for idx_flag, rows_flag in tmp_flag.iterrows():
                flag_start = int(rows_flag.peak_start)
                flag_end = int(rows_flag.peak_end)

                # annotate the overlap condition
                if flag_start <= mock_end and mock_start <= flag_end:
                    if left_peak.loc[idx_mock, 'overlap'] == 'NaN':
                        left_peak.loc[idx_mock, 'overlap'] = 'overlap'
                else:
                    if left_peak.loc[idx_mock, 'overlap'] == 'NaN':
                        left_peak.loc[idx_mock, 'overlap'] = marker
    return left_peak


def merge_overlap(sig_mock_peak, sig_flag_peak, out_mock, out_flag):
    # add 'overlap' column in to the peak table
    sig_mock_peak.insert(len(sig_mock_peak.columns), 'overlap', 'NaN', allow_duplicates=False)
    sig_flag_peak.insert(len(sig_flag_peak.columns), 'overlap', 'NaN', allow_duplicates=False)
    # get the overlap region and annotation
    sig_mock_anno = get_overlap(sig_mock_peak, sig_flag_peak, 'mock-IP')
    sig_flag_anno = get_overlap(sig_flag_peak, sig_mock_peak, 'flag-IP')
    # output
    sig_mock_anno.to_csv(out_mock + '.overlap.txt', sep='\t', header=True, index=False)
    sig_flag_anno.to_csv(out_flag + '.overlap.txt', sep='\t', header=True, index=False)

    return sig_mock_anno, sig_flag_anno


def plot_venn2(sig_mock_anno, sig_flag_anno, out_mock, out_flag):
    sig_mock_gene = set(sig_mock_anno.transcripts)
    sig_flag_gene = set(sig_flag_anno.transcripts)

    sig_peak_num = pd.concat([sig_mock_anno.overlap, sig_flag_anno.overlap]).value_counts()
    # sometimes there will be multiple crossovers in the overlap region between two samples.
    # so use the average instead
    sig_peak_num.overlap = int(sig_peak_num.overlap / 2)

    pdf_file = out_mock + '_vs_' + out_flag + '_venn.pdf'
    png_file = out_mock + '_vs_' + out_flag + '_venn.png'

    matplotlib.use('AGG')
    fig = plt.figure(figsize=(8, 4), dpi=300)
    # draw the gene overlap figure
    ax1 = fig.add_subplot(1, 2, 1)
    venn2([sig_mock_gene, sig_flag_gene], set_labels=['mock-IP', 'flag-IP'])
    plt.title('gene overlap')
    # draw the peak region overlap figure
    ax2 = fig.add_subplot(1, 2, 2)
    venn2([sig_peak_num['mock-IP'], sig_peak_num['flag-IP'], sig_peak_num['overlap']], set_labels=['mock-IP', 'flag-IP'])
    plt.title('peak region overlap')
    # output the figure
    fig.tight_layout()
    # plt.show()
    fig.savefig(fname=pdf_file)
    fig.savefig(fname=png_file, dpi=300)


def main():
    ArgsParser.now_time()
    print('\nRetrieve the sequence of peak region.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.serp_overlap()

    print('\nStep2: Import the peak region.', flush=True)
    sig_mock_peak, sig_flag_peak = import_peak_region(args.mock, args.flag)

    print('\nStep3: Get the overlap region.', flush=True)
    if not args.out_mock:
        out_mock = args.mock[::-4]
    else:
        out_mock = args.out_mock
    if not args.out_flag:
        out_flag = args.flag[::-4]
    else:
        out_flag = args.out_flag
    sig_mock_anno, sig_flag_anno = merge_overlap(sig_mock_peak, sig_flag_peak, out_mock, out_flag)

    print('\nStep4: Draw the venn figure.', flush=True)
    plot_venn2(sig_mock_anno, sig_flag_anno, args.out_mock, args.out_flag)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
