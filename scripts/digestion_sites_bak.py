#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_digest.py


import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns
from matplotlib.ticker import NullFormatter

from foo.ArgsParser import *
from foo.Digestion import *


def get_digest_sites(ribo_attr):
    ribo_attr.pysam_input = pysam.AlignmentFile(ribo_attr.bam, "rb")
    psite_start = {i: [0, 0, 0] for i in range(20, 100 + 1)}
    psite_end = {i: [0, 0, 0] for i in range(20, 100 + 1)}
    # psite_start = {i: [0, 0, 0] for i in range(ribo_attr.min_length, ribo_attr.max_length + 1)}
    # psite_end = {i: [0, 0, 0] for i in range(ribo_attr.min_length, ribo_attr.max_length + 1)}

    for line in ribo_attr.pysam_input.fetch(until_eof=True):
        if line.reference_name in ribo_attr.mrna_dict:
            map_start, map_end = line.get_blocks()[0]
            read_length = line.infer_read_length()
            if ribo_attr.min_length <= read_length <= ribo_attr.max_length:
                frame_start = (map_start - ribo_attr.mrna_dict[line.reference_name].utr5_length) % 3
                frame_end = (map_end - 1 - ribo_attr.mrna_dict[line.reference_name].utr5_length) % 3
                psite_start[read_length][frame_start] += 1
                psite_end[read_length][frame_end] += 1

    psite_start = pd.DataFrame(psite_start).T
    psite_end = pd.DataFrame(psite_end).T

    return psite_start, psite_end


def output_digest_sites(psite_start, psite_end, ribo_attr):
    psite_merge = pd.concat([psite_start, psite_end], join="outer", axis=1)
    psite_merge.fillna(0, inplace=True)
    psite_merge.columns = ["5end_frame1", "5end_frame2", "5end_frame3", "3end_frame1", "3end_frame2", "3end_frame3"]
    psite_merge.index.name = "length"
    psite_merge.to_csv(ribo_attr.output + "_digestion_sites.txt", sep='\t', index=True)


def digestion_plot(ribo_attr, psite_start, psite_end, scaled, output_prefix):
    out_pdf = output_prefix + "_digestion_sites_plot.pdf"
    out_png = output_prefix + "_digestion_sites_plot.png"

    if scaled:
        psite_start_norm = psite_start.sub(psite_start.mean(axis=1), axis=0)
        psite_start_norm_row = psite_start_norm.div(psite_start_norm.std(axis=1), axis=0)
        psite_start_norm_row.fillna(0, inplace=True)

        psite_end_norm = psite_end.sub(psite_end.mean(axis=1), axis=0)
        psite_end_norm_row = psite_end_norm.div(psite_end_norm.std(axis=1), axis=0)
        psite_end_norm_row.fillna(0, inplace=True)

        psite_end_norm_row["length"] = psite_end_norm_row.index
        psite_start_norm_row["length"] = psite_start_norm_row.index

        now_cmap = "Purples"
    else:
        psite_start_norm_row = psite_start
        psite_end_norm_row = psite_end

        psite_start["length"] = psite_start.index
        psite_end["length"] = psite_end.index

        now_cmap = "RdYlBu_r"

    x_size = int((ribo_attr.max_length - ribo_attr.min_length) / 6)
    y_size = int((ribo_attr.max_length - ribo_attr.min_length) / 4)

    matplotlib.use('AGG')
    # draw the periodicity of codon

    fig = plt.figure(figsize=(x_size, y_size), dpi=300)
    ax1 = fig.add_subplot(121)
    h1 = sns.heatmap(data=psite_start_norm_row.loc[ribo_attr.min_length: ribo_attr.max_length],
                     annot=None, linewidths=.5, ax=ax1, cmap=now_cmap,
                     cbar_kws={"orientation": "horizontal"})
    label_y = ax1.get_yticklabels()
    plt.setp(label_y, rotation=0, horizontalalignment='right')
    label_x = ax1.get_xticklabels()
    plt.setp(label_x, rotation=0, horizontalalignment='right')
    ax1.set_xlabel('reading frame')
    ax1.set_ylabel('RPFs length')
    ax1.set_title("5'-end digestion")

    ax2 = fig.add_subplot(122)
    h2 = sns.heatmap(data=psite_end_norm_row.loc[ribo_attr.min_length: ribo_attr.max_length],
                     annot=None, linewidths=.5, ax=ax2, cmap=now_cmap,
                     cbar_kws={"orientation": "horizontal"})
    label_y = ax2.get_yticklabels()
    plt.setp(label_y, rotation=0, horizontalalignment='right')
    label_x = ax2.get_xticklabels()
    plt.setp(label_x, rotation=0, horizontalalignment='right')
    ax2.set_xlabel('reading frame')
    ax2.set_ylabel('RPFs length')
    ax2.set_title("3'-end digestion")
#######################################################
    nullfmt = NullFormatter()
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.75
    left_h = left + width + 0.02
    bottom_h = bottom + height + 0.02
    hist_x_height = 0.1
    hist_y_width = 0.2

    rect_heatmap = [left, bottom, width, height]
    rect_hist_x = [left, bottom_h, width, hist_x_height]
    rect_hist_y = [left_h, bottom, hist_y_width, height]

    fig = plt.figure(figsize=(4, 6), dpi=300)
    ax_heatmap = plt.axes(rect_heatmap)
    ax_hist_x = plt.axes(rect_hist_x)
    ax_hist_y = plt.axes(rect_hist_y)

    # cbar_kws = {"orientation": "horizontal"}
    cbar_ax = fig.add_axes([0.87, 0.87, 0.015, 0.1])
    sns.heatmap(data=psite_end_norm_row.loc[ribo_attr.min_length: ribo_attr.max_length],
                annot=None, linewidths=.5, ax=ax_heatmap, cmap="Blues", cbar_ax=cbar_ax,
                cbar_kws={"orientation": "horizontal"})
    label_y = ax_heatmap.get_yticklabels()
    plt.setp(label_y, rotation=0, horizontalalignment='right')

    ax_hist_x.xaxis.set_major_formatter(nullfmt)
    ax_hist_y.yaxis.set_major_formatter(nullfmt)

    # df1 = pd.melt(df, id_vars=['length'], value_vars=['0', '1', '2'])
    df1 = pd.DataFrame(psite_end.sum(axis=0)).T
    sns.barplot(data=df1, color="#018a5c", alpha=.6, ax=ax_hist_x)
    ax_hist_x.yaxis.set_ticks_position("right")
    # ax_hist_x.yaxis.set_label_position("right")
    # sns.axes_style(rc={'ytick.left': False, 'ytick.right': True})
    ax_hist_x.set_xticklabels([])

    df2 = pd.DataFrame(psite_end.sum(axis=1)).loc[ribo_attr.min_length: ribo_attr.max_length].T
    sns.barplot(data=df2, color="#018a5c", alpha=.6, ax=ax_hist_y, orient='h')
    ax_hist_y.set_yticklabels([])
    # sns.color_palette("light:#5A9", as_cmap=True)
    # plt.tight_layout()
    plt.show()
##################################################################
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax1 = plt.subplot2grid((5, 6), (0, 0), rowspan=1, colspan=2)  # 相当于格子分成3行3列,列跨度为3，行跨度为1
    sns.barplot(data=df1, color="#388e70", alpha=.6, ax=ax1)
    ax1.set_xticklabels([])
    ax1.yaxis.set_ticks_position("right")

    cbar_ax1 = fig.add_axes([0.43, 0.75, 0.01, 0.12])
    # ax11 = plt.subplot2grid((5, 6), (0, 2), rowspan=1, colspan=1)

    ax2 = plt.subplot2grid((5, 6), (1, 0), rowspan=4, colspan=2)
    sns.heatmap(data=psite_end_norm_row.loc[ribo_attr.min_length: ribo_attr.max_length],
                annot=None, linewidths=.5, ax=ax2, cmap="Blues", cbar_ax=cbar_ax1)
    label_y = ax2.get_yticklabels()
    plt.setp(label_y, rotation=0, horizontalalignment='right')

    ax3 = plt.subplot2grid((5, 6), (1, 2), rowspan=4, colspan=1)
    sns.barplot(data=df2, color="#388e70", alpha=.6, ax=ax3, orient='h')
    ax3.set_yticklabels([])

    ax4 = plt.subplot2grid((5, 6), (0, 3), rowspan=1, colspan=2)
    sns.barplot(data=df1, color="#388e70", alpha=.6, ax=ax4)
    ax4.set_xticklabels([])
    ax4.yaxis.set_ticks_position("right")

    cbar_ax2 = fig.add_axes([0.83, 0.75, 0.01, 0.12])
    ax5 = plt.subplot2grid((5, 6), (1, 3), rowspan=4, colspan=2)
    sns.heatmap(data=psite_end_norm_row.loc[ribo_attr.min_length: ribo_attr.max_length],
                annot=None, linewidths=.5, ax=ax5, cmap="Blues", cbar_ax=cbar_ax2)
    ax5.set_yticklabels([])

    ax6 = plt.subplot2grid((5, 6), (1, 5), rowspan=4, colspan=1)
    sns.barplot(data=df2, color="#388e70", alpha=.6, ax=ax6, orient='h')
    ax6.yaxis.set_ticks_position("right")
    # ax6.set_yticklabels([])

    # plt.tight_layout()
    plt.show()

    fig.savefig(fname=out_pdf)
    fig.savefig(fname=out_png)


def main():
    sys.stdout.writelines('\nDetect the digestion sites.\n')
    sys.stdout.writelines('Step1: Checking the input Arguments.\n')
    args = digestion_args_parser()
    ribo_attr = Ribo(args)

    sys.stdout.writelines('Step2: Import the annotation of transcripts.\n')
    ribo_attr.read_transcript()

    sys.stdout.writelines('Step3: Detect the digestion sites.\n')
    psite_start, psite_end = get_digest_sites(ribo_attr)

    sys.stdout.writelines('Step4: Output the digestion sites.\n')
    output_digest_sites(psite_start, psite_end, ribo_attr)

    sys.stdout.writelines('Step5: Draw the heatmap of digestion sites.\n')
    digestion_plot(ribo_attr, psite_start, psite_end, args.scale, args.output)


if __name__ == '__main__':
    main()
