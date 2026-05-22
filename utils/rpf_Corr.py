#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_corr.py


import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import polars as pl
import seaborn as sns

from utils.ribo import ArgsParser


def read_rpf(rpf_file):
    rpf = pl.read_csv(rpf_file, separator='\t', has_header=True)
    rpf = rpf.to_pandas()
    # cds = rpf[rpf["region"] == "cds"]
    sample_num = int((len(rpf.columns) - 6) / 3)
    sample_name = pd.Series(rpf.columns[6::]).str[:-3].drop_duplicates().tolist()
    frame_rpf = pd.DataFrame()
    frame0_rpf, frame1_rpf, frame2_rpf = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # split the frame
    for now_sp in sample_name:
        print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
        now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
        frame_rpf[now_sp] = rpf.loc[:, now_sp_index].sum(axis=1)
        frame0_rpf[now_sp] = rpf.loc[:, now_sp_index].iloc[:, 0]
        frame1_rpf[now_sp] = rpf.loc[:, now_sp_index].iloc[:, 1]
        frame2_rpf[now_sp] = rpf.loc[:, now_sp_index].iloc[:, 2]

    # calculate the rpm
    total_rpf = frame_rpf.sum()
    frame_rpm = frame_rpf.div(total_rpf) * 1e6
    frame0_rpm = frame0_rpf.div(total_rpf) * 1e6
    frame1_rpm = frame1_rpf.div(total_rpf) * 1e6
    frame2_rpm = frame2_rpf.div(total_rpf) * 1e6

    return rpf, sample_name, frame_rpm, frame0_rpm, frame1_rpm, frame2_rpm


def calc_gene_rpm(rpf, sample_name, frame_rpm, frame0_rpm, frame1_rpm, frame2_rpm, out_prefix):
    # add the annotation
    merge_f_rpm = pd.concat([rpf.iloc[:, 0:5], frame_rpm], axis=1)
    merge_f0_rpm = pd.concat([rpf.iloc[:, 0:5], frame0_rpm], axis=1)
    merge_f1_rpm = pd.concat([rpf.iloc[:, 0:5], frame1_rpm], axis=1)
    merge_f2_rpm = pd.concat([rpf.iloc[:, 0:5], frame2_rpm], axis=1)

    # summary the gene rpm
    cds_f_rpm = merge_f_rpm[merge_f_rpm["region"] == "cds"
                            ][["name"] + sample_name].groupby("name")[sample_name].apply(sum)
    cds_f0_rpm = merge_f0_rpm[merge_f0_rpm["region"] == "cds"
                              ][["name"] + sample_name].groupby("name")[sample_name].apply(sum)
    cds_f1_rpm = merge_f1_rpm[merge_f1_rpm["region"] == "cds"
                              ][["name"] + sample_name].groupby("name")[sample_name].apply(sum)
    cds_f2_rpm = merge_f2_rpm[merge_f2_rpm["region"] == "cds"
                              ][["name"] + sample_name].groupby("name")[sample_name].apply(sum)

    # calculate the correlation of samples at gene level
    gene_corr_f = cds_f_rpm.corr(method='pearson')
    gene_corr_f0 = cds_f0_rpm.corr(method='pearson')
    gene_corr_f1 = cds_f1_rpm.corr(method='pearson')
    gene_corr_f2 = cds_f2_rpm.corr(method='pearson')

    gene_corr_f.to_csv(out_prefix + '_gene_corr_frame.txt', sep='\t', index=True)
    gene_corr_f0.to_csv(out_prefix + '_gene_corr_f0.txt', sep='\t', index=True)
    gene_corr_f1.to_csv(out_prefix + '_gene_corr_f1.txt', sep='\t', index=True)
    gene_corr_f2.to_csv(out_prefix + '_gene_corr_f2.txt', sep='\t', index=True)

    return gene_corr_f, gene_corr_f0, gene_corr_f1, gene_corr_f2


def calc_rpf_rpm(frame_rpm, frame0_rpm, frame1_rpm, frame2_rpm, out_prefix):
    # calculate the correlation of samples at rpf level
    rpf_corr_f = frame_rpm.corr(method='pearson')
    rpf_corr_f0 = frame0_rpm.corr(method='pearson')
    rpf_corr_f1 = frame1_rpm.corr(method='pearson')
    rpf_corr_f2 = frame2_rpm.corr(method='pearson')

    rpf_corr_f.to_csv(out_prefix + '_rpf_corr_frame.txt', sep='\t', index=True)
    rpf_corr_f0.to_csv(out_prefix + '_rpf_corr_f0.txt', sep='\t', index=True)
    rpf_corr_f1.to_csv(out_prefix + '_rpf_corr_f1.txt', sep='\t', index=True)
    rpf_corr_f2.to_csv(out_prefix + '_rpf_corr_f2.txt', sep='\t', index=True)

    return rpf_corr_f, rpf_corr_f0, rpf_corr_f1, rpf_corr_f2


def draw_corr_plot(frame, frame0, frame1, frame2, types, out_prefix):
    # output figure name
    out_pdf = out_prefix + "_" + types + "_correlation_plot.pdf"
    out_png = out_prefix + "_" + types + "_correlation_plot.png"
    num = len(frame)

    # draw the figure
    matplotlib.use('AGG')

    fig, axes = plt.subplots(2, 2, figsize=(10 + num * 0.3, 9 + num * 0.25), dpi=300)
    ax1 = sns.heatmap(data=frame, cmap="Blues",  # vmin=-1, vmax=1, center=0,
                      xticklabels=1, yticklabels=1, ax=axes[0, 0])
    ax1.set_title('all frame rpf', fontsize=20)
    ax2 = sns.heatmap(data=frame0, cmap="Blues",  # vmin=-1, vmax=1, center=0,
                      robust=True, xticklabels=1, yticklabels=1, ax=axes[0, 1])
    ax2.set_title('frame0 rpf', fontsize=20)
    ax3 = sns.heatmap(data=frame1, cmap="Blues",  # vmin=-1, vmax=1, center=0,
                      robust=True, xticklabels=1, yticklabels=1, ax=axes[1, 0])
    ax3.set_title('frame1 rpf', fontsize=20)
    ax4 = sns.heatmap(data=frame2, cmap="Blues",  # vmin=-1, vmax=1, center=0,
                      robust=True, xticklabels=1, yticklabels=1, ax=axes[1, 1])
    ax4.set_title('frame2 rpf', fontsize=20)

    plt.tight_layout()
    # plt.show()
    fig.savefig(fname=out_pdf)
    fig.savefig(fname=out_png)


def main():
    ArgsParser.now_time()
    print('\nDraw the correlation of samples.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_corr_args_parser()

    print('\nStep2: Import the RPFs file.', flush=True)
    rpf, sample_name, frame_rpm, frame0_rpm, frame1_rpm, frame2_rpm = read_rpf(args.rpf)

    print('\nStep3: calculate the RPM.', flush=True)
    gene_corr_f, gene_corr_f0, gene_corr_f1, gene_corr_f2 = calc_gene_rpm(rpf,
                                                                          sample_name, frame_rpm, 
                                                                          frame0_rpm, frame1_rpm, frame2_rpm,
                                                                          args.output)
    rpf_corr_f, rpf_corr_f0, rpf_corr_f1, rpf_corr_f2 = calc_rpf_rpm(frame_rpm,
                                                                     frame0_rpm, frame1_rpm, frame2_rpm,
                                                                     args.output)

    print('\nStep4: Draw the correlation plot.', flush=True)
    draw_corr_plot(gene_corr_f, gene_corr_f0, gene_corr_f1, gene_corr_f2, 'gene', args.output)
    draw_corr_plot(rpf_corr_f, rpf_corr_f0, rpf_corr_f1, rpf_corr_f2, 'rpf', args.output)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
