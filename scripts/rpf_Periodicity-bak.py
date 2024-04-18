#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_periodicity.py

from collections import OrderedDict

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from foo import ArgsParser


def output_meta(high_gene_sum, out_file):
    out_txt = out_file + "_periodicity.txt"
    high_gene_sum.to_csv(out_txt, sep='\t', index=True)


def read_rpf(rpf_file, tis, tts, min_rpf):
    rpf = pd.read_csv(rpf_file, sep='\t', header=0, names=None)
    sample_num = int((len(rpf.columns) - 6) / 3)
    sample_name = pd.Series(
        rpf.columns[6::]).str[:-3].drop_duplicates().tolist()
    sample_dict = OrderedDict()
    rpf = rpf[(tis <= rpf["from_tis"]) & (rpf["from_tts"] <= -tts)]
    cds_rpf = pd.DataFrame()

    # for num in range(sample_num):
    #     cds_rpf[sample_name[num]] = rpf.filter(like=sample_name[num]).sum(axis=1)
    for now_sp in sample_name:
        print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
        now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
        cds_rpf[now_sp] = rpf.loc[:, now_sp_index].sum(axis=1)

    cds_rpf_mer = pd.concat([rpf.iloc[:, 0:5], cds_rpf], axis=1)
    cds_rpf_sum = cds_rpf_mer.groupby("name")[sample_name].apply(sum)

    for now_sp in sample_name:
        # print('Now is: {file_name}.'.format(file_name=now_sp), flush=True)
        if min_rpf <= 0:
            high_gene = cds_rpf_sum[now_sp]
            now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
            high_gene_rpf = rpf.loc[:, now_sp_index].sum(axis=0)
        else:
            high_gene = cds_rpf_sum.loc[cds_rpf_sum[now_sp] > min_rpf, now_sp]
            now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
            high_gene_rpf = rpf[rpf["name"].isin(
                high_gene.index)].loc[:, now_sp_index].sum(axis=0)

        high_gene_sum = pd.concat(
            [high_gene_rpf,
             high_gene_rpf.div(high_gene_rpf.sum()) * 100],
            axis=1)
        high_gene_sum = high_gene_sum.rename(columns={0: "RPFs", 1: "ratio"})

        high_gene_sum.index.name = 'frame'
        output_meta(high_gene_sum, now_sp)
        sample_dict[now_sp] = [high_gene_sum, len(high_gene.index)]

    return sample_dict


def draw_periodicity_plot(high_gene_sum, gene_num, out_file):
    print('Now is: {fname}'.format(fname=out_file), flush=True)
    out_pdf = out_file + "_periodicity_plot.pdf"
    out_png = out_file + "_periodicity_plot.png"

    matplotlib.rc('figure', max_open_warning=0)
    matplotlib.use('AGG')
    # draw the periodicity of codon
    fig = plt.figure(figsize=(6, 4), dpi=300)
    ax1 = fig.add_subplot(121)
    bars = ax1.bar([0, 1, 2], high_gene_sum["RPFs"], align='center')
    ax1.set_xticks([0, 1, 2])
    ax1.set_xticklabels(["frame1", "frame2", "frame3"])
    # ax.invert_xaxis()
    ax1.set_ylim([0, int(max(high_gene_sum["RPFs"]) * 1.1)])
    ax1.set_ylabel('RPFs')
    ax1.set_xlabel('Codon')
    # ax1.set_title('Periodicity of codon', fontsize=20)
    for x_loc in range(3):
        y_loc = bars[x_loc].get_height()
        if y_loc > 50:
            ax1.text(x_loc,
                     y_loc,
                     '%.2f' % high_gene_sum["RPFs"][x_loc],
                     ha="center",
                     va='top')
        else:
            ax1.text(x_loc,
                     y_loc,
                     '%.2f' % high_gene_sum["RPFs"][x_loc],
                     ha="center",
                     va='bottom')
    # ax1.bar_label(bars, fmt="%.1e")
    plt.tight_layout()

    ax2 = fig.add_subplot(122)
    bars = ax2.bar([0, 1, 2], high_gene_sum["ratio"], align='center')
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(["frame1", "frame2", "frame3"])
    # ax.invert_xaxis()
    ax2.set_ylim([0, 100])
    ax2.set_ylabel('Ratio (%)')
    ax2.set_xlabel('Codon')
    # ax2.set_title('Periodicity of codon')
    for x_loc in range(3):
        y_loc = bars[x_loc].get_height()
        if y_loc > 50:
            ax2.text(x_loc,
                     y_loc,
                     '%.2f' % high_gene_sum["ratio"][x_loc],
                     ha="center",
                     va='top')
        else:
            ax2.text(x_loc,
                     y_loc,
                     '%.2f' % high_gene_sum["ratio"][x_loc],
                     ha="center",
                     va='bottom')
    # ax2.bar_label(bars, fmt='%.2f')
    plt.suptitle("3nt periodicity ({number} genes)".format(number=gene_num))
    plt.tight_layout()

    # plt.show()
    fig.savefig(fname=out_pdf)
    fig.savefig(fname=out_png)


def main():
    ArgsParser.now_time()
    print('Draw the periodicity plot.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.periodicity_args_parser()

    print('Step2: Import the RPFs file.\n', flush=True)
    sample_dict = read_rpf(args.rpf, args.tis, args.tts, args.min)

    print('Step3: Draw the periodicity plot.\n', flush=True)
    for sample_name, mess in sample_dict.items():
        high_gene_sum, gene_num = mess[0], mess[1]
        draw_periodicity_plot(high_gene_sum, gene_num, sample_name)

    print('All done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
