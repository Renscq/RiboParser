#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_geneplot.py


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from utils.ribo import ArgsParser


def read_gene_list(gene_list):
    gene_csv = pd.read_csv(gene_list, header=None, names=None)
    gene_name = []
    for gene in gene_csv.values.tolist():
        gene_name.extend(gene)

    return gene_name


def read_rpf(rpf_file, gene_name):
    if isinstance(gene_name, str):
        gene_name = [gene_name]
    else:
        pass

    rpf = pd.read_csv(rpf_file, sep='\t', header=0, names=None)

    sample_num = int((len(rpf.columns) - 6) / 3)
    sample_name = pd.Series(rpf.columns[6::]).str[:-3].drop_duplicates().tolist()
    total_cds_num = len(rpf["region"] == "cds") * 3

    total_rpf = []
    for num in range(sample_num):
        total_rpf.append(rpf.iloc[:, [6 + num * 3, 7 + num * 3, 8 + num * 3]].sum().sum())

    gene_rpf = rpf[rpf["name"].isin(gene_name)]
    merged_rpf = pd.DataFrame()

    # flt_gene_rpf1.drop(flt_gene_rpf1.columns[[1, 3, 4, 5]], axis=1, inplace=True)

    for gene in gene_name:
        now_gene_rpf = gene_rpf[gene_rpf["name"] == gene]
        if now_gene_rpf.empty:
            print('{0} are not fount in RPFs file.'.format(gene), flush=True)
            continue
        merged_gene_rpf = pd.DataFrame()

        for num in range(sample_num):
            now_idx = [0, 2, 3] + [6 + num * 3, 7 + num * 3, 8 + num * 3]
            flt_gene_rpf = now_gene_rpf.iloc[:, now_idx]
            melt_gene_rpf = pd.melt(flt_gene_rpf, id_vars=["name", "from_tis", "from_tts"]).sort_values(
                by=["name", "from_tis", "from_tts", "variable"])
            melt_gene_rpf = melt_gene_rpf.reset_index(drop=True)
            del melt_gene_rpf["variable"]
            melt_gene_rpf.columns = ['name', 'from_tis', 'from_tts', sample_name[num]]

            if not merged_gene_rpf.empty:
                merged_gene_rpf = pd.concat([merged_gene_rpf, melt_gene_rpf.iloc[:, -1]], axis=1)
            else:
                merged_gene_rpf = melt_gene_rpf.copy()

        if not merged_rpf.empty:
            merged_rpf = pd.concat([merged_rpf, merged_gene_rpf])
        else:
            merged_rpf = merged_gene_rpf.copy()

    total_rpf = pd.DataFrame(total_rpf).T

    merged_rpm = merged_rpf.iloc[:, 0:3].copy()
    merged_dst = merged_rpf.iloc[:, 0:3].copy()

    merged_rpm[sample_name] = merged_rpf[sample_name].div(total_rpf.values) * 1e7
    merged_dst[sample_name] = merged_rpf[sample_name].div(total_rpf.values) * total_cds_num

    return sample_num, sample_name, merged_rpm, merged_dst


def re_build_list(now_gene_rpm):
    start_site, stop_site = now_gene_rpm["from_tis"].iloc[0], now_gene_rpm["from_tis"].iloc[-1]
    start_nt_site = start_site * 3
    stop_nt_site = stop_site * 3 + 2

    codon = stop_site - start_site
    cds_length = codon * 3 + 2

    mrna_index = [i for i in range(start_nt_site, stop_nt_site + 1)]
    now_gene_rpm.index = mrna_index

    stop_codon_site = now_gene_rpm[now_gene_rpm["from_tts"] == 0]["from_tis"].index.tolist()[-1]
    now_labels = [start_nt_site, 0, stop_codon_site, stop_nt_site]

    return now_labels, cds_length, now_gene_rpm


def draw_density_plot(log_type, sample_num, sample_name, merged_rpm, merged_dst, gene_name):
    if isinstance(gene_name, str):
        gene_name = [gene_name]
    else:
        pass

    for gene in gene_name:
        now_gene_rpm = merged_rpm[merged_rpm["name"] == gene]
        if now_gene_rpm.empty:
            continue
        print("Now gene: {gene}.".format(gene=str(gene)), flush=True)
        now_labels, cds_length, draw_gene_rpm = re_build_list(now_gene_rpm)

        out_pdf = str(gene) + "_density_plot.pdf"
        out_png = str(gene) + "_density_plot.png"

        # matplotlib.use('AGG')

        # draw the periodicity of codon
        fig = plt.figure(figsize=(9, 3 * sample_num), dpi=300)

        for fig_num in range(sample_num):
            ax1 = fig.add_subplot(sample_num, 1, fig_num + 1)
            ax1.fill_between(draw_gene_rpm.index, draw_gene_rpm[sample_name[fig_num]], color="#f86934")
            # sns.barplot(x=draw_gene_rpm.index, y=draw_gene_rpm[sample_name[fig_num]], color="#f86934")
            # sns.lineplot(x=draw_gene_rpm.index, y=draw_gene_rpm[sample_name[fig_num]], color="#f86934")
            sns.lineplot(x=draw_gene_rpm.index, y=0, color="w")
            ax1.set_xticks(now_labels)
            ax1.set_xticklabels(now_labels)
            ax1.set_xlim([now_labels[0], now_labels[-1]])
            # ax1.set_ylim([0, int(now_gene_rpm.max() * 1.1)])

            # add the gene annotation here, need to import the txt file
            ax1.set_ylabel('RPTM')
            ax1.set_xlabel('position (nt)')
            ax1.set_title(sample_name[fig_num])
            plt.xticks(rotation=90)

            if log_type == "Not":
                pass
            else:
                plt.yscale('log', base=int(log_type))

        plt.suptitle(str(gene))
        plt.tight_layout()
        plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()


def main():
    ArgsParser.now_time()
    print('\nDraw the gene RPFs plot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_gene_plot_args_parser()

    print('\nStep2: Import the gene name.', flush=True)

    if args.list:
        gene_name = read_gene_list(args.list)
    else:
        gene_name = args.gene

    print('\nStep3: Import the RPFs file.', flush=True)
    sample_num, sample_name, merged_rpm, merged_dst = read_rpf(args.rpf, gene_name)

    print('\nStep4: Draw the gene RPFs plot.', flush=True)
    draw_density_plot(args.log, sample_num, sample_name, merged_rpm, merged_dst, gene_name)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
