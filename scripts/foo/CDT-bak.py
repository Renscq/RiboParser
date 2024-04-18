#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : CDT.py


import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


class CodonDecodingTime(object):

    def __init__(self, args):
        self.gene = args.list
        self.rpf_file = args.rpf

        self.sample_num = 0
        self.sample_name = None
        self.rel_name = None
        self.merged_rpf = None
        self.merged_rpf_na = None

        # for the high expression gene filter
        self.rpf_num = args.min
        self.rpf = None
        self.gene_length = None
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None

        # filter the stable codon coverage range
        self.scale = args.scale
        self.site = args.site
        self.frame = args.frame
        self.tis = args.tis
        self.tts = args.tts

        # table for codon occupancy
        self.cdt = None
        self.codon_anno = {
            'AAA': ['Lys', 'K'], 'AAC': ['Asn', 'N'], 'AAG': ['Lys', 'K'], 'AAT': ['Asn', 'N'],
            'ACA': ['Thr', 'T'], 'ACC': ['Thr', 'T'], 'ACG': ['Thr', 'T'], 'ACT': ['Thr', 'T'],
            'AGA': ['Arg', 'R'], 'AGC': ['Ser', 'S'], 'AGG': ['Arg', 'R'], 'AGT': ['Ser', 'S'],
            'ATA': ['Ile', 'I'], 'ATC': ['Ile', 'I'], 'ATG': ['Met', 'M'], 'ATT': ['Ile', 'I'],
            'CAA': ['Gln', 'Q'], 'CAC': ['HIS', 'H'], 'CAG': ['Gln', 'Q'], 'CAT': ['HIS', 'H'],
            'CCA': ['Pro', 'P'], 'CCC': ['Pro', 'P'], 'CCG': ['Pro', 'P'], 'CCT': ['Pro', 'P'],
            'CGA': ['Arg', 'R'], 'CGC': ['Arg', 'R'], 'CGG': ['Arg', 'R'], 'CGT': ['Arg', 'R'],
            'CTA': ['Leu', 'L'], 'CTC': ['Leu', 'L'], 'CTG': ['Leu', 'L'], 'CTT': ['Leu', 'L'],
            'GAA': ['Glu', 'E'], 'GAC': ['Asp', 'D'], 'GAG': ['Glu', 'E'], 'GAT': ['Asp', 'D'],
            'GCA': ['Ala', 'A'], 'GCC': ['Ala', 'A'], 'GCG': ['Ala', 'A'], 'GCT': ['Ala', 'A'],
            'GGA': ['Gly', 'G'], 'GGC': ['Gly', 'G'], 'GGG': ['Gly', 'G'], 'GGT': ['Gly', 'G'],
            'GTA': ['Val', 'V'], 'GTC': ['Val', 'V'], 'GTG': ['Val', 'V'], 'GTT': ['Val', 'V'],
            'TAC': ['Tyr', 'Y'], 'TAT': ['Tyr', 'Y'], 'TCA': ['Ser', 'S'], 'TCC': ['Ser', 'S'],
            'TCG': ['Ser', 'S'], 'TCT': ['Ser', 'S'], 'TGC': ['Cys', 'C'], 'TGG': ['Trp', 'W'],
            'TGT': ['Cys', 'C'], 'TTA': ['Leu', 'L'], 'TTC': ['Phe', 'F'], 'TTG': ['Leu', 'L'],
            'TTT': ['Phe', 'F']
        }
        # except, 'TAA': ['*', '*'], 'TAG': ['*', '*'], 'TGA': ['*', '*']
        self.codon_table = None
        self.output = args.output

    def get_frame_rpf(self, rpf):
        if self.frame == 'all':
            for now_sp in self.sample_name:
                print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
                now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
                self.merged_rpf[now_sp] = rpf.loc[:, now_sp_index].sum(axis=1)
        elif self.frame == '0':
            for now_sp in self.sample_name:
                print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
                now_sp_index = [now_sp + '_f0']
                self.merged_rpf[now_sp] = rpf.loc[:, now_sp_index]
        elif self.frame == '1':
            for now_sp in self.sample_name:
                print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
                now_sp_index = [now_sp + '_f1']
                self.merged_rpf[now_sp] = rpf.loc[:, now_sp_index]
        elif self.frame == '2':
            for now_sp in self.sample_name:
                print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
                now_sp_index = [now_sp + '_f2']
                self.merged_rpf[now_sp] = rpf.loc[:, now_sp_index]

    def read_rpf(self):
        self.rpf = pd.read_csv(self.rpf_file, sep='\t', header=0, names=None)
        self.sample_num = int((len(self.rpf.columns) - 6) / 3)
        self.sample_name = pd.Series(self.rpf.columns[6::]).str[:-3].drop_duplicates().tolist()
        self.gene_length = self.rpf.loc[self.rpf['from_tts'] == 0,
                                        ['name', 'now_nt']].set_index('name', drop=True)
        self.gene_length.columns = ['length']

        # merge rpfs
        self.merged_rpf = self.rpf.iloc[:, 0:6].copy()
        self.frame_rpm = self.rpf.iloc[:, 0:6].copy()
        for now_sp in self.sample_name:
            print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
            all_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
            self.merged_rpf[now_sp] = self.rpf.loc[:, all_sp_index].sum(axis=1)

            total_rpf = self.merged_rpf[now_sp].sum()
            self.frame_rpm[all_sp_index] = self.rpf[all_sp_index] * 1e6 / total_rpf

        # self.get_frame_rpf(rpf)
        self.merged_rpf = self.merged_rpf[(self.merged_rpf["from_tis"] >= self.tis) &
                                          (self.merged_rpf["from_tts"] <= -self.tts)]
        # filter the specific genes
        if self.gene:
            gene_csv = pd.read_csv(self.gene, header=None, names=None)
            gene_name = []
            for gene in gene_csv.values.tolist():
                gene_name.extend(gene)
            self.merged_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(gene_name)]
        else:
            pass

        # filter the high expression gene
        self.total_rpf_num = self.merged_rpf[self.sample_name].sum()
        self.gene_rpf_sum = self.merged_rpf.groupby('name')[self.sample_name].sum()
        self.high_gene = self.gene_rpf_sum.loc[self.gene_rpf_sum.max(axis=1) >= self.rpf_num, :].index
        self.high_rpm = self.frame_rpm.loc[self.rpf['name'].isin(self.high_gene), :]

    @staticmethod
    def scale_method(scale, cdt):
        if scale == 'minmax':
            relative_cdt = (cdt - cdt.min()) / (cdt.max() - cdt.min())
        elif scale == 'zscore':
            relative_cdt = (cdt - cdt.mean()) / cdt.std()
        else:
            print('Unknown scale method.', flush=True)
            sys.exit()
        return relative_cdt
    

    def codon_decoding_time(self):

        # codon shift by E/P/A site
        if self.site == 'E':
            shift_codon = -1
        elif self.site == 'P':
            shift_codon = 0
        elif self.site == 'A':
            shift_codon = 1
        else:
            shift_codon = 0

        # calculate the average density of each codon
        high_rpm_group = self.high_rpm.groupby('name')
        now_num = 0
        ave_cdt_list = []
        sp_index = pd.Series(self.high_rpm.columns[6::])
        for gene, rpm in high_rpm_group:
            # counter
            now_num += 1
            if now_num % 100 == 0:
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)
            if now_num == len(high_rpm_group):
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)

            # prepare the gene df
            rpm[sp_index] = rpm.loc[:, sp_index].copy().shift(shift_codon).fillna(0)
            temp_cdt = rpm.replace(0, np.nan).groupby(['name', 'codon'])[sp_index].mean()
            ave_cdt_list.append(temp_cdt)

        # merge average density and calculate cdt
        ave_dst = pd.concat(ave_cdt_list, axis=0)
        ave_dst_title = [i + '_cdt' for i in sp_index]
        ave_dst.columns = ave_dst_title
        # high_rpf_ave_dst = pd.concat([self.high_rpm, ave_dst], axis=1)
        self.cdt = ave_dst.groupby('codon')[ave_dst_title].mean()

        for codon in self.codon_anno.keys():
            # merge the codon count
            try:
                self.codon_anno[codon].extend(self.cdt.loc[codon].values.tolist())
            except KeyError:
                self.codon_anno[codon].extend([np.nan]*self.sample_num)

        self.codon_table = pd.DataFrame(self.codon_anno).T

        self.codon_table.columns = ['AA', 'Abbr'] + ave_dst_title
        self.codon_table.sort_values('Abbr', inplace=True)
        self.codon_table.index.name = 'Codon'
        self.codon_table.to_csv(self.output + '_cdt.txt', header=True, index=True, sep='\t')

    def draw_cdt_corr(self):
        out_pdf = self.output + "_cdt_corrplot.pdf"
        out_png = self.output + "_cdt_corrplot.png"

        occ_corr = self.codon_table.filter(regex='_f0_cdt').astype('float').corr()
        occ_corr.to_csv(self.output + '_cdt_corr.txt', sep='\t', index=True)

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6 + 0.2 * self.sample_num, 7 + 0.2 * self.sample_num), dpi=300)
        sns.heatmap(data=occ_corr, annot=True, fmt=".4f", cmap="vlag", linewidths=.01,
                    cbar_kws={"orientation": "horizontal"})
        plt.tight_layout()
        plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_cdt_heat(self):
        # output the absolute occupancy heatmap
        out_pdf = self.output + "_cdt_heatplot.pdf"
        out_png = self.output + "_cdt_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6 + 0.2 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=self.codon_table.filter(regex='_f0_cdt').astype('float'),
                    annot=True, fmt=".4f", cmap="vlag", linewidths=0,
                    cbar_kws={"orientation": "horizontal", "pad": 0.05})
        plt.tight_layout()
        plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()
