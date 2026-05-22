#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Density.py


import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


class Density(object):

    def __init__(self, args):
        self.gene = args.list
        self.rpf = args.rpf
        self.sample_num = 0
        self.sample_name = None
        self.rel_name = None
        self.merged_rpf = None
        self.merged_rpf_na = None

        # for the high expression gene filter
        self.rpf_num = args.min
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None
        # filter the stable codon coverage range
        self.site = args.site
        self.frame = args.frame
        self.tis = args.tis
        self.tts = args.tts

        # table for codon occupancy
        self.occupancy = None
        self.occ_title = None
        self.rel_title = None
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
        rpf = pd.read_csv(self.rpf, sep='\t', header=0, names=None)
        self.sample_num = int((len(rpf.columns) - 6) / 3)
        self.sample_name = pd.Series(rpf.columns[6::]).str[:-3].drop_duplicates().tolist()
        # merge the RPFs of nucleotide to codon
        # front 6 columns ['name', 'now_nt', 'from_tis', 'from_tes', 'region', 'codon']
        self.merged_rpf = rpf.iloc[:, 0:6].copy()
        self.get_frame_rpf(rpf)

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

    def codon_occupancy(self):
        # filter the high expression gene
        self.total_rpf_num = self.merged_rpf[self.sample_name].sum()
        self.gene_rpf_sum = self.merged_rpf.groupby('name')[self.sample_name].sum()
        self.high_gene = self.gene_rpf_sum.loc[self.gene_rpf_sum.max(axis=1) >= self.rpf_num, :].index
        self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene), :]

        # calc the average density of each gene
        ave_rpf = self.high_rpf.replace(0, np.nan).groupby('name')[self.sample_name].mean()
        ave_dst_title = [i + '_ave_dst' for i in self.sample_name]
        np.seterr(divide='ignore', invalid='ignore')
        # ave_dst = self.high_rpf.groupby('name').apply(lambda x: np.divide(x[self.sample_name], ave_rpf.loc[x.name]))
        # ave_dst.columns = ave_dst_title
        # high_rpf_ave_dst = pd.concat([self.high_rpf, ave_dst], axis=1)
        # self.occupancy = high_rpf_ave_dst.replace(0, np.nan).groupby('codon')[ave_dst_title].mean()

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
        high_rpf_group = self.high_rpf.groupby('name')
        now_num = 0
        ave_dst_list = []
        for gene, rpfs in high_rpf_group:
            # counter
            now_num += 1
            if now_num % 100 == 0:
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)
            if now_num == len(high_rpf_group):
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)

            # prepare the gene df
            temp_rpfs = rpfs.loc[:, self.sample_name].copy().shift(shift_codon).fillna(0)
            temp_ave_dst = np.divide(temp_rpfs, ave_rpf.loc[gene])
            ave_dst_list.append(temp_ave_dst)

        # merge average density and calculate occupancy
        ave_dst = pd.concat(ave_dst_list, axis=0)
        ave_dst.columns = ave_dst_title
        high_rpf_ave_dst = pd.concat([self.high_rpf, ave_dst], axis=1)
        self.occupancy = high_rpf_ave_dst.replace(0, np.nan).groupby('codon')[ave_dst_title].mean()

        for codon in self.codon_anno.keys():
            # merge the codon count
            try:
                self.codon_anno[codon].extend(self.occupancy.loc[codon].values.tolist())
            except KeyError:
                self.codon_anno[codon].extend([np.nan]*self.sample_num)

        self.codon_table = pd.DataFrame(self.codon_anno).T

        # z-score standard
        rel_occ = (self.codon_table.iloc[:, 2::] -
                   self.codon_table.iloc[:, 2::].mean())/self.codon_table.iloc[:, 2::].std()
        self.codon_table = pd.concat([self.codon_table, rel_occ], axis=1)
        self.occ_title = [i + '_occupancy' for i in self.sample_name]
        self.rel_title = [i + '_relative' for i in self.sample_name]
        self.codon_table.columns = ['AA', 'Abbr.'] + self.occ_title + self.rel_title
        self.codon_table.sort_values('Abbr.', inplace=True)
        # self.codon_table = self.codon_table.round(2)
        self.codon_table.index.name = 'codon'
        self.codon_table.to_csv(self.output + '_codon_occupancy.txt', header=True, index=True, sep='\t')

    def draw_occupancy_corr(self):
        out_pdf = self.output + "_occupancy_corrplot.pdf"
        out_png = self.output + "_occupancy_corrplot.png"

        occ_corr = self.codon_table[self.occ_title].astype('float').corr()
        occ_corr.to_csv(self.output + '_occupancy_corr.txt', sep='\t', index=True)

        matplotlib.use('AGG')
        sns.heatmap(data=occ_corr, annot=True, fmt=".4f", cmap="vlag", linewidths=.01,
                    cbar_kws={"orientation": "horizontal"})
        plt.tight_layout()
        plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_occupancy_heat(self):
        # output the absolute occupancy heatmap
        out_pdf = self.output + "_occupancy_heatplot.pdf"
        out_png = self.output + "_occupancy_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6 + 0.2 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=self.codon_table[self.occ_title].astype('float'),
                    annot=True, fmt=".4f", cmap="vlag", linewidths=.01,
                    cbar_kws={"orientation": "horizontal", "pad": 0.05})
        plt.tight_layout()
        plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_occupancy_relative_heat(self):

        # output the relative occupancy heatmap
        out_pdf = self.output + "_occupancy_relative_heatplot.pdf"
        out_png = self.output + "_occupancy_relative_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6 + 0.2 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=self.codon_table[self.rel_title].astype('float'),
                    annot=True, fmt=".4f", cmap="vlag", linewidths=.01,
                    cbar_kws={"orientation": "horizontal", "pad": 0.05})
        plt.tight_layout()
        plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()
