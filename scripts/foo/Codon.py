#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Codon.py


import CAI
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from . import RPFs


class Codon(object):

    def __init__(self, args):
        # input and output
        self.gene = args.list
        self.rpf = args.rpf
        self.output = args.output

        # for the high expression gene filter
        self.sample_num = 0
        self.sample_name = None
        self.rel_name = None
        self.merged_rpf = None

        self.rpf_num = args.min
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None
        self.total_rpf_num = None

        # filter the stable codon coverage range
        self.tis = args.tis
        self.tts = args.tts
        self.site = args.site
        self.frame = args.frame
        # self.mrna_seq = args.sequence
        # self.mrna_dict = OrderedDict()
        # self.longest = str(args.longest)
        self.codon_dict, self.codon_table = RPFs.codon_table()
        self.codon_title = ['AA', 'Abbr.']

    def read_rpf(self):
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf, sites=self.site, frame=self.frame,
                                      # tis=self.tis, tts=self.tts,
                                      sample_num=None, sample_name=None,
                                      gene=self.gene, rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0]
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        self.gene_rpf_sum = rpf_results[5]
        # self.high_gene = rpf_results[6]
        self.high_rpf = rpf_results[7]
        self.high_rpf[self.sample_name] = self.high_rpf[self.sample_name] * 10000000 / self.total_rpf_num

    def merge_codon_table(self, count, freq, rscu):
        # synonymous codon


        # merge all codon usage
        for codon in self.codon_dict.keys():
            # merge the codon count
            try:
                self.codon_dict[codon].append(round(count[codon], 3))
            except KeyError:
                self.codon_dict[codon].append(np.nan)
            # merge the frequency of the codon counts of per 1000
            try:
                self.codon_dict[codon].append(round(freq[codon], 3))
            except KeyError:
                self.codon_dict[codon].append(np.nan)
            # merge the RSCU
            try:
                self.codon_dict[codon].append(round(rscu[codon], 3))
            except KeyError:
                self.codon_dict[codon].append(np.nan)

    def get_all_codon_usage(self):
        print('Whole gene sequence.', flush=True)
        all_cds = self.merged_rpf.loc[self.merged_rpf['region'] == 'cds', 'codon']
        all_cds_codon_count = all_cds.value_counts()
        all_codon_freq = round(all_cds_codon_count * 1000 / all_cds_codon_count.sum(), 3)
        all_cds_codon_count_dict = all_cds_codon_count.to_dict()
        all_codon_freq_dict = all_codon_freq.to_dict()
        all_rscu_dict = CAI.RSCU(all_cds.to_list())
        self.merge_codon_table(all_cds_codon_count_dict, all_codon_freq_dict, all_rscu_dict)
        self.codon_title.extend(['All_Count', 'All_Freq', 'All_RSCU'])

    def get_high_codon_usage(self, high_rpf_cds, sample):
        cds = high_rpf_cds
        cds_codon_count = cds.value_counts()
        codon_freq = round(cds_codon_count * 1000 / cds_codon_count.sum(), 3)
        cds_codon_count_dict = cds_codon_count.to_dict()
        codon_freq_dict = codon_freq.to_dict()
        rscu_dict = CAI.RSCU(cds.to_list())
        self.merge_codon_table(cds_codon_count_dict, codon_freq_dict, rscu_dict)
        self.codon_title.extend([sample + '_Count', sample + '_Freq', sample + '_RSCU'])

    def get_codon_usage(self):

        for sample in self.sample_name:
            print('Now sample is: {sample}.'.format(sample=sample), flush=True)
            # filter the high expression genes
            high_gene = self.gene_rpf_sum[self.gene_rpf_sum[sample] > self.rpf_num]
            gene_ids = list(high_gene.index)
            high_rpf = self.merged_rpf[self.merged_rpf['name'].isin(gene_ids)].copy().reset_index(drop=True)
            high_rpf_cds = high_rpf.loc[high_rpf['region'] == 'cds', 'codon']
            self.get_high_codon_usage(high_rpf_cds, sample)

        self.codon_table = pd.DataFrame(self.codon_dict).T
        self.codon_table.columns = self.codon_title
        self.codon_table.sort_values('Abbr.', inplace=True)
        self.codon_table.index.name = 'codon'
        self.codon_table = self.codon_table.reset_index(drop=False)

        self.codon_table.to_csv(self.output + '_codon_usage.txt', header=True, index=False, sep='\t')

    def draw_whole_gene(self, codon_table):
        counts = codon_table['All_Count']
        freq = codon_table['All_Freq']
        rscu = codon_table['All_RSCU']

        # draw the figure
        matplotlib.use('AGG')

        # output figure name
        out_pdf = self.output + "_codon_usage.pdf"
        out_png = self.output + "_codon_usage.png"
        fig = plt.figure(figsize=(12, 6), dpi=300)

        ax1 = fig.add_subplot(311)
        sns.barplot(x=np.arange(0, 61), y=counts, ax=ax1)
        ax1.set_title('codon counts')
        ax1.set_xticks(np.arange(0, 61))
        ax1.set_xticklabels(codon_table['xticks'], rotation=90)

        ax2 = fig.add_subplot(312)
        sns.barplot(x=np.arange(0, 61), y=freq, ax=ax2)
        ax2.set_title('codon frequency')
        ax2.set_xticks(np.arange(0, 61))
        ax2.set_xticklabels(codon_table['xticks'], rotation=90)

        ax3 = fig.add_subplot(313)
        sns.barplot(x=np.arange(0, 61), y=rscu, ax=ax3)
        ax3.set_title('RSCU')
        ax3.set_xticks(np.arange(0, 61))
        ax3.set_xticklabels(codon_table['xticks'], rotation=90)

        plt.tight_layout()
        plt.show()
        plt.close()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_codon_usage(self):
        codon_table = self.codon_table.copy()
        codon_table['codon'] = codon_table.index
        codon_table.sort_values(by=['Abbr.', 'codon'], ascending=[True, True], inplace=True)
        codon_table['xticks'] = codon_table.index + '[' + codon_table['Abbr.'] + ']'

        self.draw_whole_gene(codon_table)

        for num in range(self.sample_num):
            counts = codon_table[self.sample_name[num] + '_Count']
            freq = codon_table[self.sample_name[num] + '_Freq']
            rscu = codon_table[self.sample_name[num] + '_RSCU']

            # draw the figure
            matplotlib.use('AGG')

            # output figure name
            out_pdf = self.sample_name[num] + "_codon_usage.pdf"
            out_png = self.sample_name[num] + "_codon_usage.png"
            fig = plt.figure(figsize=(12, 6), dpi=300)

            ax1 = fig.add_subplot(311)
            sns.barplot(x=np.arange(0, 61), y=counts, ax=ax1)
            ax1.set_title('codon counts')
            ax1.set_xticks(np.arange(0, 61))
            ax1.set_xticklabels(codon_table['xticks'], rotation=90)

            ax2 = fig.add_subplot(312)
            sns.barplot(x=np.arange(0, 61), y=freq, ax=ax2)
            ax2.set_title('codon frequency')
            ax2.set_xticks(np.arange(0, 61))
            ax2.set_xticklabels(codon_table['xticks'], rotation=90)

            ax3 = fig.add_subplot(313)
            sns.barplot(x=np.arange(0, 61), y=rscu, ax=ax3)
            ax3.set_title('RSCU')
            ax3.set_xticks(np.arange(0, 61))
            ax3.set_xticklabels(codon_table['xticks'], rotation=90)

            plt.tight_layout()
            plt.show()
            plt.close()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
