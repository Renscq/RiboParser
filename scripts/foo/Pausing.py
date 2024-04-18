#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Pausing.py
# @Author  : Rensc

import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from . import RPFs


class Pausing(object):

    def __init__(self, args):
        # input and output files
        self.rpf = args.rpf
        self.output_prefix = args.output
        self.figure = args.figure
        self.all_pausing = args.all

        # opts for rpf density file parsering
        self.site = args.site
        self.frame = args.frame
        self.norm = args.normal

        self.tis = args.tis
        self.tts = args.tts
        self.rpf_num = args.min
        self.gene = args.list
        self.background = args.background

        self.sample_num = None
        self.sample_name = None

        self.raw_rpf = None
        self.merged_rpf = None

        self.gene_rpf_sum = None
        self.total_rpf_num = None
        self.high_gene = None
        self.high_rpf = None

        # opts for pausing
        self.scale = args.scale
        self.pausing = None
        self.merge_pausing = None
        self.merge_pausing_flt = None
        self.cds_pausing = None
        self.cds_codon_pausing = None
        self.codon_pausing = None

        # codon list
        self.stop = args.stop
        self.codon_dict, self.codon_table = RPFs.codon_table()
        
        if self.stop:
            del self.codon_dict['TAA']
            del self.codon_dict['TAG']
            del self.codon_dict['TGA']
            
            self.codon_table = self.codon_table[self.codon_table['Abbr'] != 'Stop']

    def read_rpf(self):
        '''
        1. import the rpf density
        2. convert the rpf density to rpm
        '''
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sites=self.site,
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis = self.tis,
                                      tts = self.tts,
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0]
        # rpf_results[0] = ''
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3]
        # rpf_results[3] = ''
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5].copy()
        self.high_gene = rpf_results[6].copy()
        # self.high_rpf = rpf_results[7].copy()
        del rpf_results

        self.merged_rpf = self.merged_rpf[self.merged_rpf['codon'].isin(self.codon_dict.keys())]

        if self.norm:
            self.merged_rpf[self.sample_name] = self.merged_rpf[self.sample_name] * 1e6 / self.total_rpf_num
            # self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
            # self.gene_rpf_sum = self.high_rpf.loc[:, ["name"] + self.sample_name].groupby("name")[self.sample_name].sum()

        self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
        del self.merged_rpf

    def get_pausing_score(self):
        '''
        for each gene and codon
        pausing score means:
        1. use mean codon rpf density as the background
        (1). filter the high expression genes
        (2). pausing score = rpf density of each codon / mean of gene density
        (3). merge all gene pausing

        or

        2. use up/down-stream codon density as the background
        (1). filter the high expression genes
        (2). pausing score = codon density / mean of upstream and downstream [number > 0] codon density
        (3). merge all gene pausing

        '''

        high_rpf_group = self.high_rpf.groupby('name')

        if self.background == 0:
            merge_pausing_list = []
            # use the mean rpf for division
            now_num = 0

            for gene in high_rpf_group:
                # counter
                now_num += 1
                if now_num % 500 == 0:
                    print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene[0]), flush=True)                    

                # prepare the gene df
                cds_idx = gene[1]['region'] == 'cds'
                gene_pausing = gene[1].iloc[:, 0:6]
                gene_rpf = gene[1].loc[:, self.sample_name].copy()
                cds_rpf = gene_rpf[cds_idx].copy()

                # calculate the mean rpf
                if self.tts != 0:
                    denominator = cds_rpf.iloc[self.tis:-self.tts].mean()
                else:
                    denominator = cds_rpf.iloc[self.tis::].mean()

                # if 0 in denominator.values:
                #     continue

                # numerator = np.asarray(gene_rpf, dtype=float)

                # calculate the PS
                pausing_score = np.where(denominator == 0, 0, np.divide(gene_rpf, denominator))

                gene_pausing[self.sample_name] = pd.DataFrame(pausing_score).values
                # gene_pausing = pd.concat([gene_mess.iloc[:, 0:6], pd.DataFrame(pausing_score)], axis=1)
                merge_pausing_list.append(gene_pausing)

            self.merge_pausing = pd.concat(merge_pausing_list, axis=0, ignore_index=True)
            self.merge_pausing.replace([np.inf, -np.inf], np.nan)

            print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene[0]), flush=True)

        else:
            temp_rpf = pd.DataFrame(np.zeros((self.background * 2, self.sample_num)).astype(int))
            temp_rpf.columns = self.sample_name.copy()
            merge_pausing_list = []
            now_num = 0
            
            # calculate the ribosome pausing score with front and behind rpf
            for gene in high_rpf_group:
                now_num += 1
                if now_num % 500 == 0:
                    print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene[0]), flush=True)
                # if now_num == len(high_rpf_group):
                #     print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene[0]), flush=True)

                gene_rpf = gene[1][self.sample_name].copy()
                gene_pausing = gene[1].iloc[:, 0:6]

                # calculate the background mean rpf
                gene_rpf_rolling = pd.concat([gene_rpf, temp_rpf], ignore_index=True).rolling(self.background * 2 + 1).mean()
                # gene_rpf_rolling = temp_rpf.append([gene_rpf, temp_rpf],ignore_index=True).rolling(self.background * 2 + 1).mean()
                gene_rpf_rolling = gene_rpf_rolling.loc[self.background * 2::]

                # numerator = np.asarray(gene_rpf, dtype=float)
                denominator = np.asarray(gene_rpf_rolling, dtype=float)
                pausing_score = np.where(denominator == 0, 0, np.divide(gene_rpf, denominator))

                gene_pausing[self.sample_name] = pd.DataFrame(pausing_score).values

                # gene_pausing = pd.concat([gene_mess.iloc[:, 0:6], pd.DataFrame(pausing_score)], axis=1)
                merge_pausing_list.append(gene_pausing)
            
            print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene[0]), flush=True)

            self.merge_pausing = pd.concat(merge_pausing_list, axis=0, ignore_index=True)
            self.merge_pausing.replace([np.inf, -np.inf], np.nan)
        del high_rpf_group

    def output_cds_pausing(self):
        # cds = self.merge_pausing[self.merge_pausing['region'] == 'cds']
        cds_group = self.merge_pausing.groupby('name')
        cds_length = cds_group[['name', 'from_tis']].tail(1).set_index('name', drop=True)
        cds_length.columns = ['aa_num']
        cds_pausing = cds_group[self.sample_name].sum()
        self.cds_pausing = pd.concat([cds_length, cds_pausing], axis=1)

        out_gene_txt = self.output_prefix + "_cds_pausing_score.txt"
        self.cds_pausing.to_csv(out_gene_txt, sep='\t', index=True)

        # utr5_pausing = self.merge_pausing[self.merge_pausing['region'] == '5utr']
        # gene_utr5_pausing = utr5_pausing.groupby('name')[self.sample_name].sum()
        # utr3_pausing = self.merge_pausing[self.merge_pausing['region'] == '3utr']
        # gene_utr3_pausing = utr3_pausing.groupby('name')[self.sample_name].sum()

    def output_cds_codon_pausing(self):
        out_cds_codon_txt = self.output_prefix + "_cds_codon_pausing_score.txt"
        titles = ['Name', 'Codon'] + [i + '_sum' for i in self.sample_name] + [i + '_mean' for i in self.sample_name]
        titles = pd.DataFrame(titles).T
        titles.to_csv(out_cds_codon_txt, sep='\t', index=False, header=False)

        # cds = self.merge_pausing[self.merge_pausing['region'] == 'cds']
        self.cds_codon_pausing = pd.pivot_table(self.merge_pausing, index=['name', 'codon'],
                                                values=self.sample_name, aggfunc=[np.sum, np.mean])
        self.cds_codon_pausing.to_csv(out_cds_codon_txt, sep='\t', mode='a', header=False)


    @staticmethod
    def scale_method(scale, pausing_score):
        if scale == 'minmax':
            # relative_pausing_score = (pausing_score - pausing_score.min()) / (pausing_score.max() - pausing_score.min())
            relative_pausing_score = pausing_score / pausing_score.max()
        elif scale == 'zscore':
            relative_pausing_score = (pausing_score - pausing_score.mean()) / pausing_score.std()
        else:
            print('Unknown scale method.', flush=True)
            sys.exit()
        return relative_pausing_score


    def output_sum_codon_pausing(self):
        '''
        @Message  : output the sum of codon pausing score
        @Input    : self.high_rpf --> high express gene rpf
        @Return   : self.codon_pausing --> codon pausing score
        @Flow     : step1 --> summary the total and valid codon number
                    step2 --> summary the pausing score
                    step3 --> calculate the relative pausing score
                    step4 --> merge pausing score of total and valid codon
        '''

        # self.codon_table = self.codon_table.reset_index(drop=False).sort_values('codon')

        # sum the valid codon num and count
        cds_rpf = self.high_rpf[self.high_rpf['region'] == 'cds']
        total_codon = cds_rpf['codon'].value_counts()

        val_rpf = cds_rpf[self.sample_name] > 0
        val_rpf.loc[:, 'codon'] = cds_rpf['codon'].copy()
        valid_codon = val_rpf.groupby('codon')[self.sample_name].sum()
        rpf_count = cds_rpf.groupby('codon')[self.sample_name].sum()

        # rename or format
        total_codon.columns = 'total_codon'
        valid_codon.columns = [i + '_valid_codon' for i in self.sample_name]
        rpf_count.columns = [i + '_rpf_count' for i in self.sample_name]

        # summary the total codon pausing score
        total_codon_pausing = self.merge_pausing.groupby('codon')[self.sample_name].mean()
        total_codon_pausing.columns = [i + '_absolute_total_ps' for i in self.sample_name]

        valid_pausing = self.merge_pausing[self.merge_pausing[self.sample_name].max(axis=1) > 0]
        valid_codon_pausing = valid_pausing.groupby('codon')[self.sample_name].mean()
        valid_codon_pausing.columns = [i + '_absolute_valid_ps' for i in self.sample_name]

        # calculate the relative pausing score
        relative_total_codon_pausing = self.scale_method(self.scale, total_codon_pausing)
        relative_valid_codon_pausing = self.scale_method(self.scale, valid_codon_pausing)

        relative_total_codon_pausing.columns = [i + '_relative_total_ps' for i in self.sample_name]
        relative_valid_codon_pausing.columns = [i + '_relative_valid_ps' for i in self.sample_name]

        # ouput the total codon pausing score
        out_codon_txt = self.output_prefix + "_sum_codon_pausing_score.txt"
        self.codon_pausing = pd.concat([self.codon_table, total_codon, valid_codon, rpf_count, total_codon_pausing,
                                       valid_codon_pausing, relative_total_codon_pausing, relative_valid_codon_pausing], 
                                       axis=1, join='inner')
        
        self.codon_pausing.index.name = 'Codon'
        self.codon_pausing = self.codon_pausing.sort_values('Abbr').reset_index(drop=False)
        self.codon_pausing.sort_values(by=['Abbr', 'Codon'], ascending=True, inplace=True)
        self.codon_pausing.to_csv(out_codon_txt, sep='\t', index=False)

    def output_all_pausing(self):
        if self.all_pausing:
            all_pausing_txt = self.output_prefix + "_all_pausing_score.txt"
            self.merge_pausing.to_csv(all_pausing_txt, sep='\t', index=False)

    def draw_codon_total_pausing_heat(self):
        '''
        @Message  : draw the codon total pausing score heatmap
        @Input    : self.codon_pausing --> codon pausing score
                    self.output_prefix --> prefix of output filename
        @Return   : output --> codon total pausing score heatmap
        @Flow     : step1 --> select the column of total pausing score
                    step2 --> draw the heatmap
        '''
        
        out_pdf = "{prefix}_total_pausing_heatplot.pdf".format(prefix=self.output_prefix)
        out_png = "{prefix}_total_pausing_heatplot.png".format(prefix=self.output_prefix)

        ps_title = [i + '_absolute_total_ps' for i in self.sample_name]

        # select the column of total pausing score
        total_codon_pausing = self.codon_pausing[['Codon', 'Abbr'] + ps_title].copy()

        # draw the rpf pausing score
        matplotlib.use('AGG')
        fig = plt.figure(figsize=(10 + 0.5 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=total_codon_pausing[ps_title].astype('float16'),
                    yticklabels=total_codon_pausing['Codon'] + '[' + total_codon_pausing['Abbr'] + ']',
                    annot=True, fmt=".3f", cmap="vlag", linewidths=.01)

        # plt.suptitle("codon pausing")
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_codon_valid_pausing_heat(self):
        '''
        @Message  : draw the codon total pausing score heatmap
        @Input    : self.codon_pausing --> codon pausing score
                    self.output_prefix --> prefix of output filename
        @Return   : output --> codon valid pausing score heatmap
        @Flow     : step1 --> select the column of valid pausing score
                    step2 --> draw the heatmap
        '''

        out_pdf = "{prefix}_valid_pausing_heatplot.pdf".format(prefix=self.output_prefix)
        out_png = "{prefix}_valid_pausing_heatplot.png".format(prefix=self.output_prefix)

        ps_title = [i + '_absolute_valid_ps' for i in self.sample_name]

        # select the column of valid pausing score
        valid_codon_pausing = self.codon_pausing[['Codon', 'Abbr'] + ps_title].copy()

        # draw the rpf pausing score
        matplotlib.use('AGG')
        fig = plt.figure(figsize=(10 + 0.5 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=valid_codon_pausing[ps_title].astype('float16'),
                    yticklabels=valid_codon_pausing['Codon'] + '[' + valid_codon_pausing['Abbr'] + ']',
                    annot=True, fmt=".3f", cmap="vlag", linewidths=.01)

        # plt.suptitle("codon pausing")
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_codon_pausing_plot(self, category):
        out_pdf = "{prefix}_{category}_pausing_plot.pdf".format(prefix=self.output_prefix, category=category)
        out_png = "{prefix}_{category}_pausing_plot.png".format(prefix=self.output_prefix, category=category)

        # draw the rpf pausing score
        matplotlib.use('AGG')
        fig = plt.figure(figsize=(9, 3 * self.sample_num), dpi=300)

        for fig_num in range(self.sample_num):
            ax = fig.add_subplot(self.sample_num, 1, fig_num + 1)
            sns.barplot(x=self.codon_pausing.index,
                        y=self.codon_pausing[self.sample_name[fig_num] + '_' + category + '_ps'],
                        ax=ax)
            ax.set_xticks(self.codon_pausing.index)
            ax.set_xticklabels(self.codon_pausing['Codon'] + " [" + self.codon_pausing['Abbr'] + "]",
                               size=8,
                               va='top')
            ax.set_xlim([-1, 65])

            # add the gene annotation here, need to import the txt file
            ax.set_ylabel('Pausing score')
            ax.set_xlabel('Codon')
            ax.set_title(self.sample_name[fig_num])
            plt.xticks(rotation=90)

        plt.suptitle("{category} codon pausing".format(category=category))
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    @staticmethod
    def get_labels(now_gene):
        start_site = now_gene["from_tis"].iloc[0]
        start_codon_site = 0
        stop_codon_site = int(now_gene[now_gene["from_tts"] == 0]["from_tis"])
        stop_site = now_gene["from_tis"].iloc[-1]

        now_labels = [start_site, start_codon_site, stop_codon_site, stop_site]

        return now_labels

    def draw_rpf_pausing_plot(self):

        merge_pausing_group = self.merge_pausing.groupby('name')

        for gene in merge_pausing_group:
            print("Now gene: {gene}.".format(gene=gene[0]), flush=True)

            now_gene_pausing = gene[1]
            now_labels = self.get_labels(now_gene_pausing)

            # draw the rpf pausing score
            matplotlib.use('AGG')
            fig = plt.figure(figsize=(9, 3 * self.sample_num), dpi=300)

            for fig_num in range(self.sample_num):
                ax1 = fig.add_subplot(self.sample_num, 1, fig_num + 1)
                ax1.bar(now_gene_pausing.from_tis,
                        now_gene_pausing[self.sample_name[fig_num]],
                        width=2,
                        color="#f86934")

                ax1.set_xticks(now_labels)
                ax1.set_xticklabels(now_labels)
                ax1.set_xlim([now_labels[0], now_labels[-1]])

                # add the gene annotation here, need to import the txt file
                ax1.set_ylabel('Pausing score')
                ax1.set_xlabel('Position (nt)')
                ax1.set_title(self.sample_name[fig_num])
                plt.xticks(rotation=90)

            plt.suptitle(str(gene[0]))
            plt.tight_layout()
            # plt.show()

            if self.figure == 'pdf':
                out_pdf = str(gene[0]) + "_pausing_plot.pdf"
                fig.savefig(fname=out_pdf)
            elif self.figure == 'png':
                out_png = str(gene[0]) + "_pausing_plot.png"
                fig.savefig(fname=out_png)
            else:
                print('Unknown figure format.', flush=True)
                sys.exit()
            plt.close()
