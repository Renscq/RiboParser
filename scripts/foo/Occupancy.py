
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Occupancy.py


import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from . import RPFs


class Occupancy(object):

    def __init__(self, args):
        # input and output
        self.gene = args.list
        self.rpf = args.rpf
        self.output = args.output
        self.output_all = args.all

        # filter the stable codon coverage range
        self.norm = args.normal
        self.site = args.site
        self.frame = args.frame
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts

        # titles
        self.sample_num = 0
        self.sample_name = None
        self.dst_title = None
        self.occ_title = None
        self.rel_title = None

        # for the high expression gene filter
        self.merged_rpf = None
        self.raw_rpf = None
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None

        # table for codon occupancy
        self.scale = args.scale
        self.total_density = None
        self.occupancy = None

        self.stop = args.stop
        self.codon_dict, self.codon_table = RPFs.codon_table()
        
        if self.stop:
            del self.codon_dict['TAA']
            del self.codon_dict['TAG']
            del self.codon_dict['TGA']
            
            self.codon_table = self.codon_table[self.codon_table['Abbr'] != 'Stop']

    def import_rpf(self):
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sites=self.site,
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis=self.tis,
                                      tts=self.tts,
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0]
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
            # self.gene_rpf_sum = self.high_rpf.loc[:, ["name"] + self.sample_name].groupby("name")[self.sample_name].sum()

        self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
        del self.merged_rpf

        # self.high_rpf.loc[:, self.sample_name] = self.high_rpf[self.sample_name] * 1e6 / self.total_rpf_num

    @staticmethod
    def scale_method(scale, occupancy):
        if scale == 'minmax':
            # relative_occupancy = (occupancy - occupancy.min()) / (occupancy.max() - occupancy.min())
            relative_occupancy = occupancy / occupancy.max()
        elif scale == 'zscore':
            relative_occupancy = (occupancy - occupancy.mean()) / occupancy.std()
        else:
            print('Unknown scale method.', flush=True)
            sys.exit()
            
        return relative_occupancy
    
    def codon_occupancy(self):
        # calc the average density of each gene
        # ave_rpf = self.high_rpf.replace(0, np.nan).groupby('name')[self.sample_name].mean()
        ave_rpf = self.high_rpf.groupby('name')[self.sample_name].mean()

        # calculate the average density of each codon
        high_rpf_group = self.high_rpf.groupby('name')
        now_num = 0
        ave_dst_list = []
        for gene, rpfs in high_rpf_group:
            # counter
            now_num += 1
            if now_num % 500 == 0:
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)

            # prepare the gene df
            temp_rpfs = rpfs.loc[:, self.sample_name].copy()
            # temp_ave_dst = np.divide(temp_rpfs, ave_rpf.loc[gene].to_list(), where=ave_rpf.loc[gene] != 0)
            # where is discarded in the new version of numpy
            temp_ave_dst = np.divide(temp_rpfs, ave_rpf.loc[gene].to_list())
            ave_dst_list.append(temp_ave_dst)

        print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)
        
        # merge average density
        ave_dst = pd.concat(ave_dst_list, axis=0)
        ave_dst_title = [i + '_ave_dst' for i in self.sample_name]
        ave_dst.columns = ave_dst_title
        high_rpf_ave_dst = pd.concat([self.high_rpf, ave_dst], axis=1)
        codon_density = pd.pivot_table(high_rpf_ave_dst,
                                       index=['name', 'codon'],
                                       values=self.sample_name + ave_dst_title,
                                       aggfunc=np.sum)
        if self.output_all:
            high_rpf_ave_dst.to_csv(self.output + '_rpf_density.txt', sep='\t', index=False)

        codon_density.to_csv(self.output + '_codon_density.txt', sep='\t', index=True)

        # calculate the total density and occupancy (average density)
        self.total_density = high_rpf_ave_dst.replace(0, np.nan).groupby('codon')[ave_dst_title].sum()
        self.occupancy = high_rpf_ave_dst.replace(0, np.nan).groupby('codon')[ave_dst_title].mean()

        # merge the codon count
        for codon in self.codon_dict.keys():
            try:
                self.codon_dict[codon].extend(self.total_density.loc[codon].values.tolist())
                self.codon_dict[codon].extend(self.occupancy.loc[codon].values.tolist())
            except KeyError:
                self.codon_dict[codon].extend([np.nan] * self.sample_num)
                self.codon_dict[codon].extend([np.nan] * self.sample_num)

        self.codon_table = pd.DataFrame(self.codon_dict).T
        self.dst_title = [i + '_density' for i in self.sample_name]
        self.occ_title = [i + '_absolute_occupancy' for i in self.sample_name]
        self.codon_table.columns = ['AA', 'Abbr'] + self.dst_title + self.occ_title

        # z-score standard
        # rel_occ = (self.codon_table[self.occ_title] - self.codon_table[self.occ_title].mean()) / self.codon_table[self.occ_title].std()
        rel_occ = self.scale_method(self.scale, self.codon_table[self.occ_title])
        self.rel_title = [i + '_relative_occupancy' for i in self.sample_name]
        rel_occ.columns = self.rel_title
        self.codon_table = pd.concat([self.codon_table, rel_occ], axis=1)
        self.codon_table.sort_values('Abbr', inplace=True)

        # self.codon_table = self.codon_table.round(2)
        self.codon_table.index.name = 'Codon'
        self.codon_table = self.codon_table.reset_index(drop=False)
        self.codon_table.sort_values(by=['Abbr', 'Codon'], inplace=True)
        self.codon_table.to_csv(self.output + '_codon_occupancy.txt', header=True, index=False, sep='\t')

    def draw_occupancy_corr(self):
        out_pdf = self.output + "_occupancy_corrplot.pdf"
        out_png = self.output + "_occupancy_corrplot.png"

        occ_corr = self.codon_table[self.occ_title].astype('float').corr()
        occ_corr.to_csv(self.output + '_occupancy_corr.txt', sep='\t', index=True)

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(7 + 0.2 * self.sample_num, 6 + 0.2 * self.sample_num), dpi=300)
        sns.heatmap(data=occ_corr, annot=True, fmt=".4f", cmap="vlag", linewidths=.01)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_occupancy_heat(self):
        # output the absolute occupancy heatmap
        out_pdf = self.output + "_occupancy_heatplot.pdf"
        out_png = self.output + "_occupancy_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(10 + 0.5 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=self.codon_table[self.occ_title].astype('float'),
                    yticklabels=self.codon_table['Codon'] + '[' + self.codon_table['Abbr'] + ']',
                    annot=True,
                    fmt=".3f",
                    cmap="vlag",
                    linewidths=.01)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_occupancy_relative_heat(self):

        # output the relative occupancy heatmap
        out_pdf = self.output + "_occupancy_relative_heatplot.pdf"
        out_png = self.output + "_occupancy_relative_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(10 + 0.5 * self.sample_num, 20), dpi=300)
        sns.heatmap(data=self.codon_table[self.rel_title].astype('float'),
                    yticklabels=self.codon_table['Codon'] + '[' + self.codon_table['Abbr'] + ']',
                    annot=True,
                    fmt=".3f",
                    cmap="vlag",
                    linewidths=.01)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_occupancy_line(self):

        codon_occp = pd.melt(self.codon_table,
                             id_vars=['Codon', 'AA', 'Abbr'],
                             var_name='Groups',
                             value_name='Occupancy')
        codon_occp.loc[:, 'Codon'] = codon_occp['Codon'] + '[' + codon_occp['Abbr'] + ']'
        codon_occp = codon_occp.loc[codon_occp['Groups'].str.contains('relative'), ]

        out_pdf = self.output + "_occupancy_relative_lineplot.pdf"
        out_png = self.output + "_occupancy_relative_lineplot.png"

        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(figsize=(11, 5.5), dpi=300)

        flag = 1
        for aa in codon_occp['AA'].unique().tolist():
            tmp = codon_occp.loc[codon_occp['AA'] == aa, :]
            if flag == 1:
                sns.lineplot(data=tmp,
                             x='Codon',
                             y='Occupancy',
                             linewidth=1.5,
                             hue='Groups',
                             style='Groups',
                             markers=True,
                             legend=True,
                             ax=ax)
                flag += 1
            else:
                sns.lineplot(data=tmp,
                             x='Codon',
                             y='Occupancy',
                             linewidth=1.5,
                             hue='Groups',
                             style='Groups',
                             markers=True,
                             legend=False,
                             ax=ax)

        # add the gene annotation here, need to import the txt file
        ax.legend(loc=8, ncol=4, bbox_to_anchor=(0.5, -0.5), fontsize=6)

        # ax.set_ylabel('occupancy')
        # ax.set_xlabel('codon')

        plt.xticks(rotation=90)

        plt.title('Relative codon occupancy')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
