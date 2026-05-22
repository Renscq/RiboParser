#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Percentage.py


from concurrent.futures import ThreadPoolExecutor
from scipy.stats import zscore
from math import ceil
from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from . import RPFs


class Percentage(object):

    def __init__(self, args):
        # input and output
        self.rpf = args.rpf
        self.gene = args.transcript
        self.gene_table = None
        self.output = args.output

        # filter the stable codon coverage range
        self.tis = None
        self.tts = None
        self.site = 'P'
        self.frame = args.frame
        self.rpf_num = args.min

        self.norm = args.normal
        self.sample_name = None

        # for the high expression gene filter
        self.raw_rpf = None
        # self.merged_rpf = None
        # self.high_gene = None
        self.high_rpf = None
        self.gene_num = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None

        # gene coverage
        self.gene_coverage = None


    def read_rpf(self):
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
        self.sample_name = rpf_results[1].copy()
        self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        self.gene_num = len(self.high_gene)
        self.high_rpf = rpf_results[7].copy()
        del rpf_results

        if self.norm:
            self.high_rpf[self.sample_name] = self.high_rpf[self.sample_name].astype(float)

            self.high_rpf[self.sample_name] = self.high_rpf[self.sample_name].div(self.total_rpf_num) * 1e6


    def import_gene(self):
        '''
        @Message  : import the gene file
        @Input    : self.gene --> gene file
        @Return   : self.gene --> gene dataframe
        @Flow     : step1 --> import the gene file
        '''

        self.gene_table = pd.read_csv(self.gene, sep='\t', index_col=False)
        self.gene_table = self.gene_table.loc[self.gene_table['transcript_id'].isin(self.high_gene), ]


    def calc_density_percent(self):
        '''
        @Message  : calculate the gene rpf density percentage of each gene
        @Input    : gene_rpf --> gene rpf dataframe
        @Return   : gene_rpf --> gene rpf density dataframe
        @Flow     : step1 --> calculate the gene rpf density
        '''
        # remove the UTR5 and UTR3 region
        self.high_rpf = self.high_rpf.loc[self.high_rpf['region'] == 'cds', ]
        
        # group the gene by the gene name and sum the value
        self.gene_rpf_sum = self.high_rpf.groupby('name', observed=False)[self.sample_name].sum()

        # convert the rpf value to 1, if the value > 0
        self.high_rpf[self.sample_name] = self.high_rpf[self.sample_name].map(lambda x: 1 if x > 0 else 0)

        # group by the gene and sum the value, then divide the gene length
        gene_p_site_num = self.high_rpf.groupby('name', observed=False)[self.sample_name].sum()
        gene_length = self.high_rpf.groupby('name', observed=False)['from_tis'].max()
        self.gene_coverage = gene_p_site_num.div(gene_length, axis=0) * 100


    def draw_rpf_histogram(self):
        '''
        @Message  : draw the rpf histogram of the gene
        @Input    : self.gene_coverage --> the gene coverage dataframe
        @Return   : None
        @Flow     : step1 --> draw the rpf histogram of the gene
        '''

        # draw the gene coverage histogram of the sample
        for sp in self.sample_name:
            print('Draw percentage histogram of {sp}.'.format(sp=sp), flush=True)

            out_pdf = self.output + "_" + sp + "_coverage_histogram.pdf"
            out_png = self.output + "_" + sp + "_coverage_histogram.png"

            matplotlib.use('AGG')

            fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), gridspec_kw={'width_ratios': [1, 1]})

            # draw the gene coverage histogram of the sample
            sns.histplot(self.gene_coverage[sp], bins=100, shrink=0.9, 
                         kde=True, edgecolor='white', color='#2b79da', ax=axes[0])
            axes[0].set_ylabel('Gene count', fontsize=12)
            axes[0].set_xlabel("Percentage (%)", fontsize=12)
            axes[0].set_title(('Coverage histogram ({number})').format(number=self.gene_num), fontsize=14)

            sns.histplot(np.log2(self.gene_rpf_sum[sp] + 1), bins=100, shrink=0.9, 
                         kde=True, edgecolor='white', color='#2b79da', ax=axes[1])
            axes[1].set_ylabel('Gene count', fontsize=12)
            axes[1].set_xlabel("log2(Expr + 1)", fontsize=12)
            axes[1].set_title(('Abundance histogram ({number})').format(number=self.gene_num), fontsize=14)

            fig.tight_layout()
            # plt.show()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()


    def draw_rpf_boxplot(self):
        '''
        @Message  : draw the boxplot of the percentage
        @Input    : self.gene_coverage --> the gene coverage dataframe
        @Return   : None
        @Flow     : step1 --> draw the rpf boxplot of the gene
        '''

        gene_coverage = pd.melt(self.gene_coverage.reset_index(), id_vars=['name'], var_name='Sample', value_name='Percentage')

        gene_rpf_sum_log2 = np.log2(self.gene_rpf_sum + 1)
        gene_expression = pd.melt(gene_rpf_sum_log2.reset_index(), id_vars=['name'], var_name='Sample', value_name='Expression')

        # draw the gene coverage boxplot of the sample
        print('Draw percentage boxplot of RPFs coverage.', flush=True)

        out_pdf = self.output + "_coverage_boxplot.pdf"
        out_png = self.output + "_coverage_boxplot.png"

        width = 0 if self.sample_num < 6 else (self.sample_num - 6) * 0.3

        matplotlib.use('AGG')

        fig, axes = plt.subplots(1, 2, figsize=(7 + width, 5), gridspec_kw={'width_ratios': [1, 1]})

        # draw the gene coverage boxplot of the sample
        sns.boxplot(data = gene_coverage, x = 'Sample', y='Percentage', 
                    fliersize = 1, linewidth = 0.5,
                    color='#2b79da', ax=axes[0])
        axes[0].set_ylabel('Percentage (%)', fontsize=12)
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=90, ha='right', va='top', fontsize=12, color='black')
        axes[0].set_xlabel("Sample", fontsize=12)
        axes[0].set_title(('Coverage boxplot ({number})').format(number=self.gene_num), fontsize=14)

        sns.boxplot(data = gene_expression, x = 'Sample', y='Expression', 
                    fliersize = 1, linewidth = 0.5,
                    color='#2b79da', ax=axes[1])
        axes[1].set_ylabel('Expression (log2)', fontsize=12)
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=90, ha='right', va='top', fontsize=12, color='black')
        axes[1].set_xlabel("Sample", fontsize=12)
        axes[1].set_title(('Abundance boxplot ({number})').format(number=self.gene_num), fontsize=14)

        fig.tight_layout()
        # plt.show()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()


    def output_density_percent(self):
        '''
        @Message  : output the gene rpf density percentage
        @Input    : self.gene_coverage --> the gene coverage dataframe
        @Return   : None
        @Flow     : step1 --> output the gene rpf density percentage
        '''

        # output the gene coverage of the sample
        out_file = self.output + "_gene_coverage_percent.txt"
        self.gene_coverage.to_csv(out_file, sep='\t', index=True, header=True)
        