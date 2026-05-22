#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : Shift.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2025/02/18 14:58:41
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


from . import RPFs
from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import polars as pl
import pandas as pd
import numpy as np
import seaborn as sns
import math


class Shift(object):

    def __init__(self, args):
        # input and output
        self.transcript = args.transcript
        self.rpf = args.rpf
        self.output = args.output

        # filter the range around tis/tts
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts

        # for the high expression gene filter
        self.sample_num = 0
        self.sample_name = None
        self.raw_rpf = None
        self.high_gene = None

        self.period = args.period
        self.period_count = None
                
        self.gene_frame_shift = pd.DataFrame(columns=['Sample', 'Frame', 'Count'])

    def import_rpf(self):
        '''
        @Message  : function for rpf import.
        @Input    : self.ribo --> rpf density file
                    self.tis --> tis range
                    self.tts --> tts range
                    self.rpf_num --> rpf number filter
        @Return   : output --> rpf density dataframe
        '''

        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sites='P',
                                      frame='all',
                                      sample_num=None,
                                      sample_name=None,
                                      tis = self.tis,
                                      tts = self.tts,
                                      gene=self.transcript,
                                      rpf_num=self.rpf_num)
        
        self.raw_rpf = rpf_results[0].to_pandas()
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        # self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        # self.high_rpf = rpf_results[7]

        del rpf_results

    def calc_3nt_period(self):
        '''
        @Message  : function for each gene 3nt periodicity calculation.
        @Input    : param --> rpf density file
        @Return   : output --> summary the 3nt periodicity of each gene
        @Flow     : step1 --> calculate the count of each frame of samples
                    step2 --> calculate the ratio of each frame of samples
        '''

        raw_rpf = self.raw_rpf.drop(columns=["now_nt", "from_tis", "from_tts", "region", "codon"])

        # calculate the count of each frame of samples
        self.period_count  = raw_rpf.groupby('name').apply(np.sum, axis=0).drop(columns=["name"])
        # self.period_count = self.period_count .reset_index()

    def output_meta(self):
        out_txt = self.output + "_gene_periodicity.txt"
        self.period_count.to_csv(out_txt, sep='\t', index=True)

    def filter_frame_shift(self):
        '''
        @Message  : filter the frame shifting genes.
        @Input    : param --> rpf density file
                          --> frame shifting threshold
                          --> sample names
        @Return   : output --> frame shifting genes
        @Flow     : step1 --> calculate the frame shifting ratio of each gene
                    step2 --> filter the frame shifting genes of each sample
        '''
        
        def filter_frame_ratio(row, period):
            # filter the frame_shifting, two conditions:
            # if the frame1 > frame 0 ,frame1 > frame2, frame1 > self.period, then the gene is frame1 shifting
            # if the frame2 > frame 0 ,frame2 > frame1, frame2 > self.period, then the gene is frame2 shifting
            if row[0] > row[1] and row[0] > row[2] and row[0] > period:
                return 'frame0'
            if row[1] > row[0] and row[1] > row[2] and row[1] > period:
                return 'frame1'
            elif row[2] > row[0] and row[2] > row[1] and row[2] > period:
                return 'frame2'
            else:
                return 'fuzzy'

        # calculate the ratio of each frame of samples
        for sp in self.sample_name:

            output = self.output + "_" + sp + "_gene_frame_shift.txt"
            sp_colname = [sp + '_f0', sp + '_f1', sp + '_f2']

            period_ratio = self.period_count[sp_colname].apply(lambda x: x / x.sum() * 100, axis=1)
            period_ratio.loc[:, 'frame'] = period_ratio.apply(filter_frame_ratio, period=self.period, axis=1)
            
            period_ratio.to_csv(output, sep='\t', index=False)

            gene_frame_count = period_ratio['frame'].value_counts().reindex(['frame0', 'frame1', 'frame2'], fill_value=0).reset_index()
            gene_frame_count.insert(0, 'Sample', sp)
            gene_frame_count.rename(columns={'frame': 'Frame', 'count': 'Count'}, inplace=True)

            self.gene_frame_shift = pd.concat([self.gene_frame_shift, gene_frame_count], axis=0)

    def output_frame_shift(self):
        out_txt = self.output + "_gene_frame_shift_count.txt"
        self.gene_frame_shift.to_csv(out_txt, sep='\t', index=False)

    def get_sub_plot_num(self):
        factor1 = int(self.sample_num ** 0.5)
        factor2 = int(self.sample_num ** 0.5)
        square = factor1 * factor2

        while not square >= self.sample_num:
            if factor2 > factor1:
                factor1 += 1
            else:
                factor2 += 1
            square = factor1 * factor2

        return factor1, factor2
    
    def draw_frame_shift_count(self):

        out_pdf = self.output + "_gene_frame_shift_count_plot.pdf"
        out_png = self.output + "_gene_frame_shift_count_plot.png"

        nrow, ncol = self.get_sub_plot_num()
        
        matplotlib.use('AGG')
        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(nrow * 4, ncol * 3), sharey=True)

        for (sample, group), ax in zip(self.gene_frame_shift.groupby('Sample'), axes.flatten()):
            bars = group.plot(x='Frame', y='Count', kind='bar', ax=ax, title=sample, legend=False)
            ax.set_ylabel('Count')
            ax.set_xlabel('')

            for bar in bars.patches:
                height = bar.get_height()
                ax.annotate(f'{height:.0f}', 
                            xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')
        
        plt.tight_layout()

        # plt.show()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png, dpi=300)
