#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : Periodicity.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/01/22 14:13:22
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

class Periodicity(object):

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
        self.raw_rpf = None
        self.high_gene = None

        self.period = None

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
        # self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        # self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        # self.high_rpf = rpf_results[7]

        del rpf_results

    def calc_3nt_period(self):
        '''
        @Message  : function for 3nt periodicity calculation.
        @Input    : param --> rpf density file
        @Return   : output --> summary of 3nt periodicity
        @Flow     : step1 --> run
        '''

        raw_rpf = self.raw_rpf.drop(columns=["name", "now_nt", "from_tis", "from_tts", "region", "codon"])
        
        period = raw_rpf.apply(np.sum, axis=0).to_frame()
        period = period.reset_index()
        period.columns = ["Name", "Count"]
        period["Sample"] = period['Name'].str[:-3]
        period["Frame"] = period['Name'].str[-1::]
        period["Ratio"] = period["Count"] / period.groupby("Sample")["Count"].transform("sum") * 100

        self.period = period[["Sample", "Frame", "Count", "Ratio"]]

    def output_meta(self):
        out_txt = self.output + "_periodicity.txt"
        self.period.to_csv(out_txt, sep='\t', index=False)

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
    
    def draw_3nt_period_count(self):

        out_pdf = self.output + "_count_periodicity_plot.pdf"
        out_png = self.output + "_count_periodicity_plot.png"

        nrow, ncol = self.get_sub_plot_num()
        
        matplotlib.use('AGG')
        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(nrow * 4, ncol * 3), sharey=True)

        # FIX: ensure axes is always an array
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        axes = axes.flatten()

        for (sample, group), ax in zip(self.period.groupby('Sample'), axes.flatten()):
            bars = group.plot(x='Frame', y='Count', kind='bar', ax=ax, title=sample, legend=False)
            ax.set_ylabel('Count')
            ax.set_xlabel('Frame')

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

    def draw_3nt_period_ratio(self):

        out_pdf = self.output + "_ratio_periodicity_plot.pdf"
        out_png = self.output + "_ratio_periodicity_plot.png"

        nrow, ncol = self.get_sub_plot_num()

        matplotlib.use('AGG')
        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(nrow * 4, ncol * 3), sharey=True)

        # FIX: ensure axes is always an array
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        axes = axes.flatten()
        
        for (sample, group), ax in zip(self.period.groupby('Sample'), axes.flatten()):
            bars = group.plot(x='Frame', y='Ratio', kind='bar', ax=ax, title=sample, legend=False)
            ax.set_ylabel('Ratio (%)')
            ax.set_xlabel('Frame')

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
