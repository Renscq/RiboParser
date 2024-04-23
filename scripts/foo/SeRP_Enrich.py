#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : SeRP_enrich.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/04/20 17:00:28
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from . import RPFs


class SeRP_Enrich(object):

    def __init__(self, args):
        # input and output files
        self.rpf = args.rpf
        self.ip_samples = args.ip
        self.wt_samples = args.wt

        self.output_prefix = args.output
        self.all_enrich = args.all

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

        # opts for enrich analysis
        self.smooth = args.smooth
        self.enrich = None
        self.meta_enrich = None


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
        self.high_rpf = rpf_results[7].copy()
        
        del rpf_results

        self.merged_rpf = self.merged_rpf[self.merged_rpf['codon'].isin(self.codon_dict.keys())]

        if self.norm:
            self.merged_rpf[self.sample_name] = self.merged_rpf[self.sample_name] * 1e6 / self.total_rpf_num
            # self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
            # self.gene_rpf_sum = self.high_rpf.loc[:, ["name"] + self.sample_name].groupby("name")[self.sample_name].sum()

        self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
        del self.merged_rpf
