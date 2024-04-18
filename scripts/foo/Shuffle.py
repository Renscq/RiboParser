#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : Shuffle.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/04/09 17:44:20
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import pandas as pd
import polars as pl
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
import argparse
from . import RPFs


def Shuffle(args):

    def __init__(self, args):
        self.rpf_file = args.rpf
        
        self.gene = args.list

        self.shuffle_rpf = None
        self.seed = args.seed

        self.output = args.output

    def import_rpf(self):
        '''
        @Message  : import the rpf table.
        @Input    : self.rpf --> rpf table
                    self.site --> [E, P, A] site in ribosome
                    self.frame --> [0, 1, 2] frame in [E, P, A] site
                    self.gene --> gene list
                    self.rpf_num --> minimum rpf value
                    norm --> normalize the rpf value to RPM
        @Return   : self.raw_rpf --> rpf table contain all rpf value
                    self.sample_name --> samples name
                    self.sample_num --> groups sample number
                    self.merged_rpf --> merged rpf table
                    self.total_rpf_num --> total rpf number of each sample
        '''
        
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file,
                                      sample_num=None,
                                      sample_name=None,
                                      sites='P',
                                      frame='all',
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)

        self.raw_rpf = rpf_results[0].to_pandas()
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]

        self.high_rpf = self.raw_rpf.loc[self.raw_rpf['name'].isin(self.high_gene), ]

        del rpf_results

    def shuffle_rpf(self):
        '''
        @Message  : shuffle the rpf table.
        @Input    : self.raw_rpf --> rpf table contain all rpf value
                    self.seed --> random seed
        @Return   : self.shuffle_rpf --> shuffled rpf table
        '''
        np.random.seed(self.seed)
        self.shuffle_rpf = self.raw_rpf.sample(frac=1).reset_index(drop=True)
