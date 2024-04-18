#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : SeRP_Seq.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/03/22 22:05:46
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import pandas as pd
import polars as pl
# import numpy as np
# from collections import OrderedDict
from Bio import SeqIO


class SeRPSeq(object):
    def __init__(self, args):
        self.rpf_file = args.rpf
        self.rpf_seq = None

        self.peak_file = args.peak
        self.peaks = None
        self.peaks_gene = None

        self.upstream = args.upstream
        self.downstream = args.downstream
        
        self.pvalue = args.pvalue
        self.mean_fold = args.fold
        
        self.format = args.format
        self.output = args.output


    def import_peak_region(self):
        '''
        @Message  : function for import the bed file and filter the significant peaks
        @Input    : self.rpf_file --> input rpf file
                    self.pvalue --> pvalue threshold
                    self.mean_fold --> mean fold threshold
        @Return   : output --> sig_peak and sig_peak_gene
        @Flow     : step1 --> read the peaks file
                    step2 --> filter the significant peaks
                    step3 --> return the significant peaks and gene
        '''
        
        self.peaks = pd.read_csv(self.rpf_file, sep='\t')
        self.peaks = self.peaks.loc[self.peaks.iloc['P_Value'] != '-', ]
        self.peaks = self.peaks.apply(pd.to_numeric, errors='ignore')

        self.peaks = self.peaks[self.peaks['P_Value'] < self.pvalue]
        self.peaks = self.peaks[self.peaks['mean_fold'].replace('-', 'nan') >= self.mean_fold]

        self.peaks['upstream'] = self.peaks['peak_start'] - self.upstream
        self.peaks['downstream'] = self.peaks['peak_end'] + self.downstream

        self.peaks = self.peaks.reset_index(drop=True)
        self.peaks_gene = self.peaks['transcripts'].unique()


    def import_rpf(self):
        '''
        @Message  : function for import the rpf file
        @Input    : self.rpf_file --> input rpf file
        @Return   : self.rpf_seq --> rpf file only contains the sequence and gene message
        @Flow     : step1 --> import the rpf txt file
                    step2 --> filter the gene in peaks
        '''
        
        self.rpf_seq = pl.read_csv(self.rpf_file)
        self.rpf_seq = self.rpf_seq.to_pandas()
        
        self.rpf_seq = self.rpf_seq.iloc[:, 0:6]
        self.rpf_seq = self.rpf_seq[self.rpf_seq['name'].isin(self.peaks_gene)]
    

    def retrieve_seq(self):
        '''
        @Message  : function for retrieve the sequence of peak region
        @Input    : self.peaks --> input peak file
                    self.rpf_seq --> input rpf file
                    self.upstream --> upstream length
                    self.downstream --> downstream length
        @Return   : self.output --> output file name
        @Flow     : step1 --> add the [upstream, peak_start, peak_end, downstream] region of peak to the rpf table
                    step2 --> groupby the region and genes
                    step3 --> retrieve the sequence of the region
                    step4 --> output the sequence
        '''
        
        # add the [upstream, peak_start, peak_end, downstream] region of peak to the rpf table
        '''
        apply the function to the dataframe, to add a column contains annotation ['upstream', 'peak', 'downstream']
        if upstream < from_tis < peak_start:
            annotation is 'upstream'
        if peak_start < from_tis < peak_end:
            annotation is 'peak'
        if peak_end < from_tis < downstream:
            annotation is 'downstream'

        '''
        
        self.rpf_seq['annotation'] = self.rpf_seq.apply(lambda x: 'upstream' if x['from_tis'] < x['peak_start'] else 'peak' 
                                                        if x['peak_start'] < x['from_tis'] < x['peak_end'] else 'downstream', axis=1)


    def output_seq(self):

        pass
