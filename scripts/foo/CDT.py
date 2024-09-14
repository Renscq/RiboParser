#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : CDT.py


import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from . import RPFs


class CodonDecodingTime(object):

    def __init__(self, args):
        self.mrna_file = args.list
        self.gene_table = None

        self.rpf_file = args.rpf
        self.rna_file = args.rna

        self.sample_num = 0
        self.sample_name = None
        self.rel_name = None
        self.merged_rpf = None
        self.merged_rpf_na = None

        # for the high expression gene filter
        # self.norm = args.normal
        self.min = args.min

        # rpf table
        self.rpf_sample = None
        self.total_rpf_num = None
        self.rpf_gene_sum = None
        self.rpf_gene = None
        self.high_rpf = None

        # rna table
        self.rna_sample = None
        self.total_rna_num = None
        self.rna_gene_sum = None
        self.rna_gene = None
        self.high_rna = None

        # filter the stable codon coverage range
        self.scale = args.scale
        self.site = args.site
        self.frame = args.frame
        self.tis = args.tis
        self.tts = args.tts

        # table for cdt
        self.rpf_norm_dst = None
        self.cdt = None

        self.stop = args.stop
        self.codon_dict, self.codon_table = RPFs.codon_table()

        if self.stop:
            # remove the TAA TAG TGA from the codon table
            del self.codon_dict['TAA']
            del self.codon_dict['TAG']
            del self.codon_dict['TGA']

        self.output = args.output


    def import_gene(self):
        '''
        @Message  : import the gene file
        @Input    : self.gene --> gene file
        @Return   : self.gene --> gene dataframe
        @Flow     : step1 --> import the gene file
        '''

        self.gene_table = pd.read_csv(self.mrna_file, sep='\t', index_col=False)
        self.gene_table['length'] = self.gene_table['utr5_length'] + self.gene_table['cds_length'] + self.gene_table['utr3_length']
        self.gene_table.set_index('transcript_id', inplace=True)
        self.gene_table.index.name = 'name'


    def import_density(self):
        """
        import the Ribo-seq and RNA-seq density file.
        By default, the two data names are different.

        return the main results as follows:
        samples name,
        high expression gene density
        high expression gene name
        total rpf/reads
        """

        # import the rpf density
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file, sample_num=None,
                                      sample_name=None, sites=self.site, frame=self.frame,
                                      gene=self.mrna_file, rpf_num=self.min, tis=self.tis, tts=self.tts)
        self.rpf_sample = rpf_results[1]
        self.total_rpf_num = rpf_results[4]
        self.rpf_gene_sum = rpf_results[5]
        self.rpf_gene = rpf_results[6]
        self.high_rpf = rpf_results[7].copy()

        del rpf_results

        print('', flush=True)

        # import the rna density
        rna_results = RPFs.import_rpf(rpf_file=self.rna_file, sample_num=None,
                                      sample_name=None, sites=self.site, frame=self.frame,
                                      gene=self.mrna_file, rpf_num=self.min, tis=self.tis, tts=self.tts)
        self.rna_sample = rna_results[1]
        self.total_rna_num = rna_results[4]
        self.rna_gene_sum = rna_results[5]
        self.rna_gene = rna_results[6]
        self.high_rna = rna_results[7].copy()

        del rna_results


    def merge_rpf_rna(self):
        """
        get the overlap gene from rpf and mrna density file
        only these genes will submit to follow analysis
        """
        
        self.overlap_gene = set(self.rpf_gene) & set(self.rna_gene)

        self.high_rpf = self.high_rpf.loc[self.high_rpf['name'].isin(self.overlap_gene)]
        self.high_rna = self.high_rna.loc[self.high_rna['name'].isin(self.overlap_gene)]


    def normalize_density(self):
        """
        normalize the density to rpkm. rpkm = (10^9 * C) / (N * L)
        N: total reads number, L: gene length, C: reads number of the gene
        """
        # filter the rpf_gene_sum index with the self.overlap_gene
        gene_len = self.gene_table.loc[list(self.overlap_gene),'length']

        self.rpf_gene_sum = self.rpf_gene_sum.loc[list(self.overlap_gene)]
        self.rpf_gene_rpm = 1e6 * self.rpf_gene_sum / self.total_rpf_num
        self.rpf_gene_rpkm = 1e9 * self.rpf_gene_sum.div(gene_len, axis=0) / self.total_rpf_num

        self.rna_gene_sum = self.rna_gene_sum.loc[list(self.overlap_gene)]
        self.rna_gene_rpm = 1e6 * self.rna_gene_sum / self.total_rna_num
        self.rna_gene_rpkm = 1e9 * self.rna_gene_sum.div(gene_len, axis=0) / self.total_rna_num


    @staticmethod
    def scale_method(scale, cdt):
        if scale == 'minmax':
            # relative_cdt = (cdt - cdt.min()) / (cdt.max() - cdt.min())
            relative_cdt = cdt / cdt.max()
        elif scale == 'zscore':
            relative_cdt = (cdt - cdt.mean()) / cdt.std()
        else:
            print('Unknown scale method.', flush=True)
            sys.exit()
        return relative_cdt
    

    def calc_rpf_norm(self):
        '''
        @Message  : function for rpf normalization.
        @Input    : self.high_rpf --> high expression gene rpf density
        @Return   : self.rpf_norm_dst --> normalized rpf density with RNA-seq RPKM
        @Flow     : step1 --> calculate the average density and normalized rpf density of each codon
                    step2 --> merge average density 
        '''
        
        # calculate the average density of each codon
        self.high_rpf = self.high_rpf[self.high_rpf['region'] == 'cds']
        high_rpf_group = self.high_rpf.groupby('name')
        now_num = 0
        rpf_norm_list = []

        self.count_column = [i + '_count' for i in self.rpf_sample]
        self.norm_column = [i + '_norm' for i in self.rpf_sample]

        for gene, rpf in high_rpf_group:
            # counter
            now_num += 1
            if now_num % 500 == 0:
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)
            if now_num == len(high_rpf_group):
                print('Row: {rows}, {gene}.'.format(rows=now_num, gene=gene), flush=True)

            # sum the rpf and each codon number
            rpf_sum = rpf.groupby('codon')[self.rpf_sample].sum()
            rpf_norm = rpf_sum.div(self.rna_gene_rpkm.loc[gene].to_list(), axis=1)

            rpf[self.rpf_sample] = rpf[self.rpf_sample].map(lambda x: 1 if x != 0 else 0)
            codon_sum = rpf.groupby('codon')[self.rpf_sample].sum()

            rpf_sum.columns = self.count_column
            rpf_norm.columns = self.norm_column

            temp_rpf_norm = pd.concat([codon_sum, rpf_sum, rpf_norm], axis=1)
            rpf_norm_list.append(temp_rpf_norm)

        self.rpf_norm_dst = pd.concat(rpf_norm_list, axis=0)
                
        # check the stop codon
        if self.stop:
            # if cdt table contains the TAA TAG TGA, remove them
            self.rpf_norm_dst = self.rpf_norm_dst.drop(['TAA', 'TAG', 'TGA'], errors='ignore')

        # drop the na and inf values
        self.rpf_norm_dst.replace([np.inf, -np.inf], np.nan, inplace=True)

        self.rpf_norm_dst.reset_index('codon', inplace=True)

    def codon_decoding_time(self):
        '''
        @Message  : function for codon decoding time calculatioin.
        @Input    : param --> description
        @Return   : output --> description
        @Flow     : step1 --> merge all norm
                    step2 --> calculate cdt
                    step3 --> calculate relative cdt
                    step4 --> output the cdt table
        '''

        # merge average density and calculate cdt
        sum_rpf_norm = self.rpf_norm_dst.groupby('codon').sum()

        # calculate the absolute cdt
        cdt = sum_rpf_norm[self.count_column] / sum_rpf_norm[self.rpf_sample].rename(columns=dict(zip(self.rpf_sample, self.count_column)))
        norm_cdt = sum_rpf_norm[self.norm_column] / sum_rpf_norm[self.rpf_sample].rename(columns=dict(zip(self.rpf_sample, self.norm_column)))

        cdt.columns = [i + '_absolute_cdt' for i in self.rpf_sample]
        norm_cdt.columns = [i + '_norm_cdt' for i in self.rpf_sample]

        # calculate the relative cdt
        relative_cdt = self.scale_method(self.scale, norm_cdt)
        relative_norm_cdt = self.scale_method(self.scale, relative_cdt)

        relative_cdt.columns = [i + '_relative_cdt' for i in self.rpf_sample]
        relative_norm_cdt.columns = [i + '_norm_relative_cdt' for i in self.rpf_sample]

        # merge all the cdt table
        self.codon_table = pd.DataFrame(self.codon_dict).T
        self.codon_table.columns = ['AA', 'Abbr']
        self.codon_table = pd.concat([self.codon_table, sum_rpf_norm, cdt, relative_cdt, norm_cdt, relative_norm_cdt], axis=1)
        self.codon_table.sort_values('Abbr', inplace=True)
        self.codon_table.index.name = 'Codon'
        self.codon_table.reset_index(inplace=True)

        self.codon_table.to_csv(self.output + '_cdt.txt', header=True, index=False, sep='\t')


    def draw_cdt_corr(self):
        out_pdf = self.output + "_cdt_corrplot.pdf"
        out_png = self.output + "_cdt_corrplot.png"

        occ_corr = self.codon_table.filter(regex='_norm_cdt').astype('float').corr()
        occ_corr.to_csv(self.output + '_cdt_corr.txt', sep='\t', index=True)

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6 + 0.2 * len(self.rpf_sample),
                                  8 + 0.2 * len(self.rpf_sample)),
                                  dpi=300)
        sns.heatmap(data=occ_corr, 
                    annot=True, 
                    fmt=".4f", 
                    cmap="vlag",
                    linewidths=.01)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()

    def draw_cdt_heat(self):
        # output the absolute occupancy heatmap
        out_pdf = self.output + "_cdt_heatplot.pdf"
        out_png = self.output + "_cdt_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(10 + 0.5 * len(self.rpf_sample), 20), dpi=300)
        sns.heatmap(data=self.codon_table.filter(regex='_norm_cdt').astype('float'),
                    yticklabels=self.codon_table['Codon'] + '[' + self.codon_table['Abbr'].astype(str) + ']',
                    annot=True, 
                    fmt=".4f", 
                    cmap="vlag", 
                    linewidths=0)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()
