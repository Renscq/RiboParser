#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : Cumulative_CoV.py


# import pandas as pd
# import polars as pl
# import numpy as np
# from collections import OrderedDict
# from Bio import SeqIO
# import argparse

import pandas as pd
from . import RPFs


class CumulativeCoV(object):

    def __init__(self, args):
        # input and output
        self.rpf = args.rpf
        self.out = args.output

        self.list = args.list
        self.gene_rpf = None

        # parameters for rpfs filter
        self.rpf_num = args.min
        self.trim = args.trim

        self.norm = args.normal
        self.zero = args.zero

        # format rpf table
        self.split = args.split
        self.format = True
        
        # rpf table
        self.sample_name = None
        self.sample_num = None
        self.total_rpf_num = None
        self.raw_rpf = None


    def retrieve_rpf(self):
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
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sample_num=None,
                                      sample_name=None,
                                      sites='P',
                                      frame='all',
                                      gene=self.list,
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


    def rpf_to_rpm(self):
        '''
        @Message  : convert rpf to rpm value.
        @Input    : high_rpf --> rpf table in pandas dataframe
                    norm --> default [False]
        @Return   : high_rpf --> rpm table
        '''
        
        if self.norm:
            total_rpf_num = [x for x in self.total_rpf_num for _ in range(3)]
            self.high_rpf.iloc[:, 6::] = self.high_rpf.iloc[:, 6:] * 1000000 / total_rpf_num
        else:
            pass

        self.high_rpf = self.high_rpf[self.high_rpf['region'] == 'cds']

        # set the cds start position to 1
        if self.zero:
            self.high_rpf['now_nt'] = self.high_rpf.groupby('name')['now_nt'].transform(lambda x: x - x.min() + 1)


    def melt_rpf_table(self):
        '''
        @Message  : melt the data into a format that is convenient for drawing.
        @Input    : self.high_rpf --> specific genes rpf table 
        @Return   : self.melt_rpf --> melt the three frame data to long data format

        from
            name ... wt_f0 wt_f1 wt_f2
            ORF1 ... 10 3   2
            ORF1 ... 7 2   2
        to
            name ... wt
            ORF1 ... 10
            ORF1 ... 3
            ORF1 ... 2
            ...
        '''

        ids_vars = self.high_rpf.columns[0:6].to_list()
        sample_title = [x + y for x in self.sample_name for y in ['_f0', '_f1', '_f2']]
        self.gene_rpf = pd.melt(self.high_rpf, id_vars=ids_vars, value_vars=sample_title, var_name='sample', value_name='rpf')
        
        self.gene_rpf['frame'] = self.gene_rpf['sample'].str[-1]
        self.gene_rpf['sample'] = self.gene_rpf['sample'].str[:-3]

        self.gene_rpf = self.gene_rpf.sort_values(by=['sample', 'name', 'now_nt', 'frame']).reset_index(drop=True)
        self.gene_rpf = self.gene_rpf.pivot(index=ids_vars + ['frame'], columns='sample', values='rpf').reset_index()
        self.gene_rpf['now_nt'] = self.gene_rpf['now_nt'].astype(int) + self.gene_rpf['frame'].astype(int)
        self.gene_rpf['codon'] = self.gene_rpf.apply(lambda row: row['codon'][int(row['frame'])], axis=1)


    def calc_cov(self):
        '''
        @Message  : calculate the coefficient of variation.
        @Input    : self.gene_rpf --> rpf table
        @Return   : self.gene_cov --> coefficient of variation table

        from

        gene	site	a	b
        x	1	4	5
        x	2	5	5
        x	3	10	1
        x	4	3	0
        x	5	10	3
        x	6	1	4
        x	7	9	1
        x	8	4	2
        x	9	7	3
        x	10	2	2
        y	1	5	3
        y	2	2	1
        y	3	3	4
        y	4	4	4
        y	5	5	1
        y	6	6	2
        y	7	4	2
        y	8	6	4
        y	9	7	1
        y	10	5	3
        y	11	3	0
        y	12	5	1

        to

        gene	site	a	b
        x	1	nan	0
        x	2	0.15713484	0
        x	3	0.507560566	0.514259477
        x	4	0.5652957	0.828221234
        x	5	0.52524176	0.728431359
        x	6	0.677867341	0.638284739
        x	7	0.608580619	0.702192845
        x	8	0.600656765	0.685118789
        x	9	0.553155291	0.637377439
        x	10	0.601497977	0.624926031
        y	1	0	0
        y	2	0.333333333	0.428571429
        y	3	0.40824829	0.374165739
        y	4	0.447213595	0.319438282
        y	5	0.471404521	0.306892205
        y	6	0.487950036	0.32249031
        y	7	0.5	0.30061372
        y	8	0.509175077	0.301018679
        y	9	0.516397779	0.319438282
        y	10	0.522232968	0.301647806
        y	11	0.527046277	0.316227766
        y	12	0.531085005	0.301511345

        '''
        
        # group the data by gene name
        gene_group = self.gene_rpf.groupby('name')

        # calculate the coefficient of variation
        def calc_cov(df, column):
            return df[column].expanding().std() / df[column].expanding().mean()

        # calculate the coefficient of variation for each gene
        gene_cov = gene_group.apply(calc_cov, self.sample_name)

        gene_cov = gene_cov.reset_index()

        self.gene_rpf[self.sample_name] = gene_cov[self.sample_name]


    def merge_cov_table(self):
        '''
        @Message  : merge the coefficient of variation table
        @Input    : self.gene_rpf --> rpf table
                    self.gene_cov --> coefficient of variation table
        '''
        # delete the gene with max(now_nt) < self.trim
        gene_rpf_list = self.gene_rpf[self.gene_rpf['now_nt'] > self.trim]['name'].unique().tolist()

        # trim the now_nt < self.trim
        gene_rpf_trim = self.gene_rpf[self.gene_rpf['name'].isin(gene_rpf_list)]
        gene_rpf_trim = gene_rpf_trim[self.gene_rpf['now_nt'] <= self.trim]
        gene_rpf_trim = gene_rpf_trim[['now_nt'] + self.sample_name]

        gene_rpf_meta = gene_rpf_trim.groupby('now_nt').mean().reset_index()
        gene_rpf_meta = pd.melt(gene_rpf_meta, id_vars='now_nt', value_vars=self.sample_name, 
                                var_name='Sample', value_name='Density')

        gene_rpf_meta['Frame'] = gene_rpf_meta['now_nt'].apply(lambda x: (x - 1) % 3 )
        gene_rpf_meta['Codon'] = gene_rpf_meta['now_nt'].apply(lambda x: (x - 1) // 3 + 1)
        gene_rpf_meta['Nucleotide'] = gene_rpf_meta['now_nt']
        gene_rpf_meta['Meta'] = "TIS"
        
        gene_rpf_meta = gene_rpf_meta[['Sample', 'Meta', 'Nucleotide', 'Codon', 'Frame', 'Density']]

        gene_rpf_meta.to_csv(self.out + "_meta" + str(self.trim) + "_cumulative_CoV.txt", sep="\t", index=False)


    def output_rpf_table(self):
        '''
        @Message  : split merged rpf table to each file.
        @Input    : self.melt_rpf --> description
        @Return   : output --> output rpf table in txt format
        '''
        
        if self.split:
            gene_group = self.gene_rpf.groupby('name')
            for gene, expr in gene_group:
                filename = gene
                expr.to_csv(filename + "_cumulative_CoV.txt", sep="\t", index=False)
        else:
            self.gene_rpf.to_csv(self.out + "_cumulative_CoV.txt", sep="\t", index=False)
