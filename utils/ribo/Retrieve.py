#!/usr/bin/env python
# @Project : riboparser
# @Script  : Retrieve.py

import pandas as pd
from . import RPFs


class Retrieve(object):

    def __init__(self, args):
        # input and output
        self.rpf = args.rpf
        self.out = args.output

        self.list = args.list
        self.gene_rpf = None

        # parameters for rpfs filter
        self.rpf_num = args.min
        self.norm = args.normal

        # format rpf table
        self.split = args.split
        self.format = args.format
        
        # rpf table
        self.sample_name = None
        self.sample_num = None
        self.total_rpf_num = None
        self.raw_rpf = None

    def import_gene_list(self):
        '''
        @Message  : import the gene list.
        @Input    : self.list --> gene list in one column
        @Return   : self.gene --> gene list, [YDL246C, YDL243C ...]
        '''
        if self.list is None:
            self.gene = None
        else:
            gene_csv = pd.read_csv(self.list, header=None, names=None)
            self.gene = gene_csv.iloc[:, 0].to_list()

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
            #　print the class of data type of eaech column
            now_cols = self.high_rpf.columns[6:].to_list()
            self.high_rpf[now_cols] = self.high_rpf[now_cols].apply(pd.to_numeric, downcast='float')
            self.high_rpf.iloc[:, 6:] = (self.high_rpf.iloc[:, 6:] * 1000000) / total_rpf_num
        else:
            pass

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
        
        if self.format:
            ids_vars = self.high_rpf.columns[0:6].to_list()
            sample_title = [x + y for x in self.sample_name for y in ['_f0', '_f1', '_f2']]
            self.gene_rpf = pd.melt(self.high_rpf, id_vars=ids_vars, value_vars=sample_title, var_name='sample', value_name='rpf')
            
            self.gene_rpf['frame'] = self.gene_rpf['sample'].str[-1]
            self.gene_rpf['sample'] = self.gene_rpf['sample'].str[:-3]

            self.gene_rpf = self.gene_rpf.sort_values(by=['sample', 'name', 'now_nt', 'frame']).reset_index(drop=True)
            self.gene_rpf = self.gene_rpf.pivot(index=ids_vars + ['frame'], columns='sample', values='rpf').reset_index()
            self.gene_rpf['now_nt'] = self.gene_rpf['now_nt'].astype(int) + self.gene_rpf['frame'].astype(int)
            self.gene_rpf['codon'] = self.gene_rpf.apply(lambda row: row['codon'][int(row['frame'])], axis=1)
            
        else:
            self.gene_rpf = self.high_rpf

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
                expr.to_csv(filename + "_retrieve.txt", sep="\t", index=False)
        else:
            self.gene_rpf.to_csv(self.out + "_retrieve.txt", sep="\t", index=False)
