#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : Shuffle.py


import pandas as pd
import polars as pl
import numpy as np
from . import RPFs


class Shuffle(object):

    def __init__(self, args):
        self.rpf_file = args.rpf
        
        self.gene = args.list

        self.shuffle_rpf = None
        self.seed = args.seed
        self.individual = args.individual

        self.output = args.output + '_shuffle.txt'

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
                                      rpf_num=0)

        self.raw_rpf = rpf_results[0].to_pandas()
        self.sample_name = rpf_results[1]
        # self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        # self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        # self.high_gene = rpf_results[6]

        del rpf_results


    def shuffle_rpfs(self):
        '''
        @Message  : shuffle the rpf table.
        @Input    : self.raw_rpf --> rpf table contain all rpf value
                    self.seed --> random seed
        @Return   : self.shuffle_rpf --> shuffled rpf table
        '''

        sample_name = []
        for now_sp in self.sample_name:
            sample_name.extend([now_sp + '_f0', now_sp + '_f1', now_sp + '_f2'])

        raw_rpf_group = self.raw_rpf.groupby('name')
        
        if self.individual:
            raw_rpf_group_list = []

            now_num = 0
            for idx, now_group in raw_rpf_group:
                now_num += 1
                if now_num % 100 == 0:
                    print(f'Processing {now_num} genes.', flush=True)

                rpf_melt = pd.melt(now_group, id_vars=['name', 'now_nt', 'from_tis', 'from_tts', 'region', 'codon'], 
                                   value_vars=sample_name, 
                                   var_name='Frame', value_name='Count')
                rpf_melt['Sample'] = rpf_melt['Frame'].str[:-3]

                np.random.seed(self.seed)
                rpf_melt['Count'] = rpf_melt.groupby('Sample')['Count'].transform(lambda x: np.random.permutation(x))
                rpf_melt.drop('Sample', axis=1, inplace=True)

                rpf_melt = rpf_melt.pivot_table(index=['name', 'now_nt', 'from_tis', 'from_tts', 'region', 'codon'], 
                                     columns='Frame', values='Count', aggfunc='first').reset_index()

                raw_rpf_group_list.append(rpf_melt)

            print(f'Processing {now_num} genes.', flush=True)

        else:
            raw_rpf_group_list = []

            now_num = 0
            for idx, now_group in raw_rpf_group:
                now_num += 1
                if now_num % 100 == 0:
                    print(f'Processing {now_num} genes.', flush=True)

                rpf_melt = pd.melt(now_group, id_vars=['name', 'now_nt', 'from_tis', 'from_tts', 'region', 'codon'], 
                                   value_vars=sample_name, 
                                   var_name='Frame', value_name='Count')
                rpf_melt['Sample'] = rpf_melt['Frame'].str[:-3]
                rpf_melt['site'] = rpf_melt['now_nt'] + rpf_melt['Frame'].str[-1].astype(int)

                rpf_melt_count = rpf_melt.pivot_table(index='site', columns='Sample', values='Count', aggfunc='first').reset_index()

                # shuffle the rpf value
                np.random.seed(self.seed)
                rpf_melt_count_index = np.random.permutation(rpf_melt_count.index)
                rpf_melt_count = rpf_melt_count.reindex(rpf_melt_count_index).reset_index(drop=True)
                rpf_melt_count.drop('site', axis=1, inplace=True)
                rpf_melt_count = rpf_melt_count.melt(var_name='Sample', value_name='Count')
                
                rpf_melt['Count'] = rpf_melt_count['Count']
                rpf_melt.drop(['Sample', 'site'], axis=1, inplace=True)

                rpf_melt = rpf_melt.pivot_table(index=['name', 'now_nt', 'from_tis', 'from_tts', 'region', 'codon'], 
                                     columns='Frame', values='Count', aggfunc='first').reset_index()

                raw_rpf_group_list.append(rpf_melt)

            print(f'Processing {now_num} genes.', flush=True)

        self.shuffle_rpf = pd.concat(raw_rpf_group_list, axis=0)
        self.shuffle_rpf.reset_index(drop=True, inplace=True)


    def output_rpfs(self):
        '''
        @Message  : output the shuffled rpf table.
        @Input    : self.shuffle_rpf --> shuffled rpf table
                    self.output --> output file name
        @Return   : None
        '''
        
        self.shuffle_rpf = pl.DataFrame(self.shuffle_rpf)
        self.shuffle_rpf.write_csv(self.output, separator='\t')
