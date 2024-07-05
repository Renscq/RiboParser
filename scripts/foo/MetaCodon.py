#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : MetaCodon.py


# from keyword import kwlist
# from concurrent.futures import ThreadPoolExecutor

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from collections import OrderedDict
import seaborn as sns

from scipy.signal import savgol_filter

from . import RPFs


class MetaCodon(object):

    def __init__(self, args):
        # input and output
        self.ribo = args.rpf
        self.gene = args.list
        self.output = args.output

        # parameters for rpfs filter
        self.frame = args.frame
        self.rpf_num = args.min
        self.tis = args.tis
        self.tts = args.tts
        self.norm = args.normal

        # specified codon for meta-codon plot
        # self.thread = args.thread

        self.codon = args.codon
        self.codon_list = None
        self.around = args.around

        self.density_df = OrderedDict()
        self.sequence_df = OrderedDict()

        # format the meta-codon
        self.scale = args.scale
        self.unique = args.unique

        # smooth the meta-codon
        self.smooth = args.smooth

        # output the figures
        self.fig = args.fig

        # groups
        self.sample_name = None
        self.sample_num = None
        self.total_rpf_num = None
        self.high_rpf = None

    def import_codon(self):
        if self.codon:
            self.codon_list = pd.read_csv(self.codon,
                                          index_col=False,
                                          header=None,
                                          sep='\t')
            self.codon_list.columns = ['codon']
            self.codon_list['codon'] = self.codon_list['codon'].str.upper()
            self.codon_list['codon'] = self.codon_list['codon'].str.replace('U', 'T')
        else:
            codon_anno, codon_df = RPFs.codon_table()
            self.codon_list = codon_df.reset_index()

    def import_rpf(self):
        '''
        1. import the rpf density
        2. filter the high expressiong genes
        3. convert the rpf density to rpm density
        4. scale the rpf density with each gene expression

        '''
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.ribo,
                                      sites='A',
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis = self.tis,
                                      tts = self.tts,
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0]
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        self.high_rpf = rpf_results[7]

        del rpf_results

        if self.norm:
            self.merged_rpf[self.sample_name] = self.merged_rpf[self.sample_name] * 1e6 / self.total_rpf_num
            self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
            self.gene_rpf_sum = self.high_rpf.loc[:, ["name"] + self.sample_name].groupby("name")[self.sample_name].sum()

        if self.scale:
            # merge the gene expression with rpf density
            self.gene_rpf_sum = self.gene_rpf_sum.add_suffix("_sum")
            self.high_rpf = pd.merge(self.high_rpf, self.gene_rpf_sum, on = "name")

            # create the rpm dataframe and rpm_sum dataframe
            np.seterr(divide='ignore', invalid='ignore')
            rpm = self.high_rpf[self.sample_name]
            rpm_sum = self.high_rpf[self.gene_rpf_sum.columns]
            rpm_sum.columns = self.sample_name
            
            # scale the rpf density with each gene expression
            self.high_rpf[self.sample_name] = np.divide(rpm, rpm_sum)
            self.high_rpf = self.high_rpf.drop(columns = self.gene_rpf_sum.columns)

    def smooth_rpf_density(self):
        '''
        Smooth the rpf density with the mean of rpf density.

        RPFs density is a discrete distribution, which is not suitable for meta-codon plot.

        1. set the windows and power of smooth
        2. smooth the rpf with the mean of rpf density of each gene and sample
        '''
        
        self.smoothed_rpf = pd.DataFrame()

        if self.smooth:
            self.window = int(self.smooth.split(',')[0])
            self.power = int(self.smooth.split(',')[1])

            rpf_grouped = self.high_rpf.groupby('name')

            for gene, group in rpf_grouped:
                for colname in self.sample_name:
                    group[colname] = savgol_filter(group[colname], self.window, self.power)
                
                self.smoothed_rpf = pd.concat([self.smoothed_rpf, group])
        
            del self.high_rpf
            self.high_rpf = self.smoothed_rpf

        else:
            pass

    def delete_repeat_codon(self, codon_idx_true):
        '''
        Check the repeat codon in the range of specific codon window,
        Because, if the same codon in the same window, it will introduce the voice in the background.
        A clean background makes the results more reliable.

        1. delete the repeat codon in same window:
        [..., AAA, TGT, CAC, AAA, GGG, ATG, TGC, CTG, ...] ==> delete the 'AAA' codon site

        '''

        codon_idx_delta = (codon_idx_true.index[1:] - codon_idx_true.index[:-1]) > self.around * 2
        codon_idx_delta1_list = [True] + codon_idx_delta.tolist()
        codon_idx_delta2_list = codon_idx_delta.tolist() + [True]

        uniq_codon_list = [a and b for a, b in zip(codon_idx_delta1_list, codon_idx_delta2_list)]
        if len(uniq_codon_list) < 0:
            return 0, 0
        else:
            uniq_codon_idx = codon_idx_true.loc[uniq_codon_list]
            uniq_site_num = uniq_codon_idx.sum()

            return uniq_codon_idx, uniq_site_num

    def delete_out_of_range_site(self, uniq_codon_idx):
        '''
        reterieve the codon windows may out of the range of CDS region
        so we need to check the limits of codon windows.
        like the codon [GCT]:

        ACACAGTCGTACGTAGCA  GCT  AGCGAG[TGA]

        The GCT codon is located close to the stop codon, making it impossible to take a complete width sliding window centered on it.
        '''

        left_limit = uniq_codon_idx.index - self.around - 4
        right_limit = uniq_codon_idx.index + self.around + 4

        left_limit = left_limit[left_limit.isin(self.high_rpf.index)] + self.around + 4
        right_limit = right_limit[right_limit.isin(self.high_rpf.index)] - self.around - 4

        uniq_codon_idx = left_limit.intersection(right_limit)

        return uniq_codon_idx

    def scale_density(self, density_merge):
        '''
        scale the density with the mean of density
        '''
        if self.scale:
            density_merge = np.divide(density_merge, density_merge.mean())
            return density_merge
        else:
            return density_merge

    def single_codon(self, codon):
        '''
        match single codon and retrieve the density
        1. reterieve window: [-20, 0, 20]
        
        2. delete the repeat codon in same window:
        [..., AAA, TGT, CAC, AAA, GGG, ATG, TGC, CTG, ...] ==> delete the 'AAA' codon site
        
        3. reterieve the codon density by the window
        
        4. merge all results

        '''

        # retrieve the codon index and exclude the codon site out of range
        codon_idx = (self.high_rpf['codon'] == codon) & (self.high_rpf['from_tis'] > self.around) & (self.high_rpf['from_tts'] < -self.around)
        codon_idx_true = codon_idx[codon_idx == True]
        site_num = sum(codon_idx)
        
        # exit if none fitted results
        if site_num == 0:
            return

        # delete the Cross repetition codon in same window
        if self.unique:
            uniq_codon_idx, uniq_site_num = self.delete_repeat_codon(codon_idx_true)
        else:
            uniq_codon_idx = codon_idx_true
            uniq_site_num = site_num
        
        # exit if none fitted results
        if uniq_site_num == 0:
            return

        # delete the codon site out of range
        # uniq_codon_idx = self.delete_out_of_range_site(uniq_codon_idx)

        # exit if none fitted results
        if len(uniq_codon_idx) == 0:
            return

        # retrieve the codon density codon by codon
        density_df = []
        sequence_df = []
        
        # retrieve the codon in upstream/downstream window [-20, 0, 20]
        # for posi in range(-self.around, self.around + 1): # this is wrong, should be -21, -1, 19
        for posi in range(-self.around, self.around + 1):
            idx_tmp = uniq_codon_idx.index + posi
            tmp_dst = self.high_rpf.loc[idx_tmp, self.sample_name].mean()
            density_df.append(tmp_dst)

            tmp_seq = self.high_rpf.loc[idx_tmp, ['name', 'codon']]
            tmp_seq.index -= posi
            tmp_seq.columns = ['name', posi]
            sequence_df.append(tmp_seq)

        # merge meta RPFs
        density_merge = pd.concat(density_df, axis = 1).T
        density_merge = self.scale_density(density_merge)
        density_merge.index = np.arange(-self.around, self.around + 1)
        density_merge.insert(loc=0, column='Frame', value=self.frame)
        if self.frame == 'all':
            nt_index = density_merge.index * 3
            density_merge.insert(loc=0, column='Nucleotide', value=nt_index)
        else:
            nt_index = density_merge.index * 3 + int(self.frame)
            density_merge.insert(loc=0, column='Nucleotide', value=nt_index)

        self.density_df[codon] = [site_num, uniq_site_num, density_merge]

        # merge meta sequences
        sequence_merge = pd.concat(sequence_df, axis=1)
        sequence_merge = sequence_merge.loc[:, ~sequence_merge.columns.duplicated()]
        sequence_merge.reindex(columns=np.array(range(-self.around, self.around + 1)).tolist())
        self.sequence_df[codon] = [site_num, uniq_site_num, sequence_merge]

    def multiple_codon(self, codon, codon_num):
        '''
        match multiple codon and retrieve the density
        1. 
        reterieve first codon index: [4, 26, 79, 148, 205]

        reterieve second codon index - 1 : [5-1, 27-1, 89-1, 135-1, 206-1]
        
        2. multiple codons should then be aligned to the same index:
        
        aligend index is : [4, 26, 205]

        3. reterieve the codon window

        4. delete the repeat codon in same window:
        [..., AAAAAA, TGTTGT, CACCAC, AAAAAA, GGGGGG, ATGATG, TGCTGC, CTGCTG, ...] ==> delete the 'AAAAAA' codon site
        
        5. reterieve the codon density by the window
        
        6. merge all results

        '''

        print('Now is: {codon}.'.format(codon=codon), flush=True)

        # match first and second codon index
        idx_list = []
        shift_codon = 0
        for f0 in range(0, len(codon), 3):
            single_codon = codon[f0:f0 + 3]

            idx = (self.high_rpf['codon'] == single_codon) & (self.high_rpf['from_tis'] > self.around) & (self.high_rpf['from_tts'] < -self.around)
            idx = idx.shift(shift_codon, fill_value=False)
            shift_codon -= 1

            idx_list.append(idx)

        # merge multiple codon index
        idx_mer = pd.concat(idx_list, axis=1)
        idx = idx_mer.sum(axis=1)
        idx = idx.apply(lambda x: True if x == codon_num else False)
        site_num = sum(idx)

        # exit if none fitted results
        if site_num == 0:
            return

        # delete the Cross repetition codon in same window
        if self.unique:
            uniq_codon_idx, uniq_site_num = self.delete_repeat_codon(idx[idx])
        else:
            uniq_codon_idx = idx
            uniq_site_num = site_num
        
        # exit if none fitted results
        if uniq_site_num == 0:
            return

        # delete the codon site out of range
        # uniq_codon_idx = self.delete_out_of_range_site(uniq_codon_idx)
        
        # exit if none fitted results
        if len(uniq_codon_idx) == 0:
            return

        # retrieve the codon density codon by codon
        density_df = []
        sequence_df = []

        for posi in range(-self.around, self.around + codon_num):
            idx_tmp = uniq_codon_idx.index + posi
            
            tmp_dst = self.high_rpf.loc[idx_tmp, self.sample_name].mean()
            density_df.append(tmp_dst)

            tmp_seq = self.high_rpf.loc[idx_tmp, ['name', 'codon']]
            tmp_seq.index -= posi
            tmp_seq.columns = ['name', posi]
            sequence_df.append(tmp_seq)

        # merge meta RPFs
        density_merge = pd.concat(density_df, axis=1).T
        density_merge = self.scale_density(density_merge)

        density_merge.index = np.arange(-self.around, self.around + codon_num)
        density_merge.insert(loc=0, column='Frame', value=self.frame)
        if self.frame == 'all':
            nt_index = density_merge.index * 3
            density_merge.insert(loc=0, column='Nucleotide', value=nt_index)
        else:
            nt_index = density_merge.index * 3 + int(self.frame)
            density_merge.insert(loc=0, column='Nucleotide', value=nt_index)

        self.density_df[codon] = [site_num, uniq_site_num, density_merge]

        sequence_merge = pd.concat(sequence_df, axis=1)
        sequence_merge = sequence_merge.loc[:, ~sequence_merge.columns.duplicated()]
        self.sequence_df[codon] = [site_num, uniq_site_num, sequence_merge]
        
    def reterieve_codon_density(self):
        from concurrent.futures import ThreadPoolExecutor

        for codon in self.codon_list['codon']:
            print('Now is: {codon}.'.format(codon=codon), flush=True)
            # skip codon does not fit 3-nt periodicity
            if len(codon) % 3 != 0:
                print('Skip! {codon} does not fit 3-nt periodicity!'.format(codon=codon), flush=True)
                continue

            codon_num = int(len(codon) / 3)
            if codon_num == 1:
                self.single_codon(codon)
            elif codon_num > 1:
                self.multiple_codon(codon, codon_num)
            else:
                print('Skip! {codon} does not fit 3-nt periodicity!'.format(codon=codon), flush=True)

    def output_meta_codon_density(self):
        for codon, mess in self.density_df.items():
            mess[2].index.name = 'Codon'
            mess[2].to_csv(self.output + '_' + codon + '_' + str(mess[0]) + '_' + str(mess[1]) + '_meta_density.txt', sep='\t', index=True)
            mess[2] = mess[2].drop(columns = ['Nucleotide', 'Frame'])
    
    def output_meta_codon_seq(self):
        for codon, mess in self.sequence_df.items():
            mess[2].index.name = 'Codon'
            mess[2].to_csv(self.output + '_' + codon + '_' + str(mess[0]) + '_' + str(mess[1]) + '_meta_sequence.txt', sep='\t', index=True)
    
    def draw_meta_codon(self):
        # draw the line meta plot of each codon
        if self.fig:

            for codon, mess in self.density_df.items():
                out_pdf = self.output + '_' + codon + '.pdf'
                out_png = self.output + '_' + codon + '.png'

                # if self.scale:
                #     mess[2] = mess[2]/mess[2].mean()

                matplotlib.use('AGG')
                fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
                ax = sns.lineplot(data=mess[2], linewidth=1.2, ax=ax)

                # ax.set_xlim(left=-20, right=20)
                ax.autoscale(enable=True, axis='y')

                plt.legend(loc='center left',
                        bbox_to_anchor=(1, 0.5),
                        fancybox=True,
                        shadow=True,
                        ncol=1)

                ax.set_title(codon + '(' + str(mess[0]) + ', ' + str(mess[1]) + ')')
                ax.set_xlabel('Distance from codon')

                # set the ylabel
                if self.norm:
                    if self.scale:
                        ax.set_ylabel('Scaled mean RPM density')
                    else:
                        ax.set_ylabel('Mean RPM density')
                else:
                    if self.scale:
                        ax.set_ylabel('Scaled mean RPFs density')
                    else:
                        ax.set_ylabel('Mean RPFs density')

                plt.tight_layout()
                # plt.show()

                fig.savefig(fname=out_pdf)
                fig.savefig(fname=out_png)
                plt.close()
        else:
            pass
