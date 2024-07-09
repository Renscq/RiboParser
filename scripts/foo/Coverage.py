#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Coverage.py


from concurrent.futures import ThreadPoolExecutor
from math import ceil
from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from . import RPFs


class Coverage(object):

    def __init__(self, args):
        # input and output
        self.rpf = args.rpf
        self.gene = args.transcript
        self.gene_table = None
        self.output = args.output

        # filter the stable codon coverage range
        self.tis = None
        self.tts = None
        self.frame = args.frame
        self.site = 'all'
        self.rpf_num = args.min
        bins = args.bin.split(',')
        self.utr5_bin = int(bins[0])
        self.cds_bin = int(bins[1])
        self.utr3_bin = int(bins[2])

        self.norm = args.normal
        self.outlier = args.outlier
        self.thread = args.thread
        self.set = args.set
        self.sample_name = None

        # for the high expression gene filter
        self.raw_rpf = None
        self.merged_rpf = None
        # self.high_gene = None
        self.high_rpf = None
        self.gene_num = None
        self.total_rpf_num = None

        # mean of coverage
        self.utr5_dict = OrderedDict()
        self.cds_dict = OrderedDict()
        self.utr3_dict = OrderedDict()
        self.utr5_mean = None
        self.cds_mean = None
        self.utr3_mean = None

    def read_rpf(self):
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sites=self.site,
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis=self.tis,
                                      tts=self.tts,
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0]
        self.sample_name = rpf_results[1].copy()
        self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        self.gene_num = len(self.high_gene)
        self.high_rpf = rpf_results[7].copy()
        del rpf_results

        if self.norm:
            self.high_rpf.loc[:, self.sample_name] = self.high_rpf[self.sample_name] * 1e6 / self.total_rpf_num


    def import_gene(self):
        '''
        @Message  : import the gene file
        @Input    : self.gene --> gene file
        @Return   : self.gene --> gene dataframe
        @Flow     : step1 --> import the gene file
        '''

        self.gene_table = pd.read_csv(self.gene, sep='\t', index_col=False)
        self.gene_table = self.gene_table.loc[self.gene_table['transcript_id'].isin(self.high_gene), ]


    @staticmethod
    def run_fill_outliers(args):

        def iqr_fill_outliers(groups):
            q1 = groups.quantile(0.25)
            q3 = groups.quantile(0.75)
            iqr = q3 - q1
            lower_bound = q1 - 2 * iqr
            upper_bound = q3 + 2 * iqr
            groups[(groups < lower_bound) | (groups > upper_bound)] = groups.mean()
            return groups

        def zscore_fill_outliers(groups):
            z_scores = (groups - groups.mean()) / groups.std()
            threshold = 5
            groups[abs(z_scores) > threshold] = groups.mean()
            # groups[abs(z_scores) > threshold] = groups.mean()
            return groups
        
        gene_df = args
        return gene_df.groupby('Gene').apply(lambda x: x.apply(iqr_fill_outliers))
    
    
    @staticmethod
    def filter_outliers(self, gene_coverage):
        '''
        @Message  : filter the outliers
        @Input    : self.utr5_mean, self.cds_mean, self.utr3_mean --> the mean of the gene rpfs dataframe
        @Return   : self.utr5_mean, self.cds_mean, self.utr3_mean --> the mean of the gene rpfs dataframe
        @Flow     : step1 --> calculate the z-score
                    step2 --> filter the outliers
        '''

        gene_cover_df_n = gene_coverage.drop(columns=['Bins'])
        
        # stat the gene number
        gene_list = list(self.high_gene)
        gene_list_len = len(self.high_gene)

        gene_list_split = [gene_list[i:i + gene_list_len // self.thread] for i in range(0, gene_list_len, ceil(gene_list_len / self.thread))]
        gene_cover_split = [gene_cover_df_n[gene_cover_df_n['Gene'].isin(gene_list_split[i])] for i in range(self.thread)]

        args = [(gene_cover_split[i]) for i in range(self.thread)]

        # run the fill_outliers function with multi-thread
        from multiprocessing import Pool
        pool = Pool(processes=self.thread)
        gene_cover_list = pool.map(self.run_fill_outliers, args)
        pool.close()
        pool.join()

        # gene_cover_df_n = gene_cover_df_n.groupby('Gene').apply(lambda x: x.apply(fill_outliers))
        gene_cover_df_n = pd.concat(gene_cover_list, axis=0)
        gene_cover_df_n.reset_index(inplace=True)
        gene_cover_df_n.rename(columns={'level_1': 'Bins'}, inplace=True)
        gene_cover_df_n['Bins'] = gene_coverage['Bins']
        gene_coverage = gene_cover_df_n.groupby('Bins')[self.sample_name].sum() / self.gene_num
        
        return gene_coverage


    def adjust_dimension(self, rpm, bins):
        '''
        @Message  : adjust the dimension of the gene rpfs dataframe
        @Input    : rpm --> gene rpfs dataframe contain the high expression genes
        @Return   : self.utr5_dict, self.cds_dict, self.utr3_dict --> the mean of the gene rpfs dataframe
                    self.utr5_mean, self.cds_mean, self.utr3_mean --> the mean of the gene rpfs dataframe
        @Flow     : step1 --> set the bins
                    step2 --> calculate the mean of the each gene rpfs
                    step3 --> calculate the mean of the gene rpfs in 5-UTR, CDS, 3-UTR
        '''
        
        region_dict = OrderedDict()
        
        for idx, gene in rpm.groupby('name'):
            gene.loc[:, 'Bins'] = pd.cut(gene.index, bins, labels=range(bins))
            # region_dict[idx] = gene.groupby('Bins')[self.sample_name].mean()
            # FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. 
            # Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
            region_dict[idx] = gene.groupby('Bins', observed=True)[self.sample_name].mean()
        
        region_df = pd.concat(region_dict, axis=0)
        region_df.reset_index(inplace=True)
        region_df.rename(columns={'level_0': 'Gene'}, inplace=True)

        return region_df
    

    def process_utr5(self):

        print('Processing 5-UTR.', flush=True)

        # filter utr5 length
        fit_gene = self.gene_table.loc[self.gene_table['utr5_length'] > self.utr5_bin * 2, ]
        fit_gene_num = fit_gene.shape[0]

        if fit_gene_num == 0:
            # create the 0 dictionary with the same dimension
            self.utr5_dict = OrderedDict()
            for idx, gene in self.high_rpf.groupby('name'):
                self.utr5_dict[idx] = pd.DataFrame(np.zeros((self.utr5_bin, self.sample_num)), index=range(self.utr5_bin), columns=self.sample_name)

            # create the 0 dataframe for the output
            self.utr5_df = pd.concat(self.utr5_dict, axis=0).reset_index()
            self.utr5_df.rename(columns={'level_0': 'Gene', 'level_1': 'Bins'}, inplace=True)
            self.utr5_mean = pd.DataFrame(np.zeros((self.utr5_bin, self.sample_num)), index=range(self.utr5_bin), columns=self.sample_name)
            self.utr5_mean.index.name = 'Bins'

            return
        else:
            print('The number of gene with 5-UTR %d.' % fit_gene_num, flush=True)

        # calculate the mean of the gene rpfs in 5-UTR
        utr5_rpm = self.high_rpf.loc[self.high_rpf['region'] == '5utr', ]
        utr5_rpm = utr5_rpm.loc[utr5_rpm['name'].isin(fit_gene['transcript_id']), ]
        self.utr5_df = self.adjust_dimension(rpm = utr5_rpm, bins = self.utr5_bin)

        # check the outlier
        if self.outlier:
            print('Filter the outliers.', flush=True)
            self.utr5_mean = self.filter_outliers(self, self.utr5_df)
        else:
            print('No outliers filter.', flush=True)
            self.utr5_mean = self.utr5_df.groupby('Bins')[self.sample_name].sum() / self.gene_num


    def process_cds(self):

        print('Processing CDS.', flush=True)

        # filter cds length
        fit_gene = self.gene_table.loc[self.gene_table['cds_length'] > self.cds_bin * 2, ]
        fit_gene_num = fit_gene.shape[0]

        if fit_gene_num == 0:
            # create the 0 dictionary with the same dimension
            self.cds_dict = OrderedDict()
            for idx, gene in self.high_rpf.groupby('name'):
                self.cds_dict[idx] = pd.DataFrame(np.zeros((self.cds_bin, self.sample_num)), index=range(self.cds_bin), columns=self.sample_name)

            # create the 0 dataframe for the output
            self.cds_df = pd.concat(self.cds_dict, axis=0).reset_index()
            self.cds_df.rename(columns={'level_0': 'Gene', 'level_1': 'Bins'}, inplace=True)
            self.cds_mean = pd.DataFrame(np.zeros((self.cds_bin, self.sample_num)), index=range(self.cds_bin), columns=self.sample_name)
            self.cds_mean.index.name = 'Bins'

            return
        else:
            print('The number of gene with CDS %d.' % fit_gene_num, flush=True)

        # calculate the mean of the gene rpfs in CDS
        cds_rpm = self.high_rpf.loc[self.high_rpf['region'] == 'cds', ]
        cds_rpm = cds_rpm.loc[cds_rpm['name'].isin(fit_gene['transcript_id']), ]
        self.cds_df = self.adjust_dimension(rpm = cds_rpm, bins = self.cds_bin)
        self.cds_df.rename(columns={'level_0': 'Gene', 'level_1': 'Bins'}, inplace=True)

        # check the outlier
        if self.outlier:
            print('Filter the outliers.', flush=True)
            self.cds_mean = self.filter_outliers(self, self.cds_df)
        else:
            print('No outliers filter.', flush=True)
            self.cds_mean = self.cds_df.groupby('Bins')[self.sample_name].sum() / self.gene_num
        

    def process_utr3(self):

        print('Processing 3-UTR.', flush=True)

        # filter utr3 length
        fit_gene = self.gene_table.loc[self.gene_table['utr3_length'] > self.utr3_bin * 2, ]
        fit_gene_num = fit_gene.shape[0]

        if fit_gene_num == 0:
            # create the 0 dictionary with the same dimension
            self.utr3_dict = OrderedDict()
            for idx, gene in self.high_rpf.groupby('name'):
                self.utr3_dict[idx] = pd.DataFrame(np.zeros((self.utr3_bin, self.sample_num)), index=range(self.utr3_bin), columns=self.sample_name)

            # create the 0 dataframe for the output
            self.utr3_df = pd.concat(self.utr3_dict, axis=0).reset_index()
            self.utr3_df.rename(columns={'level_0': 'Gene', 'level_1': 'Bins'}, inplace=True)
            self.utr3_mean = pd.DataFrame(np.zeros((self.utr3_bin, self.sample_num)), index=range(self.utr3_bin), columns=self.sample_name)
            self.utr3_mean.index.name = 'Bins'

            return
        else:
            print('The number of gene with 3-UTR %d.' % fit_gene_num, flush=True)

        # calculate the mean of the gene rpfs in 3-UTR
        utr3_rpm = self.high_rpf.loc[self.high_rpf['region'] == '3utr', ]
        utr3_rpm = utr3_rpm.loc[utr3_rpm['name'].isin(fit_gene['transcript_id']), ]
        self.utr3_df = self.adjust_dimension(rpm = utr3_rpm, bins = self.utr3_bin)

        # check the outlier
        if self.outlier:
            print('Filter the outliers.', flush=True)
            self.utr3_mean = self.filter_outliers(self, self.utr3_df)
        else:
            print('No outliers filter.', flush=True)
            self.utr3_mean = self.utr3_df.groupby('Bins')[self.sample_name].sum() / self.gene_num


    def output_meta_gene(self):
        # output the mean coverage of
        self.utr5_mean.insert(0, 'Region', '5-UTR')
        self.cds_mean.insert(0, 'Region', 'CDS')
        self.utr3_mean.insert(0, 'Region', '3-UTR')

        self.utr5_mean.reset_index(inplace=True)
        self.cds_mean.reset_index(inplace=True)
        self.utr3_mean.reset_index(inplace=True)

        self.utr5_mean['Bins'] = list(range(self.utr5_bin))
        self.cds_mean['Bins'] = list(range(self.utr5_bin, self.utr5_bin + self.cds_bin))
        self.utr3_mean['Bins'] = list(range(self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin))

        utr5_mean_melt = pd.melt(self.utr5_mean, id_vars=['Bins', 'Region'], var_name='Sample', value_name='Density')
        cds_mean_melt = pd.melt(self.cds_mean, id_vars=['Bins', 'Region'], var_name='Sample', value_name='Density')
        utr3_mean_melt = pd.melt(self.utr3_mean, id_vars=['Bins', 'Region'], var_name='Sample', value_name='Density')

        mean_coverage = pd.concat([utr5_mean_melt, cds_mean_melt, utr3_mean_melt], ignore_index=True)
        mean_coverage.sort_values(by=['Sample', 'Bins'], inplace=True)
        mean_coverage = mean_coverage[['Sample', 'Region', 'Bins', 'Density']]
        mean_coverage.to_csv(self.output + '_utr5_cds_utr3_mean_coverage.txt', sep='\t', index=False)
        
    def draw_meta_gene_line(self):
        merge_coverage = pd.concat([self.utr5_mean, self.cds_mean, self.utr3_mean], ignore_index=True)
        merge_coverage.index = merge_coverage.index + 1

        for sp in self.sample_name:
            out_pdf = sp + "_coverage_line_plot.pdf"
            out_png = sp + "_coverage_line_plot.png"
            print('Draw {sp} line plot.'.format(sp=sp), flush=True)
            matplotlib.use('AGG')
            fig = plt.figure(figsize=(8, 5), dpi=300)
            plt.plot(merge_coverage[sp], label=sp)
            plt.legend()
            plt.xticks([1, self.utr5_bin, self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin], ['TSS', 'TIS', 'TTS', 'TES'])
            plt.title("Mean coverage of ({number} genes)".format(number=self.gene_num), fontsize=16)
            if self.norm:
                plt.ylabel('Mean coverage (RPM)', fontsize=14)
            else:
                plt.ylabel('Mean coverage (RPFs)', fontsize=14)
            plt.xlabel('position of genes', fontsize=14)
            fig.tight_layout()
            # plt.show()

            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()

    def draw_meta_gene_heat(self):
        
        utr5_df = self.utr5_df.reset_index().rename(columns={'level_0': 'Gene'})
        cds_df = self.cds_df.reset_index().rename(columns={'level_0': 'Gene'})
        utr3_df = self.utr3_df.reset_index().rename(columns={'level_0': 'Gene'})

        for sp in self.sample_name:
            print('Draw {sp} heatmap.'.format(sp=sp), flush=True)
            utr5_sp = pd.pivot_table(utr5_df.loc[:, ['Gene', "Bins", sp]],
                                     index='Gene',
                                     values=sp,
                                     columns='Bins')
            cds_sp = pd.pivot_table(cds_df.loc[:, ['Gene', "Bins", sp]],
                                     index='Gene',
                                     values=sp,
                                     columns='Bins')
            utr3_sp = pd.pivot_table(utr3_df.loc[:, ['Gene', "Bins", sp]],
                                     index='Gene',
                                     values=sp,
                                     columns='Bins')
            coverage_sp = pd.concat([utr5_sp, cds_sp, utr3_sp], axis=1, ignore_index=True).fillna(0)
            coverage_sp.index.name = 'Gene'
            coverage_sp = coverage_sp.reindex(coverage_sp.mean(axis=1).sort_values(ascending=False).index, axis=0)

            gene_bins = '{utr5}_{cds}_{utr3}'.format(utr5=str(self.utr5_bin),
                                                     cds=str(self.cds_bin),
                                                     utr3=str(self.utr3_bin))
            file_name = "{sp}_{bins}_coverage.txt".format(sp=sp, bins=gene_bins)
            coverage_sp.to_csv(file_name, sep='\t', index=True)

            # out_pdf = sp + "_heat_plot.pdf"
            out_png = sp + "_" + gene_bins + "_heat_plot.png"

            matplotlib.use('AGG')
            fig = plt.figure(figsize=(8, 4 + self.gene_num / 500), dpi=300)
            # sns.clustermap(np.log2(coverage_sp + 1), col_cluster=False, cmap="YlGnBu")
            sns.heatmap(np.log2(coverage_sp + 1), cmap="YlGnBu")
            plt.xticks([1, self.utr5_bin, self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin], ['TSS', 'TIS', 'TTS', 'TES'])
            plt.title("Mean coverage of ({number} genes)".format(number=self.gene_num))
            fig.tight_layout()
            # plt.show()

            # fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()
