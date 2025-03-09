#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Coverage.py


from concurrent.futures import ThreadPoolExecutor
from scipy.stats import zscore
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
        self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        self.gene_num = len(self.high_gene)
        self.high_rpf = rpf_results[7].copy()
        del rpf_results

        if self.norm:
            self.high_rpf[self.sample_name] = self.high_rpf[self.sample_name].astype(float)

            self.high_rpf[self.sample_name] = self.high_rpf[self.sample_name].div(self.total_rpf_num) * 1e6


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
        '''
        @Message  : function for fill the outliers of the gene coverage.
        @Input    : param --> gene rpf density dataframe
        @Return   : output --> gene rpf density dataframe
        @Flow     : step1 --> calculate the iqr value
                    step2 --> calculate the lower bound and upper bound
                    step3 --> fill the outliers with the mean value
        '''
        
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
        gene_coverage = gene_cover_df_n.groupby('Bins', observed=False)[self.sample_name].sum() / self.gene_num
        
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
        '''
        @Message  : function for process the 5-UTR coverage
        @Input    : param --> gene rpf density dataframe of 5-UTR
        @Return   : output --> the mean of the gene rpfs dataframe of 5-UTR
        @Flow     : step1 --> filter the utr5 length
                    step2 --> calculate the mean of the gene rpfs in 5-UTR
                    step3 --> check the outlier
        '''
        
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
            self.utr5_mean = self.utr5_df.groupby('Bins', observed=False)[self.sample_name].sum() / self.gene_num


    def process_cds(self):
        '''
        @Message  : function for process the CDS coverage
        @Input    : param --> gene rpf density dataframe of CDS
        @Return   : output --> the mean of the gene rpfs dataframe of CDS
        @Flow     : step1 --> filter the CDS length
                    step2 --> calculate the mean of the gene rpfs in CDS
                    step3 --> check the outlier
        '''

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
            self.cds_mean = self.cds_df.groupby('Bins', observed=False)[self.sample_name].sum() / self.gene_num.astype(float)
        

    def process_utr3(self):
        '''
        @Message  : function for process the 3-UTR coverage
        @Input    : param --> gene rpf density dataframe of 3-UTR
        @Return   : output --> the mean of the gene rpfs dataframe of 3-UTR
        @Flow     : step1 --> filter the utr3 length
                    step2 --> calculate the mean of the gene rpfs in 3-UTR
                    step3 --> check the outlier
        '''

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
            self.utr3_mean = self.utr3_df.groupby('Bins', observed=False)[self.sample_name].sum() / self.gene_num.astype(float)


    def draw_meta_gene_line(self):
        merge_coverage = pd.concat([self.utr5_mean, self.cds_mean, self.utr3_mean], ignore_index=True)
        merge_coverage.index = merge_coverage.index + 1

        for sp in self.sample_name:
            out_pdf = self.output + "_" + sp + "_coverage_line_plot.pdf"
            out_png = self.output + "_" + sp + "_coverage_line_plot.png"
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
                                     index='Gene', values=sp, columns='Bins', 
                                     observed=False)
            cds_sp = pd.pivot_table(cds_df.loc[:, ['Gene', "Bins", sp]],
                                     index='Gene', values=sp,
                                     columns='Bins', 
                                     observed=False)
            utr3_sp = pd.pivot_table(utr3_df.loc[:, ['Gene', "Bins", sp]],
                                     index='Gene', values=sp, columns='Bins', 
                                     observed=False)
            coverage_sp = pd.concat([utr5_sp, cds_sp, utr3_sp], axis=1, ignore_index=True).fillna(0)
            coverage_sp.index.name = 'Gene'
            coverage_sp = coverage_sp.reindex(coverage_sp.mean(axis=1).sort_values(ascending=False).index, axis=0)

            gene_bins = '{utr5}_{cds}_{utr3}'.format(utr5=str(self.utr5_bin),
                                                     cds=str(self.cds_bin),
                                                     utr3=str(self.utr3_bin))
            file_name = "{prefix}_{sp}_{bins}_coverage.txt".format(prefix=self.output, sp=sp, bins=gene_bins)
            coverage_sp.to_csv(file_name, sep='\t', index=True)

            # out_pdf = sp + "_heat_plot.pdf"
            out_pdf = self.output + "_" + sp + "_" + gene_bins + "_heat_plot.pdf"
            out_png = self.output + "_" + sp + "_" + gene_bins + "_heat_plot.png"

            matplotlib.use('AGG')

            fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [1, 1]})
            sns.heatmap(np.log2(coverage_sp + 1), cmap="Blues", cbar_kws={'label': 'log2(RPM + 1)'}, ax=axes[0])
            axes[0].set_title("Mean coverage of ({number})".format(number=self.gene_num))
            axes[0].set_xticks([1, self.utr5_bin, self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin])
            axes[0].set_xticklabels(['TSS', 'TIS', 'TTS', 'TES'])
            
            sns.heatmap(zscore(coverage_sp, axis = 1), cmap="Blues", cbar_kws={'label': 'Z-score'}, ax=axes[1])
            axes[1].set_title("Z-score of mean coverage of ({number})".format(number=self.gene_num))
            axes[1].set_xticks([1, self.utr5_bin, self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin])
            axes[1].set_xticklabels(['TSS', 'TIS', 'TTS', 'TES'])

            # fig = plt.figure(figsize=(8, 8), dpi=300)
            # # sns.clustermap(np.log2(coverage_sp + 1), row_cluster=False, col_cluster=False, cmap="YlGnBu", z_score=0, cbar_kws={'label': 'log2(RPM)'})
            # sns.heatmap(zscore(np.log2(coverage_sp + 1), axis = 1), cmap="Blues", cbar_kws={'label': 'log2(RPM)'})
            # plt.xticks([1, self.utr5_bin, self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin], ['TSS', 'TIS', 'TTS', 'TES'])
            # plt.title("Mean coverage of ({number} genes)".format(number=self.gene_num))

            fig.tight_layout()
            # plt.show()

            # fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()

    def draw_meta_gene_bar(self):

        utr5_df = self.utr5_df.reset_index(drop=True).rename(columns={'level_0': 'Gene'})
        cds_df = self.cds_df.reset_index(drop=True).rename(columns={'level_0': 'Gene'})
        utr3_df = self.utr3_df.reset_index(drop=True).rename(columns={'level_0': 'Gene'})
        
        # replace the all value > 0 with 1
        utr5_df.loc[:, self.sample_name] = utr5_df[self.sample_name].map(lambda x: 1 if x > 0 else 0)
        cds_df.loc[:, self.sample_name] = cds_df[self.sample_name].map(lambda x: 1 if x > 0 else 0)
        utr3_df.loc[:, self.sample_name] = utr3_df[self.sample_name].map(lambda x: 1 if x > 0 else 0)

        # group the gene by the bins and sum the value
        utr5_sp = utr5_df.groupby('Bins', observed=False)[self.sample_name].sum().div(self.gene_num) * 100
        cds_sp = cds_df.groupby('Bins', observed=False)[self.sample_name].sum().div(self.gene_num) * 100
        utr3_sp = utr3_df.groupby('Bins', observed=False)[self.sample_name].sum().div(self.gene_num) * 100

        utr5_sp.insert(0, 'Region', '5-UTR')
        cds_sp.insert(0, 'Region', 'CDS')
        utr3_sp.insert(0, 'Region', '3-UTR')

        merged_sp = pd.concat([utr5_sp, cds_sp, utr3_sp], axis=0, ignore_index=True).reset_index(names=['Bins'])
        merged_sp['Bins'] += 1

        merged_sp_melt = pd.melt(merged_sp, id_vars=['Bins', 'Region'], var_name='Sample', value_name='Percentage')

        merged_sp_melt.sort_values(by=['Sample', 'Bins'], inplace=True)
        merged_sp_melt = merged_sp_melt[['Sample', 'Region', 'Bins', 'Percentage']]

        # output the coverage percentage of the gene rpfs
        merged_sp_melt.to_csv(self.output + '_utr5_cds_utr3_mean_coverage_percentage.txt', sep='\t', index=False)

        # draw the cds_sp bar plot, x= bins, y= cds_sp
        for sp in self.sample_name:
            print('Draw percentage barplot of {sp}.'.format(sp=sp), flush=True)

            out_pdf = self.output + "_" + sp + "_coverage_bar_plot.pdf"
            out_png = self.output + "_" + sp + "_coverage_bar_plot.png"

            matplotlib.use('AGG')

            fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), 
                                     gridspec_kw={'width_ratios': [self.utr5_bin, self.cds_bin, self.utr3_bin]})
            
            # barplot of the utr5
            utr5_sp = merged_sp.loc[merged_sp['Region'] == '5-UTR', ]
            axes[0].bar(utr5_sp['Bins'], utr5_sp[sp], color='#23a9f2')
            # set the xlabel and ylabel
            axes[0].set_ylabel('Percentage (%)')
            axes[0].set_xlabel("5' UTR")
            axes[0].set_ylim(0, 100)
            axes[0].set_xlim(0, self.utr5_bin + 1)

            # barplot of the cds
            cds_sp = merged_sp.loc[merged_sp['Region'] == 'CDS', ]
            axes[1].bar(cds_sp['Bins'], cds_sp[sp], color='#23a9f2')
            axes[1].set_ylabel('Percentage (%)')
            axes[1].set_xlabel("CDS")
            axes[1].set_ylim(0, 100)
            axes[1].set_xlim(self.utr5_bin, self.utr5_bin + self.cds_bin + 1)
            
            # barplot of the utr3
            utr3_sp = merged_sp.loc[merged_sp['Region'] == '3-UTR', ]
            axes[2].bar(utr3_sp['Bins'], utr3_sp[sp], color='#23a9f2')
            axes[2].set_ylabel('Percentage (%)')
            axes[2].set_xlabel("3' UTR")
            axes[2].set_ylim(0, 100)
            axes[2].set_xlim(self.utr5_bin + self.cds_bin, self.utr5_bin + self.cds_bin + self.utr3_bin + 1)

            fig.suptitle("Coverage of ({number} genes)".format(number=self.gene_num), fontsize=16)
            fig.tight_layout()
            # plt.show()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()


    def output_meta_gene(self):
        '''
        @Message  : output the mean coverage of the gene rpfs
        @Input    : self.utr5_mean, self.cds_mean, self.utr3_mean --> the mean of the gene rpfs dataframe
        @Return   : self.utr5_mean, self.cds_mean, self.utr3_mean --> the mean of the gene rpfs dataframe
        @Flow     : step1 --> output the mean coverage of the gene rpfs
        '''

        # output the mean coverage of
        self.utr5_mean.insert(0, 'Region', '5-UTR')
        self.cds_mean.insert(0, 'Region', 'CDS')
        self.utr3_mean.insert(0, 'Region', '3-UTR')

        mean_coverage = pd.concat([self.utr5_mean, self.cds_mean, self.utr3_mean], ignore_index=True).reset_index(names='Bins')
        mean_coverage['Bins'] += 1
        mean_coverage_melt = pd.melt(mean_coverage, id_vars=['Bins', 'Region'], var_name='Sample', value_name='Density')

        mean_coverage_melt.sort_values(by=['Sample', 'Bins'], inplace=True)
        mean_coverage_melt = mean_coverage_melt[['Sample', 'Region', 'Bins', 'Density']]
        mean_coverage_melt.to_csv(self.output + '_utr5_cds_utr3_mean_coverage_density.txt', sep='\t', index=False)
        