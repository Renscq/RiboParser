#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Odd_Ratio.py
# @Author  : Rensc


# from turtle import color
# from scipy.stats import fisher_exact
from joblib import Parallel, delayed
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from . import RPFs


class OddRatio(object):
    def __init__(self, args):
        # input and output files
        self.rpf_file = args.rpf
        self.output_prefix = args.output

        # opts for odd ratio
        self.gene = args.list
        self.site = args.site
        self.frame = args.frame
        self.tis = args.tis
        self.tts = args.tts
        
        self.normal = args.normal
        self.flag = args.fdr
        self.pvalue = args.value
        self.fillna = args.zero

        self.thread = args.thread

        # samples
        self.control = args.control.split(',')
        self.treat = args.treat.split(',')
        self.rpf_num = args.min

        self.sample_name = self.control + self.treat
        self.sample_num = len(self.sample_name)

        # for the high expression gene filter
        self.merged_rpf = None
        self.raw_rpf = None
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None
        self.dim2_df = None
        self.merge_odd_ratio = None

        # summary stalling codon
        self.scale = args.scale
        self.codon_odd_ratio = None
        self.codon_odd_ratio_num = None
        self.codon_odd_ratio_prop = None

        # codon list
        self.stop = args.stop
        self.codon_dict, self.codon_table = RPFs.codon_table()

        if self.stop:
            del self.codon_dict['TAA']
            del self.codon_dict['TAG']
            del self.codon_dict['TGA']
            
            self.codon_table = self.codon_table[self.codon_table['Abbr'] != 'Stop']

    def read_rpf(self):
        '''
        @Message  : import the rpf density file and convert them to rpm density
        @Input    : self.rpf_file --> RPF density file
                    self.site --> the site of the codon
                    self.frame --> the frame of the codon
                    self.tis --> the start site of the CDS
                    self.tts --> the end site of the CDS
                    self.sample_num --> the number of the samples
                    self.sample_name --> the name of the samples
                    self.gene --> the gene list
                    self.rpf_num --> the minimal number of the RPFs
        @return   : merged_rpf, sample_name, sample_num, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        @Flow     : step 1. import the rpf density file with import_rpf function,
                    step 2. convert the rpf density to rpm density
                    step 3. drop the codon without RPFs coverage
                    step 4. fill the empty codon rpm value with minimal rpm value (>0),
        '''

        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file,
                                      sites=self.site, frame=self.frame,
                                      tis=self.tis, tts=self.tts,
                                      sample_num=self.sample_num, sample_name=self.sample_name,
                                      gene=self.gene, rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0].copy()
        self.sample_name = rpf_results[1]
        self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3].copy()
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        # self.high_rpf = rpf_results[7].copy()
        del rpf_results
        
        self.merged_rpf = self.merged_rpf[self.merged_rpf['codon'].isin(self.codon_dict.keys())]

        # convert the rpf density to rpm density
        if self.normal:
            self.merged_rpf.loc[:, self.sample_name] = self.merged_rpf[self.sample_name] * 1e6 / self.total_rpf_num
            # change the float to int with ceil
            # self.merged_rpf.loc[:, self.sample_name] = self.merged_rpf[self.sample_name].apply(np.ceil)
            # self.merged_rpf.loc[:, self.sample_name] = self.merged_rpf[self.sample_name].astype(int)

        self.high_rpf = self.merged_rpf.loc[self.merged_rpf['name'].isin(self.high_gene)]
        
        # drop the codon without RPFs coverage
        print('Ignore the codon without RPFs coverage.', flush=True)
        self.high_rpf = self.high_rpf[self.high_rpf[self.sample_name].sum(axis=1) > 0]

        # fill emtpy rpf with 0.1 * minimal rpm value (> 0)
        if self.fillna:
            min_rpm = self.high_rpf[self.high_rpf[self.sample_name] !=0 ].min().min()
            self.high_rpf[self.sample_name] += min_rpm * 0.01
        else:
            self.high_rpf = self.high_rpf[self.high_rpf[self.sample_name].min(axis=1) > 0]

        # summary the gene rpm
        self.gene_rpf_sum = self.high_rpf.groupby('name')[self.sample_name].sum()

    def make_two_dimensional_table(self):
        '''
        @Message  : prepare two dimensional table for fisher exactly test.
        @Input    : self.high_rpf --> high expression rpf table
                    self.control --> control sample table
                    self.treat --> treat sample table
        @Return   : self.gene_rpf_sum --> sum of high expression gene table
        @Flow     : 1 --> sum the rpf count of each gene with control and treatment group
                    2 --> calculate the gene delta
                    delta = each gene sum - each codon of this gene
        '''

        self.high_rpf = self.high_rpf.reset_index(drop=True)

        if len(self.control) > 1:
            self.high_rpf.loc[:, 'control'] = np.array(self.high_rpf[self.control].sum(axis=1))
            self.gene_rpf_sum.loc[:, 'control_sum'] = np.array(self.gene_rpf_sum[self.control].sum(axis=1))
        else:
            self.high_rpf.loc[:, 'control'] = np.array(self.high_rpf[self.control])
            self.gene_rpf_sum.loc[:, 'control_sum'] = np.array(self.gene_rpf_sum[self.control])

        if len(self.treat) > 1:
            self.high_rpf.loc[:, 'treat'] = np.array(self.high_rpf[self.treat].sum(axis=1))
            self.gene_rpf_sum.loc[:, 'treat_sum'] = np.array(self.gene_rpf_sum[self.treat].sum(axis=1))
        else:
            self.high_rpf.loc[:, 'treat'] = np.array(self.high_rpf[self.treat])
            self.gene_rpf_sum.loc[:, 'treat_sum'] = np.array(self.gene_rpf_sum[self.treat])

        # rpf_sum = np.array([self.total_rpf_num[self.control].sum(), self.total_rpf_num[self.treat].sum()])
        # rpf_sum = np.array([self.high_rpf['control'].sum(), self.high_rpf['treat'].sum()])
        # self.high_rpf.loc[:, ['control_sum', 'treat_sum']] = rpf_sum - np.array(self.high_rpf[['control', 'treat']])

        self.gene_rpf_sum = self.gene_rpf_sum[['control_sum', 'treat_sum']].reset_index()
        self.high_rpf = pd.merge(self.high_rpf, self.gene_rpf_sum, on='name', how='left')

    def calc_odd_ratio_v1(self):
        '''
        @Message  : fisher exactly test with fisher package.
                    package fisher were developed with python version < 3.9
        @Input    : self.dim2_df --> two dimension table for fisher-exactly test
        @Return   : control --> rpm of control
                    treat --> rpm of treat
                    control_delta --> delta rpm of control, 
                    treat_delta --> delta rpm of treat samples
                    self.high_rpf --> contain the high rpf and odd-ratio
                    bhfdr --> calculate bhfdr with multipletests
        @Flow     : step1 --> run
        '''
        
        from fisher import pvalue_npy

        self.dim2_df = self.high_rpf[['control', 'treat', 'control_sum', 'treat_sum']]

        # set the table
        control = self.dim2_df['control'].to_numpy(dtype='float32')
        treat = self.dim2_df['treat'].to_numpy(dtype='float32')
        control_delta = self.dim2_df['control_sum'].to_numpy(dtype='float32')
        treat_delta = self.dim2_df['treat_sum'].to_numpy(dtype='float32')

        self.high_rpf.loc[:, ['control', 'treat']] = self.dim2_df.loc[:, ['control', 'treat']].copy()

        # calculate the odd ratio
        print('\n\n' + 'Ignore codon with the none RPFs coverage.', flush=True)
        np.seterr(divide='ignore', invalid='ignore')
        self.high_rpf.loc[:, 'odd2'] = (treat * control_delta) / (control * treat_delta)

        # calculate the pvalue
        p_left, p_right, p_two_tail = pvalue_npy(control, treat_delta, treat, control_delta)
        bhfdr = multipletests(p_two_tail, alpha=0.05, method='fdr_bh')

    def calc_odd_ratio(self):
        '''
        @Message  : calculate the odd ratio and pvalue with chi2-test
        @Input    : self.high_rpf --> high expression rpf table
        @Return   : self.high_rpf --> contain the high rpf and odd-ratio
        @Flow     : step1 --> calculate the delta of each gene
                    step2 --> calculate the odd ratio
        '''
        
        control = self.high_rpf['control'].to_numpy(dtype='float32')
        treat = self.high_rpf['treat'].to_numpy(dtype='float32')
        control_delta = self.high_rpf['control_sum'].to_numpy(dtype='float32') - control
        treat_delta = self.high_rpf['treat_sum'].to_numpy(dtype='float32') - treat
        
        np.seterr(divide='ignore', invalid='ignore')
        self.high_rpf.loc[:, 'odd'] = (treat / treat_delta) / (control / control_delta)

    def calc_chi2_test(self):
        '''
        @Message  : calculate the odd ratio and pvalue with chi2-test
        @Input    : self.high_rpf --> high expression rpf table
        @Return   : self.high_rpf --> contain the high rpf and odd-ratio and pvalue
        @Flow     : step1 --> calculate the odd ratio and pvalue
                    step2 --> merge all results into one table

        '''
        import statsmodels.api as sm

        # calculate the pvalue with chi2-test
        statistic_list = []
        pvalue_list = []

        for i in range(len(self.high_rpf)):
            if i % 5000 == 0:
                print('row number: ' + str(i), flush=True)

            n = self.high_rpf['control'][i]
            m = self.high_rpf['treat'][i]
            N = self.high_rpf['control_sum'][i]
            M = self.high_rpf['treat_sum'][i]

            # calculate the odd ratio with statmodels

            # from scipy.stats import fisher_exact
            # observed = np.array([[n, N - n], [m, M - m]])
            # statistic, p_val = fisher_exact(observed)

            # ATTENTION： fisher-exact test does not fit the float values !!!
            table = [[n, m], [N - n, M - m]]
            result = sm.stats.Table2x2(table)

            # table = sm.stats.Table([[n, m], [N - n, M - m]])
            # result = table.test_nominal_association()

            statistic_list.append(result.oddsratio)
            pvalue_list.append(result.oddsratio_pvalue())

        print('row number: ' + str(i), flush=True)

        # correct the pvalue with Benjamin-Hochberg
        bhfdr = multipletests(pvalue_list, alpha=0.05, method='fdr_bh')

        # merge all results
        self.high_rpf.loc[:, 'statistic'] = statistic_list
        self.high_rpf.loc[:, 'pvalue'] = pvalue_list
        self.high_rpf.loc[:, 'bhfdr'] = bhfdr[1]
        self.high_rpf.loc[:, 'flag'] = bhfdr[0]

    def calc_chi2_test2(self):
        '''
        @Message  : calculate the chi2-test odd ratio and pvalue with multiple progress
        @Input    : self.high_rpf --> high expression rpf table
        @Return   : self.high_rpf --> contain the high rpf and odd-ratio and pvalue
        @Flow     : step1 --> define the function to calculate the pvalue with multiple progress
                    step2 --> calculate the odd ratio and pvalue
                    step3 --> correct the pvalue with Benjamin-Hochberg
                    step4 --> merge all results into one table

        '''
        import statsmodels.api as sm

        # define the function to calculate the pvalue with multiple progress
        def calc_pvalue(n, m, N, M):
            table = [[n, m], [N - n, M - m]]
            result = sm.stats.Table2x2(table)
            return result.oddsratio_pvalue()

        # calculate the pvalue with multiple progress
        pvalue_list = Parallel(n_jobs=self.thread)(delayed(calc_pvalue)(n, m, N, M) for n, m, N, M in zip(self.high_rpf['control'],
                                                                                                    self.high_rpf['treat'],
                                                                                                    self.high_rpf['control_sum'],
                                                                                                    self.high_rpf['treat_sum']))
        
        # correct the pvalue with Benjamin-Hochberg
        bhfdr = multipletests(pvalue_list, alpha=0.05, method='fdr_bh')

        # merge all results
        self.high_rpf.loc[:, 'pvalue'] = pvalue_list
        self.high_rpf.loc[:, 'bhfdr'] = bhfdr[1]
        self.high_rpf.loc[:, 'flag'] = bhfdr[0]

    def output_odd_ratio(self):
        '''
        output the odd ratio table
        '''
        out_odd_ratio_txt = self.output_prefix + '_codon_odd_ratio.txt'
        self.merge_odd_ratio = self.high_rpf[(self.high_rpf['from_tis'] >= self.tis) &
                                             (self.high_rpf['from_tts'] <= -self.tts) &
                                             (self.high_rpf[self.flag] < self.pvalue)]
        self.merge_odd_ratio.to_csv(out_odd_ratio_txt, sep='\t', index=False)

    @staticmethod
    def scale_method(scale, oddratio):
        if scale == 'minmax':
            # relative_oddratio = (oddratio - oddratio.min()) / (oddratio.max() - oddratio.min())
            relative_oddratio = oddratio / oddratio.max()
        elif scale == 'zscore':
            relative_oddratio = (oddratio - oddratio.mean()) / oddratio.std()
        else:
            print('Unknown scale method.', flush=True)
            sys.exit()
        return relative_oddratio
    
    def summarize_odd_ratio(self):
        '''
        @Message  : summarize the odd ratio
        @Input    : self.high_rpf --> high expression rpf table
        @Return   : count_aa_odd --> the count of codon odd ratio
                    per_aa_odd --> the proportion of codon odd ratio
        @Flow     : step1 --> filter the significant codon odd ratio
                    step2 --> calculate the count of significant codon odd ratio
                    step3 --> calculate the proportion of significant codon odd ratio

        '''
        # output filename
        codon_odd_file = self.output_prefix + '_sum_codon_odd_ratio.txt'

        codon_odd_num_column = ['Codon', 'AA', 'Abbr', 
                                 'codon_sum', 'control_codon_sum', 'treat_codon_sum', 
                                 'control_mean_odd_ratio', 'treat_mean_odd_ratio',
                                 'control_number', 'treat_number']

        # filter the significant codon odd ratio
        codon_odd_sum_all = self.high_rpf.groupby('codon')['control'].count()
        codon_odd_sum_all_0 = self.high_rpf.loc[self.high_rpf['odd'] < 1].groupby('codon')['control'].count()
        codon_odd_sum_all_1 = self.high_rpf.loc[self.high_rpf['odd'] > 1].groupby('codon')['treat'].count()
        
        # summary the odd ratio of each codon
        codon_odd_mean_all_0 = self.high_rpf.loc[self.high_rpf['odd'] < 1].groupby('codon')['control'].mean()
        codon_odd_mean_all_1 = self.high_rpf.loc[self.high_rpf['odd'] > 1].groupby('codon')['treat'].mean()

        sig_high_rpf = self.high_rpf.loc[self.high_rpf[self.flag] < self.pvalue, :]
        codon_odd_num_0 = sig_high_rpf.loc[sig_high_rpf['odd'] < 1].groupby('codon')['control'].count()
        codon_odd_num_1 = sig_high_rpf.loc[sig_high_rpf['odd'] > 1].groupby('codon')['treat'].count()

        # sig oddratio number
        codon_odd_ratio = pd.concat([self.codon_table,
                                    codon_odd_sum_all, codon_odd_sum_all_0, codon_odd_sum_all_1, 
                                    codon_odd_mean_all_0, codon_odd_mean_all_1,
                                    codon_odd_num_0, codon_odd_num_1], 
                                    join='inner', axis=1)
        codon_odd_ratio.reset_index(inplace=True)
        codon_odd_ratio.columns = codon_odd_num_column

        # filter the stop codon
        if self.stop:
            codon_odd_ratio = codon_odd_ratio.loc[~codon_odd_ratio['Codon'].isin(["TGA", "TAG", "TAA"]), :]
            # codon_odd_ratio_prop = codon_odd_ratio_prop.loc[~codon_odd_ratio_prop['Codon'].isin(["TGA", "TAG", "TAA"]), :]
        else:
            pass

        # sig oddratio proportion
        codon_odd_ratio_prop = codon_odd_ratio.loc[:, ['control_number', 'treat_number']].div(codon_odd_ratio['codon_sum'], axis=0) * 100
        codon_odd_ratio_prop.columns = ['control_proportion', 'treat_proportion']
        
        # calculate the relative odd ratio
        codon_odd_ratio_num_rel = self.scale_method(self.scale, codon_odd_ratio.loc[:, ['control_number', 'treat_number']])
        codon_odd_ratio_prop_rel = self.scale_method(self.scale, codon_odd_ratio_prop)
        codon_odd_ratio_num_rel.columns = ['control_number_relative', 'treat_number_relative']
        codon_odd_ratio_prop_rel.columns = ['control_proportion_relative', 'treat_proportion_relative']

        self.codon_odd_ratio = pd.concat([codon_odd_ratio, codon_odd_ratio_prop, codon_odd_ratio_num_rel, codon_odd_ratio_prop_rel], axis=1)
        self.codon_odd_ratio.sort_values(by=['Abbr', 'Codon'], inplace=True)

        self.codon_odd_ratio['control_group'] = ",".join(self.control)
        self.codon_odd_ratio['treat_group'] = ",".join(self.treat)

        # melt the table
        self.codon_odd_ratio_num = pd.melt(self.codon_odd_ratio, 
                                    id_vars=['Codon', 'AA', 'Abbr'], value_vars=['control_number', 'treat_number'],
                                    var_name='groups', value_name='Number')
        self.codon_odd_ratio_num['Codon'] = self.codon_odd_ratio_num['Codon'] + '(' + self.codon_odd_ratio_num['Abbr'] + ')'

        self.codon_odd_ratio_prop = pd.melt(self.codon_odd_ratio, 
                                    id_vars=['Codon', 'AA', 'Abbr'], value_vars=['control_proportion', 'treat_proportion'],
                                    var_name='groups', value_name='Proportion')
        self.codon_odd_ratio_prop['Codon'] = self.codon_odd_ratio_prop['Codon'] + '(' + self.codon_odd_ratio_prop['Abbr'] + ')'

        # output the odd ratio table
        self.codon_odd_ratio.to_csv(codon_odd_file, sep='\t', index=False)


    def draw_odd_ratio_bar(self):
        '''
        @Message  : draw the odd ratio barplot
        @Input    : self.codon_odd_ratio_num --> the count of codon odd ratio
                    self.codon_odd_ratio_prop --> the proportion of codon odd ratio
        @Return   : odd_ratio_barplot.pdf

        '''

        out_pdf = self.output_prefix + "_odd_barplot.pdf"
        out_png = self.output_prefix + "_odd_barplot.png"
        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15, 12), dpi=300)
        sns.barplot(data=self.codon_odd_ratio_num, x='Codon', y='Number', hue='groups', ax=ax[0])
        sns.barplot(data=self.codon_odd_ratio_prop, x='Codon', y='Proportion', hue='groups', ax=ax[1])

        # add the gene annotation here, need to import the txt file
        ax[0].set_ylabel('Number')
        ax[1].set_ylabel('Proportion (%)')
        ax[1].set_xlabel('Codon')

        plt.sca(ax[0])
        plt.xticks(rotation=90)
        plt.sca(ax[1])
        plt.xticks(rotation=90)

        plt.suptitle('Significant differential pausing codon')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_odd_ratio_line(self):
        '''
        @Message  : draw the odd ratio lineplot
        @Input    : self.codon_odd_ratio_num --> the count of codon odd ratio
                    self.codon_odd_ratio_prop --> the proportion of codon odd ratio
        @Return   : odd_ratio_lineplot.pdf

        '''
        out_pdf = self.output_prefix + "_odd_lineplot.pdf"
        out_png = self.output_prefix + "_odd_lineplot.png"

        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 6), dpi=300)

        flag = 1
        for aa in self.codon_odd_ratio_num['AA'].unique().tolist():
            tmp = self.codon_odd_ratio_num.loc[self.codon_odd_ratio_num['AA']==aa, :]
            if flag == 1:
                sns.lineplot(data=tmp, x='Codon', y='Number', linewidth=1.5, hue='groups', style='groups',
                              markers=['d','o'], legend = True, ax=ax[0])
                flag += 1
            else:
                sns.lineplot(data=tmp, x='Codon', y='Number', linewidth=1.5, hue='groups', style='groups',
                              markers=['d', 'o'], legend = False, ax=ax[0])

        flag = 1
        for aa in self.codon_odd_ratio_prop['AA'].unique().tolist():
            tmp = self.codon_odd_ratio_prop.loc[self.codon_odd_ratio_prop['AA']==aa, :]
            if flag == 1:
                sns.lineplot(data=tmp, x='Codon', y='Proportion', linewidth=1.5, hue='groups', style='groups',
                             markers=['d', 'o'], legend = True, ax=ax[1])
                flag += 1
            else:
                sns.lineplot(data=tmp, x='Codon', y='Proportion', linewidth=1.5, hue='groups', style='groups',
                             markers=['d', 'o'], legend = False, ax=ax[1])

        # add the gene annotation here, need to import the txt file
        ax[0].legend(loc=0, ncol=2, fontsize=6)
        ax[1].legend(loc=0, ncol=2, fontsize=6)

        ax[0].set_ylabel('Number')
        ax[1].set_ylabel('Proportion (%)')
        ax[1].set_xlabel('Codon')

        plt.sca(ax[0])
        plt.xticks(rotation=90)
        plt.sca(ax[1])
        plt.xticks(rotation=90)

        plt.suptitle('Number of sig. diff. pausing codon')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)

    def draw_odd_ratio_scatter(self):
        '''
        @Message  : draw the odd ratio scatterplot
        @Input    : self.codon_odd_ratio_num --> the count of codon odd ratio
                    self.codon_odd_ratio_prop --> the proportion of codon odd ratio
        @Return   : odd_ratio_scatterplot.pdf

        '''

        out_pdf = self.output_prefix + "_odd_scatter.pdf"
        out_png = self.output_prefix + "_odd_scatter.png"

        # count_df = self.codon_odd_ratio_num.pivot_table(index=['Codon', 'AA', 'Abbr'], columns='groups', values='Number').reset_index()
        # per_df = self.codon_odd_ratio_prop.pivot_table(index=['Codon', 'AA', 'Abbr'], columns='groups', values='Proportion').reset_index()

        # draw the odd ratio number of each codon
        matplotlib.use('AGG')

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 6.5), dpi=300)
        
        count_plot = sns.scatterplot(data=self.codon_odd_ratio, x='control_number', y='treat_number', 
                                     s = 80, hue='AA', legend = 'brief', ax=ax[0])
        count_plot.legend(loc=8, bbox_to_anchor=(0.5, -0.5), ncol=7, fontsize=6)

        per_plot = sns.scatterplot(data=self.codon_odd_ratio, x='control_proportion', y='treat_proportion', 
                                   s = 80, hue='AA', legend = 'brief', ax=ax[1])
        per_plot.legend(loc=8, bbox_to_anchor=(0.5, -0.5), ncol=7, fontsize=6)

        # add the gene annotation here, need to import the txt file
        ax[0].set_title('Number')
        ax[1].set_title('Proportion (%)')

        plt.suptitle('Number of sig. diff. pausing codon')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
