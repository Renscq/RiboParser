#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : Coefficient_of_Variation.py


from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pandas as pd
import seaborn as sns
from . import RPFs


class CoV(object):
    def __init__(self, args):
        self.rpf = args.rpf

        # self.anno = args.anno
        self.output_prefix = args.output

        # parameters for the trimmed CDS region
        self.site = 'P'
        self.tis = args.tis
        self.tts = args.tts
        self.rpf_num = args.min
        self.frame = args.frame
        self.length = None
        self.norm = args.normal

        # parameters for groups file parser
        self.group_file = args.group
        self.group = None

        # parameters for rpf file parser
        self.gene = args.list
        self.all_rpf = pd.DataFrame()
        self.sample_name = None
        self.sample_num = None
        self.col_index = None

        # quantity of cds
        self.CoV_dict = OrderedDict()

        self.cds_rpf = None
        self.cds_rpm = None
        self.cds_tpm = None

        self.cds_sum = None
        self.cds_mean = None
        self.cds_CoV = None
        self.cds_std = None
        self.merge_file = None

        # ttest and ks-test
        self.test_CoV = None


    def read_rpf(self):
        '''
        import the rpf density file and convert them to rpm density

        1. import the rpf density file with import_rpf function, 
        keep all expressed gene.
        2. convert the rpf density to rpm density

        '''
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sites=self.site,
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis = self.tis,
                                      tts = self.tts,
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)
        # self.all_rpf = rpf_results[0]

        self.sample_name = rpf_results[1]
        self.col_index = ['name'] + self.sample_name
        self.sample_num = rpf_results[2]

        self.merged_rpf = rpf_results[3].copy()
        #self.length = self.merged_rpf.loc[:, ['name', 'region']].groupby("name").count()

        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        #self.high_rpf = rpf_results[7].copy()

        if self.norm:
            self.merged_rpf.loc[:, self.sample_name] = self.merged_rpf[self.sample_name] * 1e6 / self.total_rpf_num
            self.high_rpf = self.merged_rpf[self.merged_rpf['name'].isin(self.high_gene)]
        else:
            self.high_rpf = rpf_results[7].copy()

        self.length = self.high_rpf.loc[:, ['name', 'region']].groupby("name").count()
        del rpf_results

    def read_group(self):
        '''
        import the group file, wihich contain these two column
        +----------+---------+
        | name     | group  |
        +==========+=========+
        | wt1      | wt      |
        | wt2      | wt      |
        | treat1   | treat   |
        | treat2   | treat   |
        +----------+---------+

        '''
        self.group = pd.read_csv(self.group_file, sep='\t', header=0, names=None)

    def cal_CoV(self):
        '''
        CoV = (std of each gene) / (mean of each gene)
        '''
        # calculate the mean and std of per gene
        # self.CoV_dict = OrderedDict()
        # for now_sample in self.sample_name:
        self.cds_sum = self.high_rpf.loc[:, self.col_index].groupby("name")[self.sample_name].sum()
        self.cds_mean = self.high_rpf.loc[:, self.col_index].groupby("name")[self.sample_name].mean()
        self.cds_std = self.high_rpf.loc[:, self.col_index].groupby("name")[self.sample_name].std()

        # calculate the CoV of per gene
        np.seterr(divide='ignore', invalid='ignore')
        # np.where(self.cds_mean != 0, np.divide(self.cds_std, self.cds_mean), 0)
        self.cds_CoV = np.divide(self.cds_std, self.cds_mean)
        # self.cds_CoV = self.cds_CoV.fillna(0)

        # add suffix
        self.cds_sum = self.cds_sum.add_suffix('_sum')
        self.cds_mean = self.cds_mean.add_suffix('_mean')
        self.cds_CoV = self.cds_CoV.add_suffix('_CoV')

        # merge the mean file and CoV file
        self.merge_file = pd.concat([self.cds_sum, self.cds_mean, self.cds_CoV], axis=1)
        self.merge_file = self.merge_file.dropna(how='any')

    def output_CoV(self):
        # output the merge file
        file_name = self.output_prefix + '_CoV.txt'
        self.merge_file.to_csv(file_name, sep='\t', index=True)

    def re_combination(self, groups):
        '''
        combine each two groups to one comparasion.
        [a, b, c] ==> [[a, b], [a, c], [b, c]]

        '''
        combinations = []
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                combinations.append((groups[i], groups[j]))
        
        return combinations

    def draw_group_plot(self, CoV1, CoV2, mean1, mean2, group1, group2):

        out_pdf = self.output_prefix + "_CoV_pointplot.pdf"
        out_png = self.output_prefix + "_CoV_pointplot.png"

        log2_mean1 = np.log2(mean1)
        log2_CoV1 = np.log2(CoV1)

        log2_mean2 = np.log2(mean2)
        log2_CoV2 = np.log2(CoV2)

        #f = Fitter(log2_CoV1.to_list(), distributions='nbinom')
        #f.fit()
        #f.summary()

        plt.scatter(log2_mean1, log2_CoV1)
        plt.xlabel('log2(mean)')
        plt.ylabel('log2(CoV)')
        # plt.show()

    def test_group_CoV(self):
        '''
        @Message  : test the CoV of each group
        @Input    : self.group --> the group file
                    self.re_combination --> the function to reterieve the comparasion
        @Return   : self.test_CoV_dict --> the dict contain the CoV of each comparasion
                    self.test_CoV --> the dataframe contain the CoV of each comparasion
        @Flow     : step1 --> reterieve the comparasion
                    step2 --> reterieve the samples and CoV
                    step3 --> calculate the differential pvalue with ks-test and student's-test
                    step4 --> output the result

        '''


        # import the group file
        uniq_group = self.group.Group.unique()
        print('Now groups: ', uniq_group, flush=True)
        
        # split group to each comparasion
        combinations = self.re_combination(uniq_group)

        # reterieve the CoV and calculate the differential pvalue with ks-test and student's-test
        test_CoV_dict = OrderedDict()

        for comparasion in combinations:
            group1 = comparasion[0]
            group2 = comparasion[1]

            # reterieve the samples and CoV
            samples1 = self.group[self.group['Group'] == group1]
            samples2 = self.group[self.group['Group'] == group2]

            CoV1 = self.merge_file[samples1['Name'] + '_CoV'].mean(axis=1)
            CoV2 = self.merge_file[samples2['Name'] + '_CoV'].mean(axis=1)
            
            mean1 = self.merge_file[samples1['Name'] + '_mean'].mean(axis=1)
            mean2 = self.merge_file[samples2['Name'] + '_mean'].mean(axis=1)

            # calculate the differential pvalue with ks-test and student's-test
            t_statistic, t_p_value = stats.ttest_ind(CoV1, CoV2)
            ks_statistic, ks_p_value = stats.ks_2samp(CoV1, CoV2)

            comp_name = group1 + '_' + group2
            test_CoV_dict[comp_name] = [t_statistic, t_p_value, ks_statistic, ks_p_value]

            # draw the figure
            # self.draw_group_plot(CoV1, CoV2, mean1, mean2, group1, group2)

        self.test_CoV = pd.DataFrame(test_CoV_dict).transpose()
        self.test_CoV.columns = ['t_statistic', 't_p_value', 'ks_statistic', 'ks_p_value']

        # output the 
        file_name = self.output_prefix + '_compared_CoV.txt'
        self.test_CoV.to_csv(file_name, sep='\t', index=True)

    def draw_fit_plot(self):
        '''
        @Message  : draw the fit plot
        @Input    : self.merge_file --> the CoV merge file
        @Return   : out_pdf --> the pdf file of fit plot
                    self.group --> the dataframe contain the group information
        @Flow     : step1 --> split the group to each comparasion
                    step2 --> reterieve the samples and CoV
                    step3 --> create the figure with 3 columns and 1 row, containing 3 subplots
                            first subplot --> draw the scatter and fitted line plot for group1, x = log2(mean), y = log2(CoV)
                            second subplot --> draw the scatter and fitted line plot for group2, x = log2(mean), y = log2(CoV)
                            third subplot --> draw the fitted line plot for group1 and group2, x = log2(mean), y = log2(CoV)
                    step4 --> draw the scatter and fitted line plot for each group, x = log2(mean), y = log2(CoV)
                    step5 --> draw the fitted line plot for each group, x = log2(mean), y = log2(CoV)
                    step6 --> output the result
        '''
        
        # define the fitted function
        def fitting_function(mean, alpha, beta):
            # fit the scatter with function log2(CoV) = 0.5 * log2(β/mean + α)
            # y = 0.5 * log2(β/mean + α)
            return 2 ** 0.5 * np.log2(beta / mean + alpha)
        
        # set the initial guess
        initial_guess = [1, 1]

        from scipy.optimize import curve_fit
        
        # split group to each comparasion
        if self.group is None:
            self.group = pd.DataFrame({'Name': self.sample_name, 'Group': self.sample_name})

        uniq_group = self.group['Group'].unique()
        combinations = self.re_combination(uniq_group)
        
        matplotlib.use('AGG')
        # draw the scatter and fitted line plot for each group, x = log2(mean), y = log2(CoV)
        for comparasion in combinations:
            group1 = comparasion[0]
            group2 = comparasion[1]
        
            out_pdf = self.output_prefix + "_" + group1 + "_vs_" + group2 + "_CoV_fitplot.pdf"
            out_png = self.output_prefix + "_" + group1 + "_vs_" + group2 + "_CoV_fitplot.png"

            # create the figure with m + 1 columns and 1 row, containing m + 1 subplots, m = len(comparasion)
            fig, axs = plt.subplots(1, 3, figsize=(20, 5))
            
            # reterieve the samples and CoV
            samples1 = self.group[self.group['Group'] == group1]
            samples2 = self.group[self.group['Group'] == group2]

            CoV1 = self.merge_file[samples1['Name'] + '_CoV'].mean(axis=1)
            CoV2 = self.merge_file[samples2['Name'] + '_CoV'].mean(axis=1)
            
            mean1 = self.merge_file[samples1['Name'] + '_mean'].mean(axis=1)
            mean2 = self.merge_file[samples2['Name'] + '_mean'].mean(axis=1)

            # fit the data
            params1, covariance1 = curve_fit(fitting_function, mean1, CoV1, p0=initial_guess)
            params2, covariance2 = curve_fit(fitting_function, mean2, CoV2, p0=initial_guess)

            # retrieve the standard deviation of the parameters
            alpha_fit1, beta_fit1 = params1
            alpha_fit2, beta_fit2 = params2

            # retrieve the mean and CoV
            mean_fit1 = np.linspace(min(mean1), max(mean1), 100)
            mean_fit2 = np.linspace(min(mean2), max(mean1), 100)

            cov_fit1 = fitting_function(mean_fit1, alpha_fit1, beta_fit1)
            cov_fit2 = fitting_function(mean_fit2, alpha_fit2, beta_fit2)

            axs[0].scatter(np.log2(mean1), np.log2(CoV1), s=10, color='blue', marker='o', label=group1)
            axs[0].set_xlabel('log2(mean)')
            axs[0].set_ylabel('log2(CoV)')
            # axs[0].set_ylim(-1, 5)
            axs[0].set_title(group1)
            axs[0].legend(loc='upper right')
            axs[0].grid(False)

            # draw the second scatter plot and fitted plot
            axs[1].scatter(np.log2(mean2), np.log2(CoV2), s=10, color='orange', marker='o', label=group2)
            axs[1].set_xlabel('log2(mean)')
            axs[1].set_ylabel('log2(CoV)')
            # axs[1].set_ylim(-1, 5)
            axs[1].set_title(group2)
            axs[1].legend(loc='upper right')
            axs[1].grid(False)

            # draw the third plot only contain the fitted line plot
            axs[2].plot(np.log2(mean_fit1), np.log2(cov_fit1) + beta_fit1**0.5, color='blue', label='fitted1')
            axs[2].plot(np.log2(mean_fit2), np.log2(cov_fit2) + beta_fit2**0.5, color='orange', label='fitted2')
            # axs[2].set_ylim(-1, 5)
            axs[2].set_xlabel('log2(mean)')
            axs[2].set_ylabel('log2(CoV)')
            axs[2].set_title(group1 + '_vs_' + group2)
            axs[2].legend(loc='upper right')
            axs[2].grid(False)

            # output the result
            plt.savefig(out_pdf, dpi=300, format='pdf')
            plt.savefig(out_png, dpi=300, format='png')
            plt.close()
