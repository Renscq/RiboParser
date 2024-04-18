#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : TranslationEfficiency.py


import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests


class TranslationEfficiency(object):

    def __init__(self, args):
        # input and output
        self.ribo = args.ribo
        self.mrna = args.mrna

        self.ribo_design = args.ribo_design
        self.mrna_design = args.mrna_design
        self.min = args.min
        self.logfc = args.logfc
        self.pvalue = args.pvalue
        self.type = args.type
        self.output = args.output

        # gene expression
        self.ribo_df = None
        self.mrna_df = None
        self.merge_df = None
        self.te = None

        # groups
        self.ribo_design_df = None
        self.mrna_design_df = None
        self.ribo_ck = None
        self.ribo_treat = None
        self.mrna_ck = None
        self.mrna_treat = None

    def import_riboseq(self):
        self.ribo_df = pd.read_csv(self.ribo, header=0, names=None, sep='\t')
        colnames = ['name'] + list(self.ribo_df.columns[1:])
        self.ribo_df.columns = colnames
        self.ribo_df.loc[:, 'ave'] = self.ribo_df.filter(regex='_norm').mean(axis=1)
        self.ribo_df = self.ribo_df.loc[self.ribo_df['ave'] >= self.min, colnames]

    def import_rnaseq(self):
        self.mrna_df = pd.read_csv(self.mrna, header=0, names=None, sep='\t')
        colnames = ['name'] + list(self.mrna_df.columns[1:])
        self.mrna_df.columns = colnames
        self.mrna_df.loc[:, 'ave'] = self.mrna_df.filter(regex='_norm').mean(axis=1)
        self.mrna_df = self.mrna_df.loc[self.mrna_df['ave'] >= self.min, colnames]

    def merge_logfc(self):

        def get_class_anno(rows):
            if rows['mrna_logfc'] <= -self.logfc and rows['ribo_logfc'] >= self.logfc:
                return 'A: mRNA-down, Ribo-up'
            elif -self.logfc < rows['mrna_logfc'] < self.logfc and rows['ribo_logfc'] >= self.logfc:
                return 'B: mRNA-ns, Ribo-up'
            elif rows['mrna_logfc'] >= self.logfc and rows['ribo_logfc'] >= self.logfc:
                return 'C: mRNA-up, Ribo-up'

            elif rows['mrna_logfc'] <= -self.logfc and -self.logfc < rows['ribo_logfc'] < self.logfc:
                return 'D: mRNA-down, Ribo-ns.'
            elif -self.logfc < rows['mrna_logfc'] < self.logfc and -self.logfc < rows['ribo_logfc'] < self.logfc:
                return 'E: mRNA-ns, Ribo-ns.'
            elif rows['mrna_logfc'] >= self.logfc and -self.logfc < rows['ribo_logfc'] < self.logfc:
                return 'F: mRNA-up, Ribo-ns.'

            elif rows['mrna_logfc'] <= -self.logfc and rows['ribo_logfc'] <= -self.logfc:
                return 'G: mRNA-down, Ribo-down'
            elif -self.logfc < rows['mrna_logfc'] < self.logfc and rows['ribo_logfc'] <= -self.logfc:
                return 'H: mRNA-ns, Ribo-down'
            elif rows['mrna_logfc'] >= self.logfc and rows['ribo_logfc'] <= -self.logfc:
                return 'I: mRNA-up, Ribo-down'
            else:
                pass

        # detect the algorithm
        if "FDR" in self.ribo_df.columns:
            if self.type == 'pvalue':
                now_col = ['name', 'logFC', 'PValue']
            elif self.type == 'padj':
                now_col = ['name', 'logFC', 'FDR']
            else:
                print('specify the pvalue or padj')
                sys.exit()
        elif "log2FoldChange" in self.ribo_df.columns:
            if self.type == 'pvalue':
                now_col = ['name', 'log2FoldChange', 'padj']
            elif self.type == 'padj':
                now_col = ['name', 'log2FoldChange', 'pvalue']
            else:
                print('specify the pvalue or padj')
                sys.exit()
        else:
            print('unknown algorithm')

        # retrieve the logfc and padj
        ribo_df = self.ribo_df.loc[:, now_col]
        ribo_df.columns = ['name', 'ribo_logfc', 'ribo_padj']
        mrna_df = self.mrna_df.loc[:, now_col]
        mrna_df.columns = ['name', 'mrna_logfc', 'mrna_padj']

        # merge ribo and mrna
        self.merged_df = pd.merge(ribo_df, mrna_df, how='inner', on='name')
        self.merged_df.dropna(axis=0, how='any', inplace=True)

        # annotated the different class
        self.merged_df.loc[:, 'class'] = self.merged_df.apply(lambda x: get_class_anno(x), axis=1)

    def draw_co_diff_genes(self, sig):
        # filter the significant different translation genes
        if sig == 'sig':
            out_pdf = self.output + '_co_sig_DEGs.pdf'
            out_png = self.output + '_co_sig_DEGs.png'
            co_diff_genes = self.merged_df.loc[(self.merged_df['ribo_padj'] < self.pvalue) & (self.merged_df['mrna_padj'] < self.pvalue)].copy()
        else:
            out_pdf = self.output + '_co_DEGs.pdf'
            out_png = self.output + '_co_DEGs.png'
            co_diff_genes = self.merged_df.copy()

        # annotate the number of class
        counts = co_diff_genes['class'].value_counts()
        now_class = pd.DataFrame(counts)
        now_class.loc[:, "class"] = now_class.index + "(" + now_class["class"].astype(str) + ")"
        co_diff_genes['class'] = co_diff_genes['class'].apply(lambda x: now_class.loc[x])
        co_diff_genes.sort_values(['class'], inplace=True)

        # draw the figure
        matplotlib.use('AGG')

        corr_value = co_diff_genes[['ribo_logfc', 'mrna_logfc']].corr().iloc[1, 0]

        fig, ax = plt.subplots(figsize=(8, 5))
        # sns.despine(fig, left=True, bottom=True, right=True, top=True)
        clarity_ranking = co_diff_genes['class'].sort_values().unique().tolist()
        sns.scatterplot(x="mrna_logfc", y="ribo_logfc", hue="class", palette="deep", hue_order=clarity_ranking, sizes=(5, 100), linewidth=0, data=co_diff_genes, ax=ax)

        # set the figure range
        xmax = co_diff_genes['mrna_logfc'].abs().max() * 1.1
        xmin = -xmax
        ymax = co_diff_genes['ribo_logfc'].abs().max() * 1.1
        ymin = -ymax

        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='correlation of DEGs (R=' + str(corr_value) + ')')

        # draw the threshold line
        ax.vlines(-self.logfc, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.vlines(self.logfc, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.hlines(-self.logfc, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.hlines(self.logfc, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)

        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, ncol=1)

        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()

    def import_design(self):
        # import the riboseq design
        self.ribo_design_df = pd.read_csv(self.ribo_design, header=0, names=None, sep='\t')
        self.ribo_design_df.columns = ['name', 'ribo_group', 'file']
        self.ribo_design_df['name'] = [sp.replace('_norm', '_ribo') for sp in list(self.ribo_design_df.iloc[:, 0])]
        self.ribo_ck = self.ribo_design_df[self.ribo_design_df.iloc[:, 1] == 'ck']
        self.ribo_treat = self.ribo_design_df[self.ribo_design_df.iloc[:, 1] == 'treat']

        # import the rnaseq design
        self.mrna_design_df = pd.read_csv(self.mrna_design, header=0, names=None, sep='\t')
        self.mrna_design_df.columns = ['name', 'mrna_group', 'file']
        self.mrna_design_df['name'] = [sp.replace('_norm', '_mrna') for sp in list(self.mrna_design_df.iloc[:, 0])]
        self.mrna_ck = self.mrna_design_df[self.mrna_design_df.iloc[:, 1] == 'ck']
        self.mrna_treat = self.mrna_design_df[self.mrna_design_df.iloc[:, 1] == 'treat']

        # retrieve and rename the expression table
        self.ribo_df.columns = [sp.replace('_norm', '_ribo') for sp in self.ribo_df.columns]
        self.mrna_df.columns = [sp.replace('_norm', '_mrna') for sp in self.mrna_df.columns]

        self.ribo_df = self.ribo_df[['name'] + list(self.ribo_design_df.iloc[:, 0])]
        self.mrna_df = self.mrna_df[['name'] + list(self.mrna_design_df.iloc[:, 0])]

    @staticmethod
    def run_ranksum(rows, te_ck_name, te_treat_name):
        flag, pvalue = stats.ranksums(x=rows[te_ck_name], y=rows[te_treat_name])
        return pvalue

    @staticmethod
    def run_ttest(rows, te_ck_name, te_treat_name):
        flag, pvalue = stats.ttest_ind(rows[te_ck_name], rows[te_treat_name])
        return pvalue

    def calc_te(self):
        # retrieve riboseq
        ribo_ck = list(self.ribo_ck.iloc[:, 0])
        ribo_ck_df = self.ribo_df[['name'] + ribo_ck]

        ribo_treat = list(self.ribo_treat.iloc[:, 0])
        ribo_treat_df = self.ribo_df[['name'] + ribo_treat]

        # retrieve rnaseq
        mrna_ck = list(self.mrna_ck.iloc[:, 0])
        mrna_ck_df = self.mrna_df[['name'] + mrna_ck]

        mrna_treat = list(self.mrna_treat.iloc[:, 0])
        mrna_treat_df = self.mrna_df[['name'] + mrna_treat]

        # calculate the te of ck groups
        merge_ck = pd.merge(ribo_ck_df, mrna_ck_df, on='name', how='inner')
        tmp_te_ck = []
        for i in range(self.mrna_ck.__len__()):
            tmp = merge_ck[ribo_ck].div(merge_ck[mrna_ck[i]], axis=0)
            tmp_te_ck.append(tmp)

        te_ck = pd.concat(tmp_te_ck, axis=1)
        te_ck.index = merge_ck['name'].copy()

        # calculate the te of treat groups
        merge_treat = pd.merge(ribo_treat_df, mrna_treat_df, on='name', how='inner')
        tmp_te_treat = []
        for i in range(self.mrna_treat.__len__()):
            tmp = merge_treat[ribo_treat].div(merge_treat[mrna_treat[i]], axis=0)
            tmp_te_treat.append(tmp)

        te_treat = pd.concat(tmp_te_treat, axis=1)
        te_treat.index = merge_treat['name'].copy()

        # merge total te
        te_ck_len = te_ck.shape[1]
        te_treat_len = te_treat.shape[1]
        te_ck_name = ['ck_' + str(i) for i in range(te_ck_len)]
        te_treat_name = ['treat_' + str(i) for i in range(te_treat_len)]
        merge_te = pd.concat([te_ck, te_treat], axis=1)
        merge_te.columns = te_ck_name + te_treat_name
        merge_te.replace(np.inf, np.nan, inplace=True)
        # merge_te.dropna(axis=0, how='any', inplace=True)

        # ttest for different gene translation
        # merge_te.loc[:, 'pvalue'] = merge_te.apply(lambda x: self.run_ranksum(x, te_ck_name, te_treat_name), axis=1)
        merge_te.loc[:, 'pvalue'] = merge_te.apply(lambda x: self.run_ttest(x, te_ck_name, te_treat_name), axis=1)

        # calculate the log fold change
        merge_te.loc[:, 'ck_mean_te'] = merge_te[te_ck_name].mean(axis=1)
        merge_te.loc[:, 'treat_mean_te'] = merge_te[te_treat_name].mean(axis=1)
        merge_te.loc[:, 'fc'] = merge_te.loc[:, 'ck_mean_te'].div(merge_te.loc[:, 'treat_mean_te'])
        merge_te.loc[:, 'log2FoldChange'] = np.log2(merge_te.loc[:, 'fc'])

        # calculate the log fold change
        fdr_test = multipletests(merge_te.loc[:, 'pvalue'], alpha=0.05, method='bonferroni')
        merge_te.loc[:, 'padj'] = fdr_test[1]
        merge_te.loc[:, 'flag'] = fdr_test[0]

        merge_te = merge_te.reset_index()
        self.te = pd.merge(self.merged_df, merge_te.loc[:, ['name', 'ck_mean_te', 'treat_mean_te', 'fc', 'log2FoldChange', 'pvalue', 'padj', 'flag']], on='name', how='outer')

    def output_te(self):

        def get_flag(rows):
            if rows['log2FoldChange'] >= self.logfc and rows[self.type] < self.pvalue:
                return 'up'
            elif rows['log2FoldChange'] <= -self.logfc and rows[self.type] < self.pvalue:
                return 'down'
            else:
                return 'ns'

        self.te.loc[:, 'te'] = self.te.apply(lambda x: get_flag(x), axis=1)

        # output whole translation efficiency
        out_te_file = self.output + '_te.txt'
        out_te = pd.merge(self.te, self.ribo_df, how='inner', on='name')
        out_te = pd.merge(out_te, self.mrna_df, how='inner', on='name')

        out_te.to_csv(out_te_file, index=False, header=True, sep='\t')

        # output up translation efficiency genes
        out_up = self.output + '_te_up.txt'
        up_genes = self.te.loc[self.te['te'] == 'up', 'name']
        up_genes.to_csv(out_up, index=False, header=True, sep='\t')

        # output down translation efficiency genes
        out_down = self.output + '_te_down.txt'
        down_genes = self.te.loc[self.te['te'] == 'down', 'name']
        down_genes.to_csv(out_down, index=False, header=True, sep='\t')

    def draw_cdf(self):
        # drop the low expression genes
        # self.te.replace([np.inf, -np.inf], np.nan)
        # self.te.dropna(axis=0)

        # get ribo expression
        ribo_ck = self.ribo_df.loc[:, self.ribo_ck['name']].mean(axis=1)
        ribo_treat = self.ribo_df.loc[:, self.ribo_treat['name']].mean(axis=1)
        ribo = pd.concat([ribo_ck, ribo_treat], axis=1)

        # get mrna expression
        mrna_ck = self.mrna_df.loc[:, self.mrna_ck['name']].mean(axis=1)
        mrna_treat = self.mrna_df.loc[:, self.mrna_treat['name']].mean(axis=1)
        mrna = pd.concat([mrna_ck, mrna_treat], axis=1)

        # get translation efficiency
        te_ck = self.te['ck_mean_te']
        te_treat = self.te['treat_mean_te']
        te = pd.concat([te_ck, te_treat], axis=1)

        # draw the cumulative of gene expression and translation efficiency
        out_pdf = self.output + "_te_cdfplot.pdf"
        out_png = self.output + "_te_cdfplot.png"

        matplotlib.use('AGG')
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 3.3), dpi=300)
        ax1 = sns.ecdfplot(data=np.log2(mrna), ax=ax[0])
        ax2 = sns.ecdfplot(data=np.log2(ribo), ax=ax[1])
        ax3 = sns.ecdfplot(data=np.log2(te), ax=ax[2])

        ax1.legend(labels=['ck', 'treat'], loc='lower right')
        ax2.legend(labels=['ck', 'treat'], loc='lower right')
        ax3.legend(labels=['ck', 'treat'], loc='lower right')

        ax[0].set_title('RNA-seq')
        ax[1].set_title('Ribo-seq')
        ax[2].set_title('TE')

        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()

    def draw_volcano(self):
        # annotate the up/down regulated genes
        self.te.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.te.dropna(axis=0, inplace=True)
        counts = self.te['te'].value_counts()
        try:
            self.te['te'].replace('up', 'up' + '(' + str(counts['up']) + ')', inplace=True)
        except KeyError:
            pass
        try:
            self.te['te'].replace('down', 'down' + '(' + str(counts['down']) + ')', inplace=True)
        except KeyError:
            pass
        try:
            self.te['te'].replace('ns', 'ns' + '(' + str(counts['ns']) + ')', inplace=True)
        except KeyError:
            pass

        # replace padj=0 to padj = minimum none-zero padj
        min_padj = self.te[self.type].apply(lambda x: x if x > 0 else 1).min()
        self.te.loc[self.te[self.type] == 0, self.type] = min_padj
        self.te.loc[:, '-' + self.type] = -np.log10(self.te[self.type])
        self.te.loc[:, 'baseMean'] = self.te['ck_mean_te'] + self.te['treat_mean_te']

        # draw the translation efficiency volcano plot
        out_pdf = self.output + '_volcano.pdf'
        out_png = self.output + '_volcano.png'

        import seaborn as sns
        matplotlib.use('AGG')

        fig, ax = plt.subplots(figsize=(5.5, 4))
        # sns.despine(fig, left=True, bottom=True, right=True, top=True)
        clarity_ranking = self.te['te'].sort_values().unique().tolist()
        sns.scatterplot(x="log2FoldChange", y="-" + self.type, hue="te", size="baseMean", palette=("#0072b5", "grey", "#bc3c28"), hue_order=clarity_ranking, sizes=(5, 100), linewidth=0, data=self.te, ax=ax)
        plt.legend(loc=2, bbox_to_anchor=(1, 1))

        # set the figure range
        xmin = self.te['log2FoldChange'].min() * 1.1
        xmax = self.te['log2FoldChange'].max() * 1.1
        ymin = -self.te['-' + self.type].max() * 0.05
        ymax = self.te['-' + self.type].max() * 1.1
        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='translation efficiency')

        # draw the threshold line
        ax.vlines(-self.logfc, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.vlines(self.logfc, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.hlines(-np.log10(self.pvalue), xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)

        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()
