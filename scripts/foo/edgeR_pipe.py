#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : DESeq2_pipe.py


import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import rpy2.robjects as robjects
from rpy2.robjects import Formula
from rpy2.robjects import pandas2ri
import rpy2.robjects.packages as rpackages


pandas2ri.activate()
deseq = rpackages.importr('DESeq2')


class DESeq2Pipe(object):

    def __init__(self, args):
        # options for input and output
        self.rpf = args.rpf
        self.design = args.design
        self.output = args.output

        # options for filter reads and design
        self.min = args.min
        self.ck = args.control.split(',')
        self.treat = args.treat.split(',')
        self.sample_name = None
        self.sample_norm = None

        # options for deseq formula
        self.design_formula = None
        self.mat = None
        self.flt_mat = None
        self.column = None
        self.gene_id = None
        self.count_mat = None
        self.design_mat = None

        # deseq pipeline
        self.dds = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_matrix = None
        self.deseq_result_py = None
        self.output_result = None

        # volcano
        self.logfc = args.logfc
        self.padj = args.padj

    def import_design(self):
        self.design = pd.read_csv(self.design, index_col=False, header=0, sep='\t')

    def import_reads(self):
        self.mat = pd.read_csv(self.rpf, skiprows=0, sep='\t', header=0, index_col=False)
        sample_name = [i.replace('.genes.results', '').replace('.isoforms.results', '') for i in self.mat.columns[1::]]
        sample_name = [
            i.replace('_trim_cds_rpf', '').replace('_cds_rpf', '').replace('_utr5_rpf', '').replace('_utr3_rpf', '')
            for i in sample_name]
        self.mat.columns = ['name'] + sample_name
        self.mat[sample_name] = self.mat[sample_name].astype(int)

    def perepare_deseq_pipe(self):
        # make the design matrix
        design_mat = self.design[self.design.iloc[:, 1].isin(self.ck + self.treat)].iloc[:, 0:2]
        design_mat['now_groups'] = design_mat.iloc[:, 1].apply(lambda x: 'ck' if x in self.ck else 'treat')
        design_mat.drop([design_mat.columns[1]], axis=1, inplace=True)
        design_formula = '~ now_groups'
        self.design_formula = Formula(design_formula)

        # filter the reads matrix
        ck = design_mat[design_mat['now_groups'] == 'ck'].iloc[:, 0].to_list()
        treat = design_mat[design_mat['now_groups'] == 'treat'].iloc[:, 0].to_list()

        # rename for normalize reads
        self.column = self.mat.columns[0]
        self.sample_name = [self.column] + ck + treat
        self.sample_norm = [i + '_norm' for i in self.sample_name[1::]]

        # filter the specified samples and high expression genes
        flt_mat = self.mat[(self.mat[ck].sum(axis=1) > self.min) | (self.mat[treat].sum(axis=1) > self.min)].copy()
        self.flt_mat = flt_mat[self.sample_name].copy()
        self.gene_id = self.flt_mat.iloc[:, 0].copy()

        # make the deseq matrix
        self.count_mat = pandas2ri.py2rpy(self.flt_mat.drop(self.column, axis=1))
        self.design_mat = pandas2ri.py2rpy(design_mat)

        # output the design_mat
        # ck_group = '_'.join(self.ck)
        # treat_group = '_'.join(self.treat)
        # group_file = self.output + '-' + ck_group + '_vs_' + treat_group + '-design.txt'
        group_file = self.output + '_design.txt'
        design_mat.iloc[:, 0] += '_norm'
        design_mat.to_csv(group_file, index=False, header=True, sep="\t")

    def run_deseq2(self):
        # run the deseq2
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_mat, colData=self.design_mat,
                                                design=self.design_formula)
        self.dds = deseq.DESeq(self.dds, fitType='local')
        self.normalized_matrix = pd.DataFrame(deseq.counts_DESeqDataSet(self.dds, normalized=True),
                                              columns=self.sample_norm)

    def get_deseq_result(self):
        now_df = robjects.r('function(x) data.frame(x)')

        self.comparison = deseq.resultsNames(self.dds)
        self.deseq_result = deseq.results(self.dds)
        self.deseq_result_py = now_df(self.deseq_result)
        self.count_mat = now_df(self.count_mat)

    def output_deseq(self):
        self.count_mat.index = self.gene_id.values
        self.deseq_result_py.index = self.gene_id.values
        self.normalized_matrix.index = self.gene_id.values

        self.output_result = pd.concat([self.deseq_result_py, self.count_mat, self.normalized_matrix], axis=1)
        self.output_result.index.name = 'name'
        print(self.output_result)
        self.output_result.to_csv(self.output + '_deseq2.txt', index=True, header=True, na_rep="nan", sep='\t')

    def out_gene_list(self):
        # output sig different gene
        up_genes = self.output_result.loc[self.output_result['expression'] == 'up',].index.tolist()
        down_genes = self.output_result.loc[self.output_result['expression'] == 'down',].index.tolist()
        up_genes = pd.DataFrame(up_genes, columns=['name'])
        down_genes = pd.DataFrame(down_genes, columns=['name'])
        up_genes.to_csv(self.output + '_up_gene.txt', index=False, header=True)
        down_genes.to_csv(self.output + '_down_gene.txt', index=False, header=True)

    def prepare_volcano(self):
        # replace padj=0 to padj = minimum none-zero padj
        min_padj = self.output_result['padj'].apply(lambda x: x if x > 0 else 1).min()
        self.output_result.loc[self.output_result['padj'] == 0, 'padj'] = min_padj
        self.output_result.loc[:, '-log10(padj)'] = -np.log10(self.output_result['padj'])

        # annotate the up/down regulated genes
        def get_flag(rows):
            if rows['log2FoldChange'] >= self.logfc and rows['padj'] < self.padj:
                return 'up'
            elif rows['log2FoldChange'] <= -self.logfc and rows['padj'] < self.padj:
                return 'down'
            else:
                return 'ns'

        self.output_result.loc[:, 'expression'] = self.output_result.apply(lambda x: get_flag(x), axis=1)
        self.out_gene_list()

        counts = self.output_result['expression'].value_counts()
        self.output_result['expression'].replace('up', 'up' + '(' + str(counts['up']) + ')', inplace=True)
        self.output_result['expression'].replace('down', 'down' + '(' + str(counts['down']) + ')', inplace=True)
        self.output_result['expression'].replace('ns', 'ns' + '(' + str(counts['ns']) + ')', inplace=True)

    def draw_volcano(self):
        out_pdf = self.output + '_volcano.pdf'
        out_png = self.output + '_volcano.png'

        import seaborn as sns
        matplotlib.use('AGG')

        fig, ax = plt.subplots(figsize=(5.5, 4))
        # sns.despine(fig, left=True, bottom=True, right=True, top=True)
        clarity_ranking = self.output_result['expression'].sort_values().unique().tolist()
        sns.scatterplot(x="log2FoldChange", y="-log10(padj)",
                        hue="expression", size="baseMean",
                        palette=("#0072b5", "grey", "#bc3c28"),
                        hue_order=clarity_ranking,
                        sizes=(5, 100), linewidth=0,
                        data=self.output_result, ax=ax)
        plt.legend(loc=2, bbox_to_anchor=(1, 1))

        # set the figure range
        xmin = self.output_result['log2FoldChange'].min() * 1.1
        xmax = self.output_result['log2FoldChange'].max() * 1.1
        ymin = -self.output_result['-log10(padj)'].max() * 0.05
        ymax = self.output_result['-log10(padj)'].max() * 1.1
        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='_'.join(self.ck) + '_vs_' + '_'.join(self.treat))

        # draw the threshold line
        ax.vlines(-self.logfc, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.vlines(self.logfc, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
        ax.hlines(-np.log10(self.padj), xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)

        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()
