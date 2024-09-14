#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/8/4 9:59
# @Project : pipe.py
# @Script  : run_DESeq2.py
# @Version : python 3.8.5
# @product : PyCharm
# @Author  : Rensc
# @E-mail  : rensc0718@163.com


import argparse
from argparse import RawTextHelpFormatter

import pandas as pd
import rpy2.robjects as robjects
import yaml
from rpy2.robjects import Formula
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
deseq = importr('DESeq2')


# get arguments and define the scripts usage
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to perform differential analysis with featureCounts output.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )
    parser.add_argument(
        "-i", dest="input", required=True, type=str, help="input reads count file"
    )
    parser.add_argument(
        "-m", dest="min", required=True, type=int, help="specify the minimum reads count"
    )
    parser.add_argument(
        "-d", dest="design", required=True, type=str, help="input whole design file"
    )
    parser.add_argument(
        "-p", dest="pairs", required=True, type=str, help="input the two group pairs"
    )
    parser.add_argument(
        "-g", dest="group", required=True, type=str, help="output the two group pairs design file"
    )
    parser.add_argument(
        "-o", dest="output", required=True, type=str, help="output the DESeq2 results"
    )
    args = parser.parse_args()
    return args


class DESeq2Pipe:

    def __init__(self, count_mat, design_matrix, design_formula, gene_column, sample_name):
        self.dds = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_matrix = None
        self.deseq_result_py = None
        self.output_result = None
        self.gene_column = gene_column
        self.gene_id = count_mat[self.gene_column].copy()
        self.sample_id = [i + '_norm' for i in sample_name[1::]]
        self.count_matrix = pandas2ri.py2rpy(count_mat.drop(self.gene_column, axis=1))
        self.design_matrix = pandas2ri.py2rpy(design_matrix)
        self.design_formula = Formula(design_formula)

    def run_deseq(self):
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix, colData=self.design_matrix,
                                                design=self.design_formula)
        self.dds = deseq.DESeq(self.dds, fitType='local')
        self.normalized_matrix = pd.DataFrame(deseq.counts_DESeqDataSet(self.dds, normalized=True),
                                              columns=self.sample_id)
        self.normalized_matrix[self.gene_column] = self.gene_id.values

    def get_deseq_result(self, **kwargs):
        now_df = robjects.r('function(x) data.frame(x)')

        self.comparison = deseq.resultsNames(self.dds)
        self.deseq_result = deseq.results(self.dds, **kwargs)
        self.deseq_result_py = now_df(self.deseq_result)

    def output_deseq(self, count_mat, file_name):
        self.deseq_result_py[self.gene_column] = self.gene_id.values
        self.output_result = pd.merge(count_mat, self.deseq_result_py, on=self.gene_column)
        self.output_result = pd.merge(self.output_result, self.normalized_matrix, on=self.gene_column)
        self.output_result.to_csv(file_name, index=False, header=True, na_rep = "nan")


def import_design(input_design):
    design = pd.read_csv(input_design, index_col=False, header=0, sep='\t')

    return design


def import_reads(input_file):
    raw_mat = pd.read_csv(input_file, comment='#', sep='\t', header=0, index_col=False)
    #sample_name = [i.split('/')[-1].split('.')[0].replace('Aligned','') for i in raw_mat.columns[6::]]
    sample_name = [i.split('/')[-1].split('.')[0].replace('Aligned','') for i in raw_mat.columns[6::]]
    raw_mat.columns = raw_mat.columns[0:6].tolist() + sample_name

    return raw_mat


def perform_deseq_pipe(raw_mat, design, pairs, min_num, group_file, output_file):

    # get the sample list of different groups
    groups = pairs.split('_vs_')
    wt = groups[0].replace(' ', '').split('+')
    treat = groups[1].replace(' ', '').split('+')
    gene_id = raw_mat.columns[0]

    # make the design matrix
    design_matrix = design[design.iloc[:, 1].isin(wt + treat)].iloc[:, 0:2]
    design_matrix['now_groups'] = design_matrix.iloc[:, 1].apply(lambda x: 'WT' if x in wt else 'Treat')

    design_matrix.to_csv(group_file, index=False, header=True, sep="\t")

    design_matrix.drop([design_matrix.columns[1]], axis=1, inplace=True)
    design_formula = '~ now_groups'

    # filter the reads matrix
    wt_name = design_matrix[design_matrix['now_groups'] == 'WT'].iloc[:, 0].to_list()
    treat_name = design_matrix[design_matrix['now_groups'] == 'Treat'].iloc[:, 0].to_list()

    sample_name = [gene_id] + treat_name + wt_name
    flt_mat = raw_mat[(raw_mat[treat_name].sum(axis=1) > min_num) | (raw_mat[wt_name].sum(axis=1) > min_num)].copy()
    count_mat = flt_mat[sample_name].copy()

    # run the deseq2
    now_comparison = DESeq2Pipe(count_mat, design_matrix, design_formula, gene_id, sample_name)
    now_comparison.run_deseq()
    now_comparison.get_deseq_result()
    now_comparison.output_deseq(count_mat, output_file)


def main():
    # return arguments
    args = parse_args()

    # import the design file
    design = import_design(args.design)

    # import the reads count file
    raw_mat = import_reads(args.input)

    # run the deseq2 for differential analysis
    perform_deseq_pipe(raw_mat, design, args.pairs, args.min, args.group, args.output)


if __name__ == '__main__':
    main()
