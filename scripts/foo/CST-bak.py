#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : CST.py


from collections import OrderedDict
import pandas as pd
import polars as pl
import xlsxwriter
from . import RPFs


class CodonTable(object):
    def __init__(self):
        """
        commonly contains the follow elements:

        codon number, frequency, properties
        ribo-seq number, rna-seq number
        valid codon number
        codon density proportion
        """
        self.codon = None
        self.name = None
        self.aa = None
        self.number = 0
        self.frequency = 0
        self.cai = 0
        self.cbi = 0
        self.rscu = 0
        self.properties = None

        self.rpf_codon_sum = 0
        self.rpf_codon_val = 0
        self.rpf_sum = 0
        self.rpf_proportion = 0

        self.rna_codon_sum = 0
        self.rna_codon_val = 0
        self.rna_sum = 0
        self.rna_proportion = []

        self.cst = []


class GeneTable(object):
    def __init__(self, name):
        """
        gene: gene name
        codon count dict: (codon: number)
        codon_num: total codon number of each gene
        rpf, rna reads: Ribo-seq and RNA-seq reads density
        density: (TE = rpf/rna)
        elongation rate: [codon_num / Σ(count * CST))]
        initiate rate: (elongation rate * density)
        """
        self.gene = name
        self.length = 0
        self.codon_num = 0
        self.elongate_rate = 0
        self.rpf = 0
        self.rna = 0
        self.density = 0
        self.initiate_rate = 0
        self.codon_count = OrderedDict()


class CodonSelectiveTime(object):

    def __init__(self, args):
        """
        receive arguments and set the self parameters
        """
        self.gene = args.list
        self.rpf_file = args.rpf
        self.rna_file = args.rna

        self.rpf_sample = None
        self.rna_sample = None

        # for the high expression gene filter
        self.gene_list = args.list
        self.rpf_num = args.min

        self.rpf_gene_sum = None
        self.rna_gene_sum = None

        self.high_rpf = None
        self.high_rna = None

        self.rpf_gene = None
        self.rna_gene = None

        self.overlap_gene = None
        self.high_gene = None

        self.high_rpf_num = 0
        self.high_rna_num = 0

        self.total_rpf_num = 0
        self.total_rna_num = 0

        # filter the stable codon coverage range
        self.site = args.site
        self.frame = args.frame
        self.tis = args.tis
        self.tts = args.tts
        self.times = args.times

        # table for codon occupancy
        self.codon_table = OrderedDict()
        self.gene_table = OrderedDict()
        self.absolute_cst = None
        self.relative_cst = None
        self.output = args.output

    def import_density(self):
        """
        import the Ribo-seq and RNA-seq density file.
        By default, the two data names are different.

        return the main results as follows:
        samples name,
        high expression gene density
        high expression gene name
        total rpf/reads
        """

        # import the rpf density
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file, sample_num=None,
                                      sample_name=None, sites=self.site, frame=self.frame,
                                      gene=None, rpf_num=self.rpf_num, tis=self.tis, tts=self.tts)
        self.rpf_sample = rpf_results[1]
        self.total_rpf_num = rpf_results[4]
        self.rpf_gene_sum = rpf_results[5]
        self.rpf_gene = rpf_results[6]
        self.high_rpf = rpf_results[7].copy()

        del rpf_results

        # import the rna density
        rna_results = RPFs.import_rpf(rpf_file=self.rna_file, sample_num=None,
                                      sample_name=None, sites=self.site, frame=self.frame,
                                      gene=None, rpf_num=self.rpf_num, tis=self.tis, tts=self.tts)
        self.rna_sample = rna_results[1]
        self.total_rna_num = rna_results[4]
        self.rna_gene_sum = rna_results[5]
        self.rna_gene = rna_results[6]
        self.high_rna = rna_results[7].copy()

        del rna_results

    def merge_rpf_rna(self):
        """
        get the overlap gene from rpf and mrna density file
        only these genes will submit to follow analysis
        """
        self.overlap_gene = set(self.rpf_gene) & set(self.rna_gene)

        self.high_rpf = self.high_rpf.loc[self.high_rpf['name'].isin(self.overlap_gene)]
        self.high_rna = self.high_rna.loc[self.high_rna['name'].isin(self.overlap_gene)]

        self.high_rpf_num = self.high_rpf[self.rpf_sample].sum()
        self.high_rna_num = self.high_rna[self.rna_sample].sum()

        self.high_gene = pd.merge(self.high_rpf, self.high_rna, how='inner')

        self.high_rpf = pl.DataFrame(self.high_rpf)
        self.high_rna = pl.DataFrame(self.high_rna)

    def make_codon_table(self):
        """
        make the codon message table contains: codon, name, amino-acid
        """
        codon_dict = RPFs.codon_table()[0]
        for codon, mess in codon_dict.items():
            self.codon_table[codon] = CodonTable()
            self.codon_table[codon].codon = codon
            self.codon_table[codon].name = mess[0]
            self.codon_table[codon].aa = mess[1]

    def ribo_codon_frequency(self):
        """
        calculate the mean rpf density of codon
        """
        for codon in self.codon_table.keys():
            codon_flt = self.high_rpf.filter(pl.col("aa_list") == codon)
            codon_sum = codon_flt[self.rpf_sample].sum(axis=0).to_series()[0]
            codon_proportion = float(codon_sum / self.high_rpf_num)

            self.codon_table[codon].rpf_codon_sum = len(codon_flt)
            self.codon_table[codon].rpf_codon_val = (codon_flt[self.rpf_sample] > 0).sum().to_series()[0]
            self.codon_table[codon].rpf_sum = codon_sum
            self.codon_table[codon].rpf_proportion = codon_proportion

    def mrna_codon_frequency(self, times=0, initiate_rate=None):
        """
        calculate the mean reads density of codon
        """
        if initiate_rate:
            high_rna = self.high_rna.with_column(
                self.high_rna[self.rna_sample] * self.high_rna['name'].apply(
                    lambda x: self.gene_table[x].initiate_rate if x in self.gene_table.keys() else 1))
        else:
            high_rna = self.high_rna

        for codon in self.codon_table.keys():
            codon_flt = high_rna.filter(pl.col("aa_list") == codon)
            codon_sum = codon_flt[self.rna_sample].sum(axis=0).to_series()[0]
            self.codon_table[codon].rna_codon_sum = len(codon_flt)
            self.codon_table[codon].rna_codon_val = (codon_flt[self.rna_sample] > 0).sum().to_series()[0]
            self.codon_table[codon].rna_sum = codon_sum

            rna_proportion = float(codon_sum / self.high_rna_num)
            self.codon_table[codon].rna_proportion.append(rna_proportion)

            cst = self.codon_table[codon].rpf_proportion / self.codon_table[codon].rna_proportion[times]
            self.codon_table[codon].cst.append(cst)

    def make_gene_table(self):
        """
        make high expression gene table with the following parameters:
        """

        gene_length_table = self.high_gene.groupby(['name']).size()
        gene_codon_table = self.high_gene.groupby(['name', 'aa_list']).size()
        # gene_table = gene_table.reset_index()

        for gene in self.overlap_gene:
            self.gene_table[gene] = GeneTable(gene)
            self.gene_table[gene].codon_count = dict(gene_codon_table[gene])
            self.gene_table[gene].length = gene_length_table[gene] * 3
            self.gene_table[gene].codon_num = gene_length_table[gene]
            self.gene_table[gene].rpf = float(self.rpf_gene_sum.loc[gene])
            self.gene_table[gene].rna = float(self.rna_gene_sum.loc[gene])
            self.gene_table[gene].density = self.gene_table[gene].rpf / self.gene_table[gene].rna

    def codon_elongate_rate(self, codon_num, codon_dict, times=0):
        """
        calculate the elongate rate for each codon of gene
        """
        codon_density = 0
        for codon, number in codon_dict.items():
            try:
                codon_density += codon_dict[codon] * self.codon_table[codon].cst[times]
            except KeyError:
                pass

        elongate_rate = codon_num / codon_density
        return elongate_rate

    def compute_cst_iterative(self):
        """
        ## compute cst iteratively
        elongate_rate = [codon_num / Σ(count * CST))]

        ## normalize the rna density with initiate_rate
        rna_density = rna_density * initiate_rate

        ## compute cst iteratively ()
        mrna_codon_proportion = Σ(codon) / Σ(gene)
        cst = rpf_codon_proportion / mrna_codon_proportion

        """

        for iters in range(0, self.times, 1):
            print('Times of iteration: {}'.format(iters + 1), flush=True)
            for gene in self.overlap_gene:
                self.gene_table[gene].elongate_rate = self.codon_elongate_rate(self.gene_table[gene].codon_num,
                                                                               self.gene_table[gene].codon_count,
                                                                               times=iters)
                self.gene_table[gene].initiate_rate = self.gene_table[gene].elongate_rate * self.gene_table[
                    gene].density

            # self.high_rna = self.high_rna.with_column(
            #     self.high_rna[self.rna_sample] * self.high_rna['name'].apply(
            #     lambda x: self.gene_table[x].initiate_rate if x in self.gene_table.keys() else 1))

            self.mrna_codon_frequency(times=iters + 1, initiate_rate=True)

    def format_cst_results(self):
        absolute_cst = OrderedDict()

        for codon in self.codon_table.keys():
            absolute_cst[codon] = [codon,
                                   self.codon_table[codon].name,
                                   self.codon_table[codon].aa,
                                   self.codon_table[codon].rpf_codon_val,
                                   self.codon_table[codon].rpf_sum,
                                   self.codon_table[codon].rpf_proportion,
                                   self.codon_table[codon].rpf_codon_val,
                                   self.codon_table[codon].rna_sum,
                                   self.codon_table[codon].rna_proportion[0],
                                   self.codon_table[codon].rna_proportion[-1]
                                   ]
            absolute_cst[codon].extend(self.codon_table[codon].cst)
        self.absolute_cst = pd.DataFrame.from_dict(absolute_cst).transpose()

        column_name = ['codon', 'name', 'AA', 'rpf_codon_val', 'rpf_sum', 'rpf_proportion',
                       'rna_codon_val', 'rna_sum', 'rna_proportion_first', 'rna_proportion_last']
        cst_name = ['cst_' + str(iters) for iters in range(self.times + 1)]
        column_name.extend(cst_name)
        self.absolute_cst.columns = column_name

        self.relative_cst = self.absolute_cst.copy()
        self.relative_cst[cst_name] = self.relative_cst[cst_name] / self.relative_cst[cst_name].max()

    def output_cst(self):
        # cst_writer = pd.ExcelWriter(self.output + '_absolute_cst.xlsx', engine='xlsxwriter')
        with pd.ExcelWriter(self.output + '_codon_selection_time.xlsx',) as cst_writer:
            self.absolute_cst.to_excel(cst_writer, sheet_name='absolute_cst', index=False)
            self.relative_cst.to_excel(cst_writer, sheet_name='relative_cst', index=False)

    def draw_cst_line(self):
        import plotly.express as px
        out_pdf = self.output + "_cst_line.pdf"
        out_png = self.output + "_cst_line.png"

        rel_cst = self.relative_cst.melt(id_vars=['codon', 'name'],
                                         var_name='iteration',
                                         value_name='cst',
                                         value_vars=['cst_' + str(iters) for iters in range(self.times + 1)])

        fig = px.line(rel_cst, x = "iteration", y = "cst", title = 'relative codon selection time')
        fig.show()
