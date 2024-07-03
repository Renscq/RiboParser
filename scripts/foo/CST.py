#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : CST.py


import sys
from collections import OrderedDict
import pandas as pd
import polars as pl
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
# import xlsxwriter
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
        self.aa = None
        self.abbr = None
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

        self.merge_rpf = None
        self.merge_rna = None

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
        self.scale = args.scale

        self.tis = args.tis
        self.tts = args.tts
        self.times = args.times

        # table for codon selection time
        self.stop = args.stop
        self.codon_dict = RPFs.codon_table()[0]
        
        if self.stop:
            del self.codon_dict['TAA']
            del self.codon_dict['TAG']
            del self.codon_dict['TGA']

        self.codon_table = OrderedDict()

        self.gene_table = OrderedDict()
        self.gene_length_table = None
        self.gene_codon_table = None

        self.detail_cst = None
        self.merge_cst = None
        self.absolute_cst = None
        self.relative_cst = None
        self.output = args.output


    def import_density(self):
        '''
        @Message  : function for import the Ribo-seq and RNA-seq density file.
                    By default, the two data names are different.
        @Input    : rpf_file --> rpf/rna file
                    sample_num --> sample number
                    sample_name --> sample name
                    sites --> filter the stable codon coverage range
                    frame --> filter the stable codon coverage frame
                    gene --> high expression gene list
                    rpf_num --> filter the low expression gene
                    tis --> filter the TIS region start
                    tts --> filter the TTS region stop
        @Return   : rpf_sample --> samples name
                    high_rpf --> high expression gene density
                    rpf_gene --> high expression gene name
                    total_rpf_num --> total rpf/reads
        @Flow     : step1 --> import the rpf density
                    step2 --> import the rna density
                    step3 --> check the samples number
        '''
        
        # import the rpf density
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file, sample_num=None,
                                      sample_name=None, sites=self.site, frame=self.frame,
                                      gene=None, rpf_num=1, tis=self.tis, tts=self.tts)
        self.rpf_sample = rpf_results[1]
        self.merge_rpf = pl.DataFrame(rpf_results[3])
        self.total_rpf_num = rpf_results[4]
        self.rpf_gene_sum = rpf_results[5]
        self.rpf_gene = rpf_results[6]
        # self.high_rpf = rpf_results[7].copy()

        del rpf_results

        print('', flush=True)

        # import the rna density
        rna_results = RPFs.import_rpf(rpf_file=self.rna_file, sample_num=None,
                                      sample_name=None, sites=self.site, frame=self.frame,
                                      gene=None, rpf_num=1, tis=self.tis, tts=self.tts)
        self.rna_sample = rna_results[1]
        self.merge_rna = pl.DataFrame(rna_results[3])
        self.total_rna_num = rna_results[4]
        self.rna_gene_sum = rna_results[5]
        self.rna_gene = rna_results[6]
        # self.high_rna = rna_results[7].copy()

        del rna_results
        
        # check the samples number
        rpf_sample_num = len(self.rpf_sample)
        rna_sample_num = len(self.rna_sample)

        if rpf_sample_num != rna_sample_num:
            print('''
                    The number of samples in Ribo-seq and RNA-seq is not equal, please check the input density file.
                  ''', flush=True)
            sys.exit(1)


    def merge_rpf_rna(self):
        '''
        @Message  : get the overlap gene from rpf and mrna density file.
                    only these genes will submit to follow analysis.
                    PS: the gene name should be the same in the two density file.
        @Input    : rpf_gene/rna_gene --> high expression gene name
                    high_rpf/high_rna --> high expression gene density
        @Return   : high_rpf/high_rna --> filtered high expression gene density
                    high_rpf_num/high_rna_num --> total rpf/reads
                    high_gene --> merged high expression gene density
        @Flow     : step1 --> run
        '''
        
        # get the overlap gene from rpf and mrna density file
        overlap_gene = set(self.rpf_gene) & set(self.rna_gene)

        self.high_rpf = self.high_rpf.loc[self.high_rpf['name'].isin(overlap_gene)]
        self.high_rna = self.high_rna.loc[self.high_rna['name'].isin(overlap_gene)]

        self.high_rpf_num = self.high_rpf[self.rpf_sample].sum()
        self.high_rna_num = self.high_rna[self.rna_sample].sum()

        # reset the index and sort the high expression gene
        self.high_rpf.reset_index(drop=True, inplace=True)
        self.high_rna.reset_index(drop=True, inplace=True)

        # self.high_rpf = self.high_rpf.sort_values(by=['name', 'now_nt'])
        # self.high_rna = self.high_rna.sort_values(by=['name', 'now_nt'])
        # self.high_gene = pd.merge(self.high_rpf, self.high_rna, on=['name', 'now_nt', 'from_tis', 'from_tts', 'region', 'codon'])
        # self.high_gene = pd.concat([self.high_rpf, self.high_rna.iloc[:, 6::]], axis=1)

        # self.high_rpf = pl.DataFrame(self.high_rpf)
        # self.high_rna = pl.DataFrame(self.high_rna)


    def make_codon_table(self, rpf_sample=None):

        '''
        @Message  :  make the codon message table contains: codon, name, amino-acid.
        @Input    : codon_dict --> codon message
                    rpf_sample --> sample name
        @Return   : codon_table --> codon message table
        @Flow     : step1 --> create the class for each codon
                    step2 --> merge the codon message
        '''
        
        sp_codon_table = OrderedDict()

        for codon, mess in self.codon_dict.items():
            sp_codon_table[codon] = CodonTable()
            sp_codon_table[codon].codon = codon
            sp_codon_table[codon].aa = mess[0]
            sp_codon_table[codon].abbr = mess[1]

        self.codon_table[rpf_sample] = sp_codon_table


    def filter_rpf_rna_expr(self, rpf_sample=None, rna_sample=None):
        '''
        @Message  : filter the high expression gene > self.rpf_num.
        @Input    : high_rpf --> high expression gene density
                    high_rna --> high expression gene density
                    rpf_sample --> sample name
                    rna_sample --> sample name
                    self.rpf_num --> filter the low expression gene
        @Return   : high_rpf --> filtered high expression gene density
                    high_rna --> filtered high expression gene density
                    overlap_gene --> high expression gene name in the Ribo-seq and RNA-seq data
                    high_rpf_num --> total rpf/reads
                    high_rna_num --> total rna/reads
        @Flow     : step1 --> filter the high expression gene with rpf > self.rpf_num
                    step2 --> check the high expression gene
                    step3 --> filter the high expression gene density
                    step4 --> summary the total rpf/reads and rna/reads
        '''
        
        # filter the high expression gene > self.rpf_num

        rpf_genes = self.rpf_gene_sum[self.rpf_gene_sum[rpf_sample] >= self.rpf_num].index
        rna_genes = self.rna_gene_sum[self.rna_gene_sum[rna_sample] >= self.rpf_num].index
        self.overlap_gene = set(rna_genes) & set(rpf_genes)

        # check the high expression gene
        if len(self.overlap_gene) == 0:
            print('No fitted gene in the Ribo-seq and RNA-seq data.', flush=True)
            sys.exit(1)
        else:
            print('The number of gene in the Ribo-seq and RNA-seq data is: {}'.format(len(self.overlap_gene)), flush=True)

        # filter the high expression gene density
        self.high_rpf = self.merge_rpf.filter(pl.col('name').is_in(self.overlap_gene))
        self.high_rna = self.merge_rna.filter(pl.col('name').is_in(self.overlap_gene))

        # self.high_rpf = self.high_rpf[self.high_rpf['name'].isin(self.overlap_gene)]
        # self.high_rna = self.high_rna[self.high_rna['name'].isin(self.overlap_gene)]

        # summary the total rpf/reads and rna/reads
        self.high_rpf_num = self.high_rpf[self.rpf_sample].sum()
        self.high_rna_num = self.high_rna[self.rna_sample].sum()

        # summary the gene length and codon number
        self.gene_length_table = self.high_rpf.to_pandas().groupby('name').size()
        self.gene_codon_table = self.high_rpf.to_pandas().groupby(['name', 'codon']).size()


    def make_gene_table(self, rpf_sample=None, rna_sample=None):

        '''
        @Message  : make high expression gene table with the following parameters:.
        @Input    : high_rpf --> high expression gene density
                    rpf_sample --> sample name
                    rna_sample --> sample name
        @Return   : gene_table --> gene table contain the [codon count, gene length, codon number, rpf, rna, density, elongate rate, initiate rate]
        @Flow     : step1 --> create the class for each gene
                    step2 --> summary the [codon count, gene length, codon number, rpf, rna, density]
                    step3 --> calculate the density with the formula: density = rpf / rna
        '''

        # make the gene table
        self.gene_table[rpf_sample] = OrderedDict()
        for gene in self.overlap_gene:
            self.gene_table[rpf_sample][gene] = GeneTable(gene)
            self.gene_table[rpf_sample][gene].codon_count = dict(self.gene_codon_table[gene])
            self.gene_table[rpf_sample][gene].length = self.gene_length_table[gene] * 3
            self.gene_table[rpf_sample][gene].codon_num = self.gene_length_table[gene]
            self.gene_table[rpf_sample][gene].rpf = float(self.rpf_gene_sum[rpf_sample].loc[gene])
            self.gene_table[rpf_sample][gene].rna = float(self.rna_gene_sum[rna_sample].loc[gene])
            self.gene_table[rpf_sample][gene].density = self.gene_table[rpf_sample][gene].rpf / self.gene_table[rpf_sample][gene].rna


    def ribo_codon_frequency(self, rpf_sample=None):

        '''
        @Message  : calculate the mean rpf density of codon.
        @Input    : high_rpf --> high expression gene density
                    codon_dict --> codon message
                    rpf_sample --> sample name
        @Return   : codon_table --> codon message table contain the [rpf codon sum, rpf codon val, rpf sum, rpf proportion]
        @Flow     : step1 --> create the class for each codon
                    step2 --> summary the [rpf codon sum, rpf codon val, rpf sum, rpf proportion]
                                rpf_codon_proportion = Σ(codon) / Σ(gene count)
                    step3 --> merge the codon message
        '''
        
        for codon in self.codon_dict.keys():

            codon_flt = self.high_rpf.filter(pl.col("codon") == codon)
            # codon_flt = self.high_rpf[self.high_rpf["codon"] == codon]
            codon_sum = codon_flt[rpf_sample].sum()
            codon_proportion = float(codon_sum / self.high_rpf_num[rpf_sample][0])

            self.codon_table[rpf_sample][codon].rpf_codon_sum = len(codon_flt)
            self.codon_table[rpf_sample][codon].rpf_codon_val = (codon_flt[rpf_sample] > 0).sum()
            self.codon_table[rpf_sample][codon].rpf_sum = codon_sum
            self.codon_table[rpf_sample][codon].rpf_proportion = codon_proportion


    def mrna_codon_frequency(self, rpf_sample=None, rna_sample=None, times=0, initiate_rate=None):

        '''
        @Message  : calculate the mean rna density of codon.
        @Input    : high_rna --> high expression gene density
                    codon_dict --> codon message
                    rpf_sample --> sample name
                    rna_sample --> sample name
                    times --> iterative times
                    initiate_rate --> normalize the rna density with initiate_rate
        @Return   : codon_table --> codon message table contain the [rna codon sum, rna codon val, rna sum, rna proportion, cst]
        @Flow     : step1 --> create the class for each codon
                    step2 --> summary the [rna codon sum, rna codon val, rna sum, rna proportion] with the formula:
                                mrna_codon_proportion = Σ(codon) / Σ(gene count)
                    step3 --> calculate the cst with the formula: 
                                cst = rpf_codon_proportion / mrna_codon_proportion
                    step4 --> merge the codon message
        '''

        if initiate_rate:
            high_rna = self.high_rna.with_columns(
                # self.high_rna[rna_sample] * self.high_rna['name'].apply(lambda x: self.gene_table[rpf_sample][x].initiate_rate if x in self.gene_table[rpf_sample].keys() else 1))
                # Calling `map_elements` without specifying `return_dtype` can lead to unpredictable results. Specify `return_dtype` to silence this warning
                self.high_rna[rna_sample] * self.high_rna['name'].apply(
                    lambda x: self.gene_table[rpf_sample][x].initiate_rate if x in self.gene_table[rpf_sample].keys() else 1,
                    return_dtype=pl.Float32))
            
        else:
            high_rna = self.high_rna

        for codon in self.codon_dict.keys():
            codon_flt = high_rna.filter(pl.col("codon") == codon)
            # codon_flt = high_rna[high_rna["codon"] == codon]
            codon_sum = codon_flt[rna_sample].sum()
            self.codon_table[rpf_sample][codon].rna_codon_sum = len(codon_flt)
            self.codon_table[rpf_sample][codon].rna_codon_val = (codon_flt[rna_sample] > 0).sum()
            self.codon_table[rpf_sample][codon].rna_sum = codon_sum

            rna_proportion = float(codon_sum / self.high_rna_num[rna_sample][0])
            self.codon_table[rpf_sample][codon].rna_proportion.append(rna_proportion)

            try:
                cst = self.codon_table[rpf_sample][codon].rpf_proportion / self.codon_table[rpf_sample][codon].rna_proportion[times]
            except ZeroDivisionError:
                cst = 0

            self.codon_table[rpf_sample][codon].cst.append(cst)


    def codon_elongate_rate(self, rpf_sample=None, codon_num=None, codon_dict=None, times=0):

        '''
        @Message  : calculate the elongate rate for each codon of gene.
        @Input    : rpf_sample --> sample name
                    codon_num --> total codon number of each gene
                    codon_dict --> codon count
                    times --> iterative times
        @Return   : elongate_rate --> elongate rate for each codon of gene
        @Flow     : step1 --> calculate the codon density
                    step2 --> calculate the elongate rate with the formula: elongate_rate = codon_num / Σ(count * CST))
        '''
        
        codon_density = 0
        for codon, number in codon_dict.items():
            try:
                codon_density += codon_dict[codon] * self.codon_table[rpf_sample][codon].cst[times]
            except KeyError:
                pass

        elongate_rate = codon_num / codon_density
        return elongate_rate
    

    def compute_cst_iterative(self, rpf_sample=None, rna_sample=None):

        '''
        @Message  : compute cst iteratively.
        @Input    : rpf_sample --> sample name
                    rna_sample --> sample name
                    self.times --> iterative times
                    self.overlap_gene --> high expression gene name in the Ribo-seq and RNA-seq data
        @Return   : gene_table --> gene table contain the [elongate rate, initiate rate]
        @Flow     : step1 --> calculate the elongate rate for each codon of gene
                                formula: elongate_rate = codon_num / Σ(count * CST))
                    step2 --> calculate the initiate rate for each gene
                                formula: initiate_rate = elongate_rate * density
        '''
        
        for iters in range(0, self.times, 1):
            print('Times of iteration: {}'.format(iters + 1), flush=True)
            for gene in self.overlap_gene:
                self.gene_table[rpf_sample][gene].elongate_rate = self.codon_elongate_rate(
                    rpf_sample=rpf_sample,
                    codon_num=self.gene_table[rpf_sample][gene].codon_num,
                    codon_dict=self.gene_table[rpf_sample][gene].codon_count,
                    times=iters)
                
                self.gene_table[rpf_sample][gene].initiate_rate = self.gene_table[rpf_sample][gene].elongate_rate * self.gene_table[rpf_sample][gene].density

            # self.high_rna = self.high_rna.with_columns(
            #     self.high_rna[self.rna_sample] * self.high_rna['name'].apply(
            #     lambda x: self.gene_table[x].initiate_rate if x in self.gene_table.keys() else 1))

            self.mrna_codon_frequency(rpf_sample=rpf_sample, rna_sample=rna_sample, times=iters + 1, initiate_rate=True)
    

    def calc_cst(self):
        
        '''
        @Message  : calculate the mean reads density of codon.
        @Input    : high_rpf --> high expression gene density
                    high_rna --> high expression gene density
                    self.rpf_sample --> sample name
                    self.rna_sample --> sample name
        @Return   : codon_cst_table --> codon cst table contain the 
                                        [rpf codon sum, rpf codon val, rpf sum, rpf proportion]
                                        [rna codon sum, rna codon val, rna sum, rna proportion, cst]
        @Flow     : step1 --> make the codon message table
                    step2 --> calculate the mean rpf density of codon
                    step3 --> calculate the mean rna density of codon
                    step4 --> make the gene table
                    step5 --> calculate the elongate rate for each codon of gene
                    step6 --> compute cst iteratively
        '''
        
        for rpf_sample, rna_sample in zip(self.rpf_sample, self.rna_sample):
            print('\nCalculate the codon selection time for samples: {rpf} + {rna}'.format(rpf=rpf_sample, rna=rna_sample), flush=True)

            self.make_codon_table(rpf_sample=rpf_sample)
            self.filter_rpf_rna_expr(rpf_sample=rpf_sample, rna_sample=rna_sample)
            self.make_gene_table(rpf_sample=rpf_sample, rna_sample=rna_sample)

            self.ribo_codon_frequency(rpf_sample=rpf_sample)
            self.mrna_codon_frequency(rpf_sample=rpf_sample, rna_sample=rna_sample)

            if self.times > 0:
                print('Calculate the codon selection time iterative.', flush=True)
                self.compute_cst_iterative(rpf_sample=rpf_sample, rna_sample=rna_sample)
            else:
                print('Skip the iterative calculation of codon selection time.', flush=True)


    @staticmethod
    def scale_method(scale, cst):
        if scale == 'minmax':
            # relative_cst = (cst - cst.min()) / (cst.max() - cst.min())
            relative_cst = cst / cst.max()
        elif scale == 'zscore':
            relative_cst = (cst - cst.mean()) / cst.std()
        else:
            print('Unknown scale method.', flush=True)
            sys.exit()
        return relative_cst
    

    def format_cst_results(self):
        '''
        @Message  : format the codon selection time results.
        @Input    : self.codon_table --> codon message table
                    self.rpf_sample --> sample name
                    self.codon_dict --> codon message
        @Return   : self.detail_cst --> details of codon selection time results
                    self.absolute_cst --> absolute codon selection time results
                    self.relative_cst --> relative codon selection time results
                    self.merge_cst --> merged codon selection time results
        @Flow     : step1 --> merge all codon selection time
                    step2 --> retrieve the absolute and relative codon selection time
                    step3 --> retrieve the absolute and relative codon selection time
                    step4 --> merge the absolute and relative codon selection time
        '''
        
        # set the column name
        column_name = ['Codon', 'AA', 'Abbr', 'rpf_codon_val', 'rpf_sum', 'rpf_proportion',
                       'rna_codon_val', 'rna_sum', 'rna_proportion_first', 'rna_proportion_last']

        cst_name = ['absolute_cst_' + str(iters) for iters in range(self.times + 1)]
        relative_cst_name = ['relative_cst_' + str(iters) for iters in range(self.times + 1)]

        # merge all codon selection time
        detail_cst = OrderedDict()

        for rpf_sample in self.rpf_sample:
            detail_cst[rpf_sample] = OrderedDict()

            for codon in self.codon_dict.keys():
                detail_cst[rpf_sample][codon] = [self.codon_table[rpf_sample][codon].aa,
                                                 self.codon_table[rpf_sample][codon].abbr,
                                                 self.codon_table[rpf_sample][codon].rpf_codon_val,
                                                 self.codon_table[rpf_sample][codon].rpf_sum,
                                                 self.codon_table[rpf_sample][codon].rpf_proportion,
                                                 self.codon_table[rpf_sample][codon].rpf_codon_val,
                                                 self.codon_table[rpf_sample][codon].rna_sum,
                                                 self.codon_table[rpf_sample][codon].rna_proportion[0],
                                                 self.codon_table[rpf_sample][codon].rna_proportion[-1]
                                                ]
                detail_cst[rpf_sample][codon].extend(self.codon_table[rpf_sample][codon].cst)

            detail_cst_df = pd.DataFrame.from_dict(detail_cst[rpf_sample]).transpose().reset_index(drop=False, inplace=False)
            detail_cst_df.columns = column_name + cst_name

            detail_rel_cst_df = self.scale_method(self.scale, detail_cst_df.loc[:, cst_name])
            detail_rel_cst_df.columns = relative_cst_name

            detail_cst_df = pd.concat([detail_cst_df, detail_rel_cst_df], axis=1)
            detail_cst_df['sample'] = rpf_sample

            if self.detail_cst is None:
                self.detail_cst = detail_cst_df
            else:
                self.detail_cst = pd.concat([self.detail_cst, detail_cst_df], axis=0)

        # retrieve the absolute and relative codon selection time
        self.absolute_cst = self.detail_cst.pivot_table(index=['Codon', 'AA', 'Abbr'],
                                                        columns='sample',
                                                        values=cst_name[-1],
                                                        aggfunc='first').reset_index(drop=False, inplace=False)

        # retrieve the absolute and relative codon selection time
        self.relative_cst = self.detail_cst.pivot_table(index=['Codon', 'AA', 'Abbr'],
                                                        columns='sample',
                                                        values=relative_cst_name[-1],
                                                        aggfunc='first').reset_index(drop=False, inplace=False)

        abs_column_name = ['Codon', 'AA', 'Abbr'] + [i + '_absolute_cst' for i in self.rpf_sample]
        relative_column_name = ['Codon', 'AA', 'Abbr'] + [i + '_relative_cst' for i in self.rpf_sample]
        
        self.absolute_cst.columns = abs_column_name
        self.relative_cst.columns = relative_column_name

        self.merge_cst = pd.merge(left=self.absolute_cst, right=self.relative_cst, how='outer', on=['Codon', 'AA', 'Abbr'])
        self.merge_cst.sort_values(by=['Abbr', 'Codon'], ascending=True, inplace=True, ignore_index=True)


    def output_cst(self):
        # cst_writer = pd.ExcelWriter(self.output + '_absolute_cst.xlsx', engine='xlsxwriter')
        # with pd.ExcelWriter(self.output + '_codon_selection_time.xlsx',) as cst_writer:
        #     self.absolute_cst.to_excel(cst_writer, sheet_name='absolute_cst', index=False)
        #     self.relative_cst.to_excel(cst_writer, sheet_name='relative_cst', index=False)

        self.detail_cst.to_csv(self.output + '_iterative_codon_selection_time.txt', sep='\t', index=False)
        self.merge_cst.to_csv(self.output + '_codon_selection_time.txt', sep='\t', index=False)


    def draw_cst_line(self):
        import plotly.express as px
        out_pdf = self.output + "_cst_line.pdf"
        out_png = self.output + "_cst_line.png"
        
        matplotlib.use('AGG')
        rel_cst = self.relative_cst.melt(id_vars=['Codon', 'AA'],
                                         var_name='iteration',
                                         value_name='cst',
                                         value_vars=['cst_' + str(iters) for iters in range(self.times + 1)])

        fig = px.line(rel_cst, x = "iteration", y = "cst", title = 'codon selection time')
        # fig.show()
        fig.write_image(out_pdf)
        fig.write_image(out_png)


    def draw_cst_corr(self):
        out_pdf = self.output + "_cst_corrplot.pdf"
        out_png = self.output + "_cst_corrplot.png"

        cst_corr = self.merge_cst.filter(regex='_absolute_cst').astype('float').corr()
        cst_corr.to_csv(self.output + '_cst_corr.txt', sep='\t', index=True)

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6 + 0.2 * len(self.rpf_sample),
                                  8 + 0.2 * len(self.rpf_sample)),
                                  dpi=300)
        sns.heatmap(data=cst_corr, 
                    annot=True, 
                    fmt=".4f", 
                    cmap="vlag",
                    linewidths=.01)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()


    def draw_cst_heat(self):
        # output the absolute cst
        out_pdf = self.output + "_cst_heatplot.pdf"
        out_png = self.output + "_cst_heatplot.png"

        matplotlib.use('AGG')
        fig = plt.figure(figsize=(10 + 0.5 * len(self.rpf_sample), 20), dpi=300)
        sns.heatmap(data=self.merge_cst.filter(regex='_absolute_cst').astype('float'),
                    yticklabels=self.merge_cst['Codon'] + '[' + self.merge_cst['Abbr'].astype(str) + ']',
                    annot=True, 
                    fmt=".4f", 
                    cmap="vlag", 
                    linewidths=0)
        plt.tight_layout()
        # plt.show()

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)
        plt.close()
