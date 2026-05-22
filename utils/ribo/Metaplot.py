#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Metaplot.py


from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import polars as pl
import pandas as pd
import numpy as np


class Metaplot(object):

    def __init__(self, args):
        # input and output
        self.transcript = args.transcript
        self.gene = None
        self.rpf = args.rpf
        self.output = args.output

        # filter the range around tis/tts
        self.rpf_num = args.min
        self.utr5 = args.utr5
        self.cds = args.cds
        self.utr3 = args.utr3

        self.norm = args.normal
        self.sample_num = 0
        self.sample_name = None
        self.sample_dict = None

        self.tis_meta = None
        self.tts_meta = None
        self.merge_meta = pd.DataFrame()

        # for the high expression gene filter
        self.raw_rpf = None
        self.merged_rpf = None
        self.high_gene = None
        self.high_rpf = None
        self.gene_rpf_sum = None
        self.total_rpf_num = None

    def gene_anno(self):
        '''
        @Message  : function for gene annotation filter
        @Input    : self.utr5 --> input utr5 length
                    self.cds --> input cds length
                    self.utr3 --> input utr3 length
                    self.transcript --> input transcript file
        @Return   : self.gene --> output gene table
        @Flow     : step1 --> import the transcript file
                    step2 --> filter the transcript file with the input utr5/cds/utr3 length
                    step3 --> output the gene table
        '''
        
        anno = pd.read_csv(self.transcript, sep='\t', header=0, names=None)
        max_utr5 = anno['utr5_length'].max()
        max_cds = int(anno['cds_length'].max() / 2)
        max_utr3 = anno['utr3_length'].max()

        if self.utr5 > max_utr5:
            print('Specified UTR5 length more than max UTR5 length, median length instead.', flush=True)
            print('Specified UTR5 length: ' + str(self.utr5), flush=True)
            print('Max UTR5 length: ' + str(max_utr5), flush=True)
            print('Median UTR5 length: ' + str(self.utr5), flush=True)
            self.utr5 = int(anno['utr5_length'].median())
            anno_flt = anno.loc[anno['utr5_length'] >= self.utr5, ].copy()
        else:
            anno_flt = anno.loc[anno['utr5_length'] >= self.utr5, ].copy()

        if self.cds > max_cds:
            print('Specified CDS length more than half max CDS length, half median length instead.', flush=True)
            print('Specified CDS length: ' + str(self.cds), flush=True)
            print('Max CDS length: ' + str(max_utr5), flush=True)
            print('Median CDS length: ' + str(self.cds), flush=True)
            self.cds = int(anno['cds_length'].median() / 2)
            anno_flt = anno_flt.loc[anno['cds_length'] >= self.cds * 2, ].copy()
        else:
            anno_flt = anno_flt.loc[anno['cds_length'] >= self.cds * 2, ].copy()

        if self.utr3 > max_utr3:
            print('Specified UTR3 length more than max UTR3 length, median length instead.', flush=True)
            print('Specified UTR3 length: ' + str(self.utr3), flush=True)
            print('Max UTR3 length: ' + str(max_utr3), flush=True)
            print('Median UTR3 length: ' + str(self.utr3), flush=True)
            self.utr3 = int(anno['utr3_length'].median())
            anno_flt = anno_flt.loc[anno['utr3_length'] >= self.utr3, ].copy()
        else:
            anno_flt = anno_flt.loc[anno['utr3_length'] >= self.utr3, ].copy()

        self.gene = anno_flt['transcript_id']

    def output_meta(self, now_sp):
        out_txt = now_sp + "_tis_tts_metaplot.txt"

        # output the tis meta
        self.tis_meta = pd.melt(self.tis_meta, id_vars='Codon', var_name='Sample', value_name='Density')
        self.tis_meta['Frame'] = self.tis_meta['Sample'].str[-1].astype(int)
        self.tis_meta["Nucleotide"] = self.tis_meta["Codon"] * 3 + self.tis_meta["Frame"]
        self.tis_meta.sort_values(by=['Nucleotide', 'Frame'], inplace=True)
        self.tis_meta['Sample'] = now_sp
        self.tis_meta['Meta'] = "TIS"
        self.tis_meta = self.tis_meta.loc[:, ['Sample', 'Meta', 'Nucleotide', 'Codon', 'Frame', 'Density']]

        # output the tts meta
        self.tts_meta = pd.melt(self.tts_meta, id_vars='Codon', var_name='Sample', value_name='Density')
        self.tts_meta['Frame'] = self.tts_meta['Sample'].str[-1].astype(int)
        self.tts_meta["Nucleotide"] = self.tts_meta["Codon"] * 3 + self.tts_meta["Frame"]
        self.tts_meta.sort_values(by=['Nucleotide', 'Frame'], inplace=True)
        self.tts_meta['Sample'] = now_sp
        self.tts_meta['Meta'] = "TTS"
        self.tts_meta = self.tts_meta.loc[:, ['Sample', 'Meta', 'Nucleotide', 'Codon', 'Frame', 'Density']]

        tis_tts_meta = pd.concat([self.tis_meta, self.tts_meta], axis=0)
        tis_tts_meta.to_csv(out_txt, sep='\t', index=False)


    def output_merge_meta(self):
        # output the merge meta
        out_txt = self.output + "_tis_tts_metaplot.txt"
        self.merge_meta.to_csv(out_txt, sep='\t', index=False)


    def read_rpf(self):
        '''
        @Message  : function for read rpf file
        @Input    : self.rpf --> input rpf file
                    
        @Return   : output --> description
        @Flow     : step1 --> import the rpf file
                    step2 --> filter the rpf file with the input gene id
                    step3 --> output the rpf file
        '''
        
        self.raw_rpf = pl.read_csv(self.rpf, has_header = True, separator = '\t')
        self.raw_rpf = self.raw_rpf.to_pandas()

        rpf_mer = self.raw_rpf.loc[self.raw_rpf['name'].isin(list(self.gene)), ]
        self.sample_num = int((len(rpf_mer.columns) - 6) / 3)
        self.sample_name = pd.Series(rpf_mer.columns[6::]).str[:-3].drop_duplicates().tolist()
        self.sample_dict = OrderedDict()
        rpf_idx = rpf_mer.iloc[:, 0:6].copy()

        for now_sp in self.sample_name:
            print('Import the RPFs file {file_name}.'.format(file_name=now_sp), flush=True)
            now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
            now_sp_rpf = rpf_mer.loc[:, now_sp_index].copy()
            now_rpf = pd.concat([rpf_idx, now_sp_rpf], axis=1)

            # get the high expression gene id
            now_rpf["sum"] = pd.DataFrame(now_sp_rpf.sum(axis=1))
            cds_rpf = now_rpf[now_rpf["region"] == "cds"][["name", "sum"]].groupby("name")["sum"].apply(np.sum)
            self.high_gene = cds_rpf[cds_rpf > self.rpf_num]

            # filter the high expression genes
            if self.norm:
                total_rpf = now_sp_rpf.sum().sum()
                now_sp_rpm = now_sp_rpf.div(total_rpf) * 1e7
                now_rpm = pd.concat([rpf_idx, now_sp_rpm], axis=1)
                self.high_rpf = now_rpm[now_rpm["name"].isin(self.high_gene.index)]
            else:
                self.high_rpf = now_rpf[now_rpf["name"].isin(self.high_gene.index)]

            tis_rpf = self.high_rpf[self.high_rpf["from_tis"].isin(range(-self.utr5, self.cds))]
            tts_rpf = self.high_rpf[self.high_rpf["from_tts"].isin(range(-self.cds + 1, self.utr3 + 1))]
            self.tis_meta = tis_rpf.groupby("from_tis")[now_sp_index].apply(sum).div(len(self.high_gene)).reset_index()
            self.tts_meta = tts_rpf.groupby("from_tts")[now_sp_index].apply(sum).div(len(self.high_gene)).reset_index()

            self.tis_meta.rename(columns={"from_tis": "Codon"}, inplace=True)
            self.tts_meta.rename(columns={"from_tts": "Codon"}, inplace=True)

            self.output_meta(now_sp)
            self.sample_dict[now_sp] = [self.tis_meta, self.tts_meta, len(self.high_gene), now_sp_index]

            # merge all tis and tts meta
            if self.merge_meta is None:
                self.merge_meta = pd.concat([self.tis_meta, self.tts_meta], axis=0)
            else:
                self.merge_meta = pd.concat([self.merge_meta, self.tis_meta, self.tts_meta], axis=0)


    def draw_bar_metaplot(self):
        for sample_name, mess in self.sample_dict.items():
            tis_meta, tts_meta, gene_num, sample_frame = mess[0], mess[1], mess[2], mess[3]
            
            if all(value == 0 for value in tis_meta['Density']):
                print("Warning: All Density values are 0.")

            out_pdf = sample_name + "_meta_bar_plot.pdf"
            out_png = sample_name + "_meta_bar_plot.png"

            matplotlib.use('AGG')
            colors = {0: '#66c2a5', 1: '#fc8d62', 2: '#8da0cb'}

            fig = plt.figure(figsize=(15, 5), dpi=300)
            # draw the figure around the TIS region
            ax1 = fig.add_subplot(121)
            for key, group in tis_meta.groupby('Frame'):
                ax1.bar(group['Nucleotide'], group['Density'], color=colors[key], width=0.6, label=f'Frame {key}')

            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax1.set_title("Mean RPFs of TIS region ({number} genes)".format(number=gene_num), fontsize=14)
            ax1.set_ylabel('Mean RPFs', fontsize=12)
            ax1.set_xlabel('From start codon (nt)', fontsize=12)
            ax1.set_xticks([-self.utr5 * 3, 0, self.cds * 3 - 1])
            # ax1.set_xticklabels([-utr5 * 3, 0, cds * 3], fontsize=20)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            ax1.legend(fontsize=12)

            # draw the figure around the tts region
            ax2 = fig.add_subplot(122)
            for key, group in tts_meta.groupby('Frame'):
                ax2.bar(group['Nucleotide'], group['Density'], color=colors[key], width=0.6, label=f'Frame {key}')

            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax2.set_title("Mean RPFs of TTS region ({number} genes)".format(number=gene_num), fontsize=14)
            ax2.set_ylabel('Mean RPFs', fontsize=12)
            ax2.set_xlabel('From stop codon (nt)', fontsize=12)
            ax2.set_xticks([-self.cds * 3 + 1, 0, self.utr3 * 3])
            # ax2.set_xticklabels([-cds * 3 + 1, 0, utr3 * 3], fontsize=20)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            ax2.legend(fontsize=12)

            fig.tight_layout()
            # plt.show()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()

    def draw_line_metaplot(self):
        for sample_name, mess in self.sample_dict.items():
            tis_meta, tts_meta, gene_num, sample_frame = mess[0], mess[1], mess[2], mess[3]

            out_pdf = sample_name + "_meta_line_plot.pdf"
            out_png = sample_name + "_meta_line_plot.png"

            matplotlib.use('AGG')
            fig = plt.figure(figsize=(15, 5), dpi=300)
            # draw the figure around the TIS region
            ax1 = fig.add_subplot(121)
            tis_x = tis_meta['Nucleotide']
            ax1.plot(tis_x, tis_meta["Density"])
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax1.set_title("Mean RPFs of TIS region ({number} genes)".format(number=gene_num), fontsize=14)
            ax1.set_ylabel('Mean RPFs', fontsize=12)
            ax1.set_xlabel('from start codon (nt)', fontsize=12)
            ax1.set_xticks([-self.utr5 * 3, 0, self.cds * 3 - 1])
            # ax1.set_xticklabels([-utr5 * 3, 0, cds * 3], fontsize=20)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

            # draw the figure around the tts region
            ax2 = fig.add_subplot(122)
            tts_x = tts_meta['Nucleotide']
            ax2.plot(tts_x, tts_meta["Density"])
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax2.set_title("Mean RPFs of TTS region ({number} genes)".format(number=gene_num), fontsize=14)
            ax2.set_ylabel('Mean RPFs', fontsize=12)
            ax2.set_xlabel('from stop codon (nt)', fontsize=12)
            ax2.set_xticks([-self.cds * 3 + 1, 0, self.utr3 * 3])
            # ax2.set_xticklabels([-cds * 3 + 1, 0, utr3 * 3], fontsize=20)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

            fig.tight_layout()
            # plt.show()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()
