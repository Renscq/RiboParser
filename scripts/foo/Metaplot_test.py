#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Metaplot.py


from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
from foo import RPFs


class Metaplot(object):

    def __init__(self, args):
        # input and output
        self.gene = args.list
        self.rpf = args.rpf

        # filter the range around tis/tts
        self.site = 'P' # 'E', 'P' or 'A'
        self.frame =  'all' # 'all' or 0, 1, 2

        self.rpf_num = args.min
        self.utr5 = args.utr5
        self.cds = args.cds
        self.utr3 = args.utr3

        self.tis_rpf = None
        self.tts_rpf = None

        self.tis_meta = None
        self.tts_meta = None

        # for the high expression gene filter
        self.norm = args.normal
        self.raw_rpf = None

        self.sample_num = 0
        self.sample_name = None

        self.high_gene = None
        self.high_rpf = None

        self.total_rpf_num = None

    def read_rpf(self):
        '''
        @Message  : function for import the rpf file.
        @Input    : self.rpf --> rpf file
                    self.gene --> rpf gene list
                    self.rpf_num --> minimum rpf number

        @Return   : self.raw_rpf --> rpf data frame
                    self.sample_name --> sample name
                    self.total_rpf_num --> total rpf number
                    self.gene_rpf_sum --> gene rpf sum
                    self.high_gene --> high expression gene
                    self.high_rpf --> high expression rpf

        @Flow     : step1 --> import the rpf file
                    step2 --> filter the high expression gene
                    step3 --> filter the length of utr5/utr3/cds
                    step4 --> normalize the rpf

        '''
        
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf,
                                      sites=self.site,
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis=None,
                                      tts=None,
                                      gene=self.gene,
                                      rpf_num=self.rpf_num)
        
        self.raw_rpf = rpf_results[0].to_pandas()
        self.sample_name = rpf_results[1].copy()
        # self.sample_num = rpf_results[2]
        # self.merged_rpf = rpf_results[3]
        self.total_rpf_num = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        self.high_gene = rpf_results[6]
        # self.high_rpf = rpf_results[7].copy()
        del rpf_results

        # filter the high expression genes
        # fitted_gene = set(self.gene_fit['transcript_id']) & set(self.high_gene)
        self.high_rpf = self.raw_rpf[self.raw_rpf["name"].isin(self.high_gene)]

        # normalize the rpf
        self.total_rpf_num_frame = self.total_rpf_num.repeat(3)
        self.total_rpf_num_frame.index = self.high_rpf.iloc[:, 6::].columns

        if self.norm:
            self.high_rpf.iloc[:, 6::] = self.high_rpf.iloc[:, 6::] * 1e7 / self.total_rpf_num_frame

    def length_filter(self):
        '''
        @Message  : function for .
        @Input    : param --> description
        @Return   : output --> description
        @Flow     : step1 --> run
        '''

        # filter the gene utr5/cds/utr3 length
        gene_length = self.high_rpf.groupby(["name", "region"]).apply(len).reset_index().rename(columns={0: "length"})
        gene_length = gene_length.pivot(index="name", columns="region", values="length").reset_index()
        gene_length = gene_length[(gene_length["3utr"] >= self.utr3) & 
                                  (gene_length["5utr"] >= self.utr5) & 
                                  (gene_length["cds"] >= self.cds)]

        # filter the length of utr5/utr3/cds
        self.high_rpf = self.high_rpf[self.high_rpf["name"].isin(gene_length["name"])]

        self.tis_rpf = self.high_rpf[self.high_rpf["from_tis"].isin(range(-self.utr5, self.cds))]
        self.tts_rpf = self.high_rpf[self.high_rpf["from_tts"].isin(range(-self.cds + 1, self.utr3 + 1))]

    def read_rpf2(self):
        '''
        @Message  : function for read rpf file
        @Input    : self.rpf --> input rpf file
                    
        @Return   : output --> description
        @Flow     : step1 --> run
        '''
        
        self.raw_rpf = pd.read_csv(self.rpf, sep='\t', header=0, names=None)
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
            cds_rpf = now_rpf[now_rpf["region"] == "cds"][["name", "sum"]].groupby("name")["sum"].apply(sum)
            self.high_gene = cds_rpf[cds_rpf > self.rpf_num]

            # filter the high expression genes
            if self.norm:
                total_rpf = now_sp_rpf.sum().sum()
                now_sp_rpm = now_sp_rpf.div(total_rpf) * 1e8
                now_rpm = pd.concat([rpf_idx, now_sp_rpm], axis=1)
                self.high_rpf = now_rpm[now_rpm["name"].isin(self.high_gene.index)]
            else:
                self.high_rpf = now_rpf[now_rpf["name"].isin(self.high_gene.index)]

            tis_rpf = self.high_rpf[self.high_rpf["from_tis"].isin(range(-self.utr5, self.cds))]
            tts_rpf = self.high_rpf[self.high_rpf["from_tts"].isin(range(-self.cds + 1, self.utr3 + 1))]
            tis_meta = tis_rpf.groupby("from_tis")[now_sp_index].apply(sum).div(len(self.high_gene))
            tts_meta = tts_rpf.groupby("from_tts")[now_sp_index].apply(sum).div(len(self.high_gene))

            self.sample_dict[now_sp] = [tis_meta, tts_meta, len(self.high_gene), now_sp_index]
            self.output_meta(now_sp_index, tis_meta, tts_meta, now_sp)

    def output_meta(self):

        out_txt = now_sp + "_metaplot.txt"
        tis_meta["from_tis"] = tis_meta.index
        tts_meta["from_tts"] = tts_meta.index
        tis_meta.index = range(len(tis_meta))
        tts_meta.index = range(len(tts_meta))
        tis_meta = tis_meta[["from_tis"] + now_sp_index]
        tts_meta = tts_meta[["from_tts"] + now_sp_index]
        meta_data = pd.concat([tis_meta, tts_meta], axis=1)
        meta_data.to_csv(out_txt, sep='\t', index=False)



    def draw_bar_metaplot(self):
        for sample_name, mess in self.sample_dict.items():
            tis_meta, tts_meta, gene_num, sample_frame = mess[0], mess[1], mess[2], mess[3]

            out_pdf = sample_name + "_meta_bar_plot.pdf"
            out_png = sample_name + "_meta_bar_plot.png"

            matplotlib.use('AGG')
            fig = plt.figure(figsize=(12, 8), dpi=300)
            # draw the figure around the TIS region
            ax1 = fig.add_subplot(211)
            tis_x = tis_meta.from_tis
            width = 0.3
            ax1.bar(tis_x - 0.33, tis_meta[sample_frame[0]], width, label='frame0')
            ax1.bar(tis_x, tis_meta[sample_frame[1]], width, label='frame1')
            ax1.bar(tis_x + 0.33, tis_meta[sample_frame[2]], width, label='frame2')
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax1.set_title("Mean RPFs of TIS region ({number} genes)".format(number=gene_num), fontsize=20)
            ax1.set_ylabel('Mean RPFs', fontsize=18)
            ax1.set_xlabel('from start codon (AA)', fontsize=18)
            ax1.set_xticks([-self.utr5, 0, self.cds - 1])
            # ax1.set_xticklabels([-utr5 * 3, 0, cds * 3], fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            ax1.legend(fontsize=15)

            # draw the figure around the tts region
            ax2 = fig.add_subplot(212)
            tts_x = tts_meta.from_tts
            width = 0.3
            ax2.bar(tts_x - 0.33, tts_meta[sample_frame[0]], width, label='frame0')
            ax2.bar(tts_x, tts_meta[sample_frame[1]], width, label='frame1')
            ax2.bar(tts_x + 0.33, tts_meta[sample_frame[2]], width, label='frame2')
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax2.set_title("Mean RPFs of TTS region ({number} genes)".format(number=gene_num), fontsize=20)
            ax2.set_ylabel('Mean RPFs', fontsize=18)
            ax2.set_xlabel('from stop codon (AA)', fontsize=18)
            ax2.set_xticks([-self.cds + 1, 0, self.utr3])
            # ax2.set_xticklabels([-cds * 3 + 1, 0, utr3 * 3], fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            ax2.legend(fontsize=15)

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

            tis_meta = pd.melt(tis_meta, id_vars='from_tis').sort_values(by=['from_tis', 'variable'])
            tis_meta.index = list(range(-self.utr5 * 3, self.cds * 3))
            tts_meta = pd.melt(tts_meta, id_vars='from_tts').sort_values(by=['from_tts', 'variable'])
            tts_meta.index = list(range(-self.cds * 3 + 1, self.utr3 * 3 + 1))

            matplotlib.use('AGG')
            fig = plt.figure(figsize=(12, 8), dpi=300)
            # draw the figure around the TIS region
            ax1 = fig.add_subplot(211)
            tis_x = tis_meta.index
            ax1.plot(tis_x, tis_meta["value"])
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax1.set_title("Mean RPFs of TIS region ({number} genes)".format(number=gene_num), fontsize=20)
            ax1.set_ylabel('Mean RPFs', fontsize=18)
            ax1.set_xlabel('from start codon (nt)', fontsize=18)
            ax1.set_xticks([-self.utr5 * 3, 0, self.cds * 3 - 1])
            # ax1.set_xticklabels([-utr5 * 3, 0, cds * 3], fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            # draw the figure around the tts region
            ax2 = fig.add_subplot(212)
            tts_x = tts_meta.index
            ax2.plot(tts_x, tts_meta["value"])
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax2.set_title("Mean RPFs of TTS region ({number} genes)".format(number=gene_num), fontsize=20)
            ax2.set_ylabel('Mean RPFs', fontsize=18)
            ax2.set_xlabel('from stop codon (nt)', fontsize=18)
            ax2.set_xticks([-self.cds * 3 + 1, 0, self.utr3 * 3])
            # ax2.set_xticklabels([-cds * 3 + 1, 0, utr3 * 3], fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            fig.tight_layout()
            # plt.show()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)
            plt.close()
