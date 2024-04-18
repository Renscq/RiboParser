#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/6/11 10:50
# @Project : riboParser
# @Script  : Digestion.py

import os
import re
from collections import OrderedDict
from itertools import islice
from Bio import SeqIO
from .Ribo import Mrna

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pysam
import seaborn as sns

from Bio import motifs


class Ribo(object):

    def __init__(self, args):
        self.mrna_file = args.transcript
        self.mrna_seq = args.sequence
        self.mrna_dict = OrderedDict()
        self.longest = str(args.longest)

        self.bam = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]
        self.pysam_input = None

        self.min_length = args.min
        self.max_length = args.max
        self.profile = None

        self.rpf_start = None
        self.rpf_end = None

        self.align_left = []
        self.align_right = []
        self.seqlogo_5end = None
        self.seqlogo_3end = None

        self.scaled = args.scale

        self.output = args.output

    def read_transcript(self):

        record = SeqIO.parse(self.mrna_seq, "fasta")
        mrna_sequence = OrderedDict()
        for line in record:
            # sys.stdout.writelines("import gene:  {gene}\n".format(gene=line.id))
            now_sequence = line.upper()
            if re.findall(r'[^ATGC]', str(now_sequence.seq)):
                print(line.id)
                print(now_sequence.seq)
                continue
            mrna_sequence[line.id] = now_sequence.seq

        with open(self.mrna_file, 'r') as trans_file_in:
            if self.longest:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')

                    if record[9] == "True":
                        now_mrna = Mrna(record)
                        now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                        try:
                            now_mrna.seq = mrna_sequence[
                                now_mrna.transcript_id]
                            now_mrna.rpf = [0] * len(now_mrna.seq)
                            self.mrna_dict[now_mrna.transcript_id] = now_mrna
                        except KeyError:
                            continue
                    else:
                        continue
            else:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    now_mrna = Mrna(record)

                    now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                    try:
                        now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                        now_mrna.rpf = [0] * len(now_mrna.seq)
                        self.mrna_dict[now_mrna.transcript_id] = now_mrna
                    except KeyError:
                        continue

    def get_digest_sites(self):
        # import the bam file
        self.pysam_input = pysam.AlignmentFile(self.bam, "rb")
        # set the length range for the calculation
        if self.min_length < 20 and self.max_length > 100:
            rpf_start = {i: [0, 0, 0] for i in range(self.min_length, self.max_length + 1)}
            rpf_end = {i: [0, 0, 0] for i in range(self.min_length, self.max_length + 1)}
        elif self.min_length < 20:
            rpf_start = {i: [0, 0, 0] for i in range(self.min_length, 100 + 1)}
            rpf_end = {i: [0, 0, 0] for i in range(self.min_length, 100 + 1)}
        elif self.max_length > 100:
            rpf_start = {i: [0, 0, 0] for i in range(20, self.max_length + 1)}
            rpf_end = {i: [0, 0, 0] for i in range(20, self.max_length + 1)}
        else:
            rpf_start = {i: [0, 0, 0] for i in range(20, 100 + 1)}
            rpf_end = {i: [0, 0, 0] for i in range(20, 100 + 1)}

        for line in self.pysam_input.fetch(until_eof=True):
            if line.is_reverse:
                continue

            elif line.reference_name in self.mrna_dict:
                map_start, map_end = line.get_blocks()[0]
                read_length = line.infer_read_length()
                # check the location of rpf end
                mrna_length = self.mrna_dict[line.reference_name].length
                if map_start - 5 >= 0:
                    temp1 = self.mrna_dict[line.reference_name].seq[map_start - 5:map_start + 11]
                    if len(temp1) > 15:
                        self.align_left.append(temp1)
                if map_end + 5 < mrna_length:
                    # shift 1 nt for the alignment, (a,b] contain the map_end
                    temp2 = self.mrna_dict[line.reference_name].seq[map_end - 10:map_end + 6]
                    if len(temp2) > 15:
                        self.align_right.append(temp2)

                # filter the rpf with specific range
                if self.min_length <= read_length <= self.max_length:
                    frame_start = (map_start - self.mrna_dict[line.reference_name].utr5_length) % 3
                    frame_end = (map_end - 1 - self.mrna_dict[line.reference_name].utr5_length) % 3
                    rpf_start[read_length][frame_start] += 1
                    rpf_end[read_length][frame_end] += 1

        self.rpf_start = pd.DataFrame(rpf_start).T
        self.rpf_end = pd.DataFrame(rpf_end).T

    def output_digest_sites(self):
        psite_merge = pd.concat([self.rpf_start, self.rpf_end], join="outer", axis=1)
        psite_merge.fillna(0, inplace=True)
        psite_merge.columns = ["5end_frame1", "5end_frame2", "5end_frame3", "3end_frame1", "3end_frame2", "3end_frame3"]
        psite_merge.index.name = "length"
        psite_merge.to_csv(self.output + "_digestion_sites.txt", sep='\t', index=True)

    def digestion_plot(self):

        if self.scaled:
            out_pdf = self.output + "_scaled_digestion_sites_plot.pdf"
            out_png = self.output + "_scaled_digestion_sites_plot.png"

            psite_start_norm = self.rpf_start.sub(self.rpf_start.mean(axis=1), axis=0)
            psite_start_norm_row = psite_start_norm.div(psite_start_norm.std(axis=1), axis=0)
            psite_start_norm_row.fillna(0, inplace=True)

            psite_end_norm = self.rpf_end.sub(self.rpf_end.mean(axis=1), axis=0)
            psite_end_norm_row = psite_end_norm.div(psite_end_norm.std(axis=1), axis=0)
            psite_end_norm_row.fillna(0, inplace=True)

            now_cmap = "Blues"
        else:
            out_pdf = self.output + "_digestion_sites_plot.pdf"
            # out_png = self.output + "_digestion_sites_plot.png"

            psite_start_norm_row = self.rpf_start.copy()
            psite_end_norm_row = self.rpf_end.copy()

            now_cmap = "RdYlBu_r"

        codon_frame_5end = pd.DataFrame(self.rpf_start.loc[self.min_length:self.max_length].sum(axis=0)).T
        codon_frame_3end = pd.DataFrame(self.rpf_end.loc[self.min_length:self.max_length].sum(axis=0)).T
        rpfs_count_5end = pd.DataFrame(self.rpf_start.sum(axis=1)).loc[self.min_length:self.max_length].T
        rpfs_count_3end = pd.DataFrame(self.rpf_end.sum(axis=1)).loc[self.min_length:self.max_length].T

        matplotlib.use('Agg')
        # draw the periodicity of codon

        fig = plt.figure(figsize=(12, 10), dpi=300)
        ax1 = plt.subplot2grid((5, 6), (0, 0), rowspan=1, colspan=2)
        sns.barplot(data=codon_frame_5end, color="#388e70", alpha=.6, ax=ax1)
        ax1.set_xticklabels([])
        ax1.yaxis.set_ticks_position("right")
        ax1.set_ylabel("RPFs counts")
        ax1.set_title("5'-end")

        cbar_ax1 = fig.add_axes([0.43, 0.75, 0.01, 0.12])
        # ax11 = plt.subplot2grid((5, 6), (0, 2), rowspan=1, colspan=1)

        ax2 = plt.subplot2grid((5, 6), (1, 0), rowspan=4, colspan=2)
        sns.heatmap(
            data=psite_start_norm_row.loc[self.min_length:self.max_length],
            annot=None,
            linewidths=.5,
            ax=ax2,
            cmap=now_cmap,
            cbar_ax=cbar_ax1)
        label_y = ax2.get_yticklabels()
        plt.setp(label_y, rotation=0, horizontalalignment='right')
        ax2.set_ylabel("RPF length")
        ax2.set_xlabel("Codon frame")

        ax3 = plt.subplot2grid((5, 6), (1, 2), rowspan=4, colspan=1)
        sns.barplot(data=rpfs_count_5end,
                    color="#388e70",
                    alpha=.6,
                    ax=ax3,
                    orient='h')
        ax3.tick_params(bottom=True, top=False, left=False, right=False)
        ax3.set_yticklabels([])
        plt.xticks(rotation=60)
        ax3.set_xlabel("RPFs count")

        ax4 = plt.subplot2grid((5, 6), (0, 3), rowspan=1, colspan=2)
        sns.barplot(data=codon_frame_3end, color="#388e70", alpha=.6, ax=ax4)
        ax4.set_xticklabels([])
        ax4.yaxis.set_ticks_position("right")
        ax4.set_title("3'-end")

        cbar_ax2 = fig.add_axes([0.83, 0.75, 0.01, 0.12])
        ax5 = plt.subplot2grid((5, 6), (1, 3), rowspan=4, colspan=2)
        sns.heatmap(
            data=psite_end_norm_row.loc[self.min_length:self.max_length],
            annot=None,
            linewidths=.5,
            ax=ax5,
            cmap=now_cmap,
            cbar_ax=cbar_ax2)
        ax5.set_yticklabels([])
        ax5.set_xlabel("Codon frame")

        ax6 = plt.subplot2grid((5, 6), (1, 5), rowspan=4, colspan=1)
        sns.barplot(data=rpfs_count_3end,
                    color="#388e70",
                    alpha=.6,
                    ax=ax6,
                    orient='h')
        ax6.yaxis.set_ticks_position("right")
        plt.xticks(rotation=60)
        ax6.set_xlabel("RPFs count")

        # plt.tight_layout()
        plt.show()

        fig.savefig(fname=out_pdf)
        # fig.savefig(fname=out_png)

    def output_counts(self):
        self.seqlogo_5end = motifs.create(self.align_left)
        counts_5end = pd.DataFrame(self.seqlogo_5end.counts)
        pwm_5end = pd.DataFrame(self.seqlogo_5end.pwm)
        counts_5end.to_csv(self.output + "_5end_counts.txt", sep='\t', index=True)
        pwm_5end.to_csv(self.output + "_5end_pwm.txt", sep='\t', index=True)

        self.seqlogo_3end = motifs.create(self.align_right)
        counts_5end = pd.DataFrame(self.seqlogo_3end.counts)
        pwm_5end = pd.DataFrame(self.seqlogo_3end.pwm)
        counts_5end.to_csv(self.output + "_3end_counts.txt", sep='\t', index=True)
        pwm_5end.to_csv(self.output + "_3end_pwm.txt", sep='\t', index=True)

    def seq_logo_plot(self):
        import weblogo
        # get the 5'end digestion seqlogo
        seq_5end = weblogo.LogoData.from_counts(self.seqlogo_5end.alphabet,
                                                np.array(pd.DataFrame(self.seqlogo_5end.counts)))
        options = weblogo.LogoOptions()
        options.fineprint = '[-5, 10] from rpf 5-end digestion'
        format = weblogo.LogoFormat(seq_5end, options)
        pdf_5end = weblogo.pdf_formatter(seq_5end, format)
        with open(self.output + '_5end_seqlogo.pdf', 'wb') as f:
            f.write(pdf_5end)

        # get the 3'end digestion seqlogo
        seq_3end = weblogo.LogoData.from_counts(self.seqlogo_3end.alphabet,
                                                np.array(pd.DataFrame(self.seqlogo_3end.counts)))
        options = weblogo.LogoOptions()
        options.fineprint = '[-10, 5] from rpf 3-end digestion'
        format = weblogo.LogoFormat(seq_3end, options)
        pdf_3end = weblogo.pdf_formatter(seq_3end, format)
        with open(self.output + '_3end_seqlogo.pdf', 'wb') as f:
            f.write(pdf_3end)

    def seq_logo_plot2(self):
        import seqlogo
        # get the 5'end digestion seqlogo
        fname_5end = self.output + '_5end_seqlogo2.pdf'
        ppm_5end = seqlogo.Ppm(np.array([list(item) for item in self.seqlogo_5end.pwm.values()]))
        seqlogo.seqlogo(ppm_5end, ic_scale=False, format='pdf', size='medium', filename=fname_5end)
        seqlogo.seqlogo(ppm_5end, ic_scale=False, format='pdf', size='medium', filename=fname_5end)

        # get the 3'end digestion seqlogo
        fname_3end = self.output + '_3end_seqlogo2.pdf'
        ppm_3end = seqlogo.Ppm(np.array([list(item) for item in self.seqlogo_3end.pwm.values()]))
        seqlogo.seqlogo(ppm_3end, ic_scale=False, format='pdf', size='medium', filename=fname_3end)
        seqlogo.seqlogo(ppm_3end, ic_scale=False, format='pdf', size='medium', filename=fname_3end)
