#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Offset_SSCBM.py


import os.path
import sys
from collections import OrderedDict
from itertools import islice

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


class Mrna(object):

    def __init__(self, record):
        self.chromosome = record[0]
        self.gene_id = record[1]
        #self.gene_name = record[2]
        self.transcript_id = record[2]
        self.start = record[3]
        self.end = record[4]
        self.utr5_length = int(record[5])
        self.cds_length = int(record[6])
        self.utr3_length = int(record[7])
        self.strand = record[8]
        self.rep_transcript = record[9]
        self.modified = record[10]
        self.bam = []
        self.rpf = []
        self.seq = None


class Offset(object):

    def __init__(self, args):
        self.min_length = args.min
        self.max_length = args.max
        self.nt_num = self.max_length - self.min_length + 1
        self.mrna_file = args.transcript
        
        # set the offset to contain the raw reads offset, adj_tis_offset to alter the error offset.
        # self.column_name = ["length", "from_tis", "tis_5end", "tis_5end_per", "from_tts", "tts_3end", "tts_3end_per", "sum"]

        # self.tis_offset = {"start_codon": OrderedDict(), "stop_codon": OrderedDict()}
        self.align = args.align
        self.tis_offset = {"tis_5end": OrderedDict(), "tis_3end": OrderedDict(),
                           "tts_5end": OrderedDict(), "tts_3end": OrderedDict()}
        self.adj_tis_offset = OrderedDict()
        self.tis_5end = None
        self.tis_3end = None
        self.tts_5end = None
        self.tts_3end = None
        self.detail = args.detail

        self.frame_offset = OrderedDict()
        self.merge_frame_offset = pd.DataFrame(columns=["length", "frame0", "rpfs0", "frame1", "rpfs1",
                                                        "frame2", "rpfs2", "p_site", "rpfs"])
        self.frame_offset_len = {}
        # self.mode = args.mode

        self.exp_offset = None
        self.shift_nt = args.shift

        self.mrna_dict = OrderedDict()
        self.length_dict = {}

        # arguments for bam file parsing
        self.silence = args.silence
        self.sample_file = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]

        self.pysam_input = None
        self.pysam_output = None
        self.output_prefix = None

        self.bam_file_list = None
        self.sample_dict = None

        self.peak_length = args.exp_peak
        self.peak_reads = 0

        self.output_prefix = args.output

    def read_transcript(self):
        # read transcripts annotation, retrieve the TIS and TTS message for offset detection.
        trans_file = self.mrna_file

        with open(trans_file, 'r') as trans_file_in:
            for line in islice(trans_file_in, 1, None):
                record = line.strip().split('\t')
                now_mrna = Mrna(record)
                self.mrna_dict[now_mrna.transcript_id] = now_mrna

    def get_mrna_reads(self):

        # check the data type, to support the SAM/BAM file format both.
        if self.sample_format.lower() == ".bam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            read_format = 'rb'
        elif self.sample_format.lower() == ".sam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            read_format = 'r'
        else:
            print("Unknown file format, please input the correct bam or sam file.", flush=True)
            sys.exit()

        self.pysam_input = pysam.AlignmentFile(self.sample_file, read_format)

        # calculate the length distribution and eliminate the reads aligned to the negative strand.
        for line in self.pysam_input.fetch(until_eof=True):
            if line.reference_name in self.mrna_dict:
                read_length = line.infer_read_length()
                if line.is_reverse:
                    try:
                        self.length_dict[read_length][1] += 1
                    except KeyError:
                        self.length_dict[read_length] = [0, 1]
                else:
                    try:
                        self.length_dict[read_length][0] += 1
                    except KeyError:
                        self.length_dict[read_length] = [1, 0]
                    self.mrna_dict[line.reference_name].bam.append(line)
            elif self.silence:
                pass
            else:
                print(line.reference_name + ': not in reference file.', flush=True)

        peak_length, self.peak_reads = max(self.length_dict.items(), key=lambda length: length[1])
        
        # if self.peak_length:
        #     self.peak_length -= 1
        # else:
        #     self.peak_length = peak_length - 1

        self.pysam_input.close()

    def get_tis_offset(self):

        # build the offset table
        for number in range(self.min_length, self.max_length + 1):
            # self.tis_offset["tis_5end"][number] = OrderedDict({-i: 0 for i in range(self.max_length + 1, 0, -1)})
            self.tis_offset["tis_5end"][number] = OrderedDict({-i: 0 for i in range(self.max_length, -1, -1)})
            self.tis_offset["tis_3end"][number] = OrderedDict({i: 0 for i in range(self.max_length + 1)})
            # self.tis_offset["tts_5end"][number] = OrderedDict({-i: 0 for i in range(self.max_length + 1, 0, -1)})
            self.tis_offset["tts_5end"][number] = OrderedDict({-i: 0 for i in range(self.max_length, -1, -1)})
            self.tis_offset["tts_3end"][number] = OrderedDict({i: 0 for i in range(self.max_length + 1)})

        for mrna, mrna_attr in self.mrna_dict.items():
            # eliminate the transcript without read mapped
            if len(mrna_attr.bam) == 0:
                continue
            else:
                # get the TIS and TTS from mRNA dict.
                # 6:shift 3 nt for stop codon, shift 2 nt for p-site， and shift 1 for the python index start
                # cds_start, cds_end = mrna_attr.utr5_length, mrna_attr.utr5_length + mrna_attr.cds_length - 6
                cds_start, cds_end = mrna_attr.utr5_length, mrna_attr.utr5_length + mrna_attr.cds_length - 6

                for line in mrna_attr.bam:
                    map_start, map_end = line.get_blocks()[0]
                    # map start is not included, so add 1 nt.
                    # map_start = map_start + 1
                    map_start = map_start
                    read_length = line.infer_read_length()

                    # reads that match a specific length are retained.
                    # an other way is to calculate whole rpfs, and filter the min-max length in the plot module
                    if self.min_length <= read_length <= self.max_length:
                        # filter the possible range of offset lengths around the start codon.
                        if map_start <= cds_start <= map_end:
                            offset_5end = map_start - cds_start
                            offset_3end = map_end - cds_start
                            try:
                                self.tis_offset["tis_5end"][read_length][offset_5end] += 1
                            except KeyError:
                                # self.tis_offset["tis_5end"][read_length][offset_5end] = 1
                                pass
                            try:
                                self.tis_offset["tis_3end"][read_length][offset_3end] += 1
                            except KeyError:
                                pass
                                # self.tis_offset["tis_3end"][read_length][offset_3end] = 1

                        # filter the possible range of offset lengths around the stop codon.
                        elif map_start <= cds_end <= map_end:
                            offset_5end = map_start - cds_end
                            offset_3end = map_end - cds_end
                            try:
                                self.tis_offset["tts_5end"][read_length][offset_5end] += 1
                            except KeyError:
                                pass
                                # self.tis_offset["tts_5end"][read_length][offset_5end] = 1
                            try:
                                self.tis_offset["tts_3end"][read_length][offset_3end] += 1
                            except KeyError:
                                pass
                                # self.tis_offset["tts_3end"][read_length][offset_3end] = 1
                    else:
                        pass

        self.tis_5end = pd.DataFrame(self.tis_offset['tis_5end']).T
        self.tis_3end = pd.DataFrame(self.tis_offset['tis_3end']).T
        self.tts_5end = pd.DataFrame(self.tis_offset['tts_5end']).T
        self.tts_3end = pd.DataFrame(self.tis_offset['tts_3end']).T

        if self.detail:
            self.tis_5end.to_csv(self.output_prefix + '_' + 'tis_5end.txt', sep='\t', header=True, index=True)
            self.tis_3end.to_csv(self.output_prefix + '_' + 'tis_3end.txt', sep='\t', header=True, index=True)
            self.tts_5end.to_csv(self.output_prefix + '_' + 'tts_5end.txt', sep='\t', header=True, index=True)
            self.tts_3end.to_csv(self.output_prefix + '_' + 'tts_3end.txt', sep='\t', header=True, index=True)


    @staticmethod
    def offset_scale(offset):
        offset_norm = offset.sub(offset.mean(axis=1), axis=0)
        offset_norm_row = offset_norm.div(offset_norm.std(axis=1), axis=0)
        offset_norm_row.fillna(0, inplace=True)
        return offset_norm_row

    def shift_codon(self, now_psite, length, left_offset, right_offset):

        # shift the codon with 3nt to the normal p-site range
        if now_psite in range(left_offset, right_offset):
            return now_psite

        if self.adj_tis_offset:
            if self.adj_tis_offset[length - 1][-1]:
                if self.adj_tis_offset[length - 1][-1] + 3 < now_psite:
                    return self.shift_codon(now_psite - 3, length, left_offset, right_offset)
                
                elif now_psite < self.adj_tis_offset[length - 1][-1] - 3:
                    return self.shift_codon(now_psite + 3, length, left_offset, right_offset)
            else:
                return np.nan

        elif not self.adj_tis_offset:
            if now_psite > right_offset:
                return self.shift_codon(now_psite - 3, length, left_offset, right_offset)
            
            elif now_psite < left_offset:
                return self.shift_codon(now_psite + 3, length, left_offset, right_offset)

    def adjust_tis_offset(self):
        '''
        @Message  : function for .
        @Input    : param --> description
        @Return   : output --> description
        @Flow     : step1 --> merge the tis_5end and tts_5end to get the offset rpfs
                    step2 --> shift the offset to the normal p-site range
        '''
        
        # monosome data only needs to align p-site on one side.
        tis_5end = self.tis_5end.T.copy()
        tis_5end.index = abs(tis_5end.index)

        tts_5end = self.tts_5end.T.copy()
        tts_5end.index = abs(tts_5end.index)
        
        if self.align == "tis":
            offset_rpfs = tis_5end
        elif self.align == "tts":
            offset_rpfs = tts_5end
        elif self.align == "both":
            offset_rpfs = tis_5end + tts_5end
        else:
            print("Please input the correct align mode: tis, tts, both.", flush=True)
            print("Use the default mode: both.", flush=True)
            offset_rpfs = tis_5end + tts_5end

        # adjust the offset to the normal p-site range
        for length in range(self.min_length, self.max_length + 1):
            # shift 1 nt for the both 5/3 reads end (5end + 3end = 2)
            shift_nt = (length - self.peak_length) // self.shift_nt

            left_offset = int(np.floor(self.peak_length / 3)) + shift_nt
            right_offset = int(np.ceil(self.peak_length * 3 / 4)) + shift_nt + 2
            
            # candidate_offset = [i for i in range(8 + shift_nt, 16 + shift_nt + 2)]
            candidate_offset = [i for i in range(left_offset, right_offset)]

            # compare the max rpfs at same offset of 5end and 3end
            max_offset_site = offset_rpfs[length].idxmax()
            temp_offset_rpfs = offset_rpfs[length].copy()
            
            if temp_offset_rpfs[max_offset_site] == 0:
                max_offset_site = candidate_offset[5]
                
            while max_offset_site not in candidate_offset:
                del temp_offset_rpfs[max_offset_site]
                max_offset_site = temp_offset_rpfs.idxmax()

            # sum the max frame rpfs of [E-site, P-site, A-site], and the psite is the max_offset_site
            offset_frame0_range = range(max_offset_site - 3, max_offset_site + 6, 3)
            offset_frame1_range = range(max_offset_site - 2, max_offset_site + 6, 3)
            offset_frame2_range = range(max_offset_site - 1, max_offset_site + 6, 3)

            frame0_offset_rpfs = offset_rpfs[length][offset_frame0_range].sum()
            frame1_offset_rpfs = offset_rpfs[length][offset_frame1_range].sum()
            frame2_offset_rpfs = offset_rpfs[length][offset_frame2_range].sum()

            max_frame_offset_rpfs = max(frame0_offset_rpfs, frame1_offset_rpfs, frame2_offset_rpfs)

            # sum the all frame rpfs of [E-site, P-site, A-site]
            offset_sum = frame0_offset_rpfs + frame1_offset_rpfs + frame2_offset_rpfs
            max_offset_per = np.true_divide(max_frame_offset_rpfs, offset_sum)

            # psite = offset + 1
            psite = max_offset_site + 1
            psite = self.shift_codon(psite, length, left_offset, right_offset + 1)

            self.adj_tis_offset[length] = [length,
                                           max_offset_site, frame0_offset_rpfs,
                                           max_offset_site + 1, frame1_offset_rpfs,
                                           max_offset_site + 2, frame2_offset_rpfs,
                                           psite, offset_sum, max_offset_per]


    def write_tis_offset(self):

        column_name = ["length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2", "p_site", "rpfs", "periodicity"]

        adj_tis_offset = pd.DataFrame(self.adj_tis_offset).T.copy()

        adj_tis_offset.columns = column_name
        adj_tis_offset["periodicity"] = adj_tis_offset["periodicity"] * 100
        adj_tis_offset["periodicity"] = adj_tis_offset["periodicity"].astype(float).round(2)

        adj_tis_offset[["length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2", "p_site",
                        "rpfs"]] = adj_tis_offset[["length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2", "p_site", "rpfs"]].astype(int)

        adj_tis_offset.sort_values(['length'], inplace=True)
        adj_tis_offset["ribo"] = ["first"] * self.nt_num
        adj_tis_offset.to_csv(self.output_prefix + "_SSCBM_offset.txt", sep='\t', index=False)

    def draw_tis_heatmap(self):
        # draw the offset heatmap
        def draw_figure(tis_5end, tis_3end, tts_5end, tts_3end, out_pdf, out_png):
            matplotlib.use('Agg')
            now_cmap = 'Blues'

            fig = plt.figure(figsize=(12, 8), dpi=300)
            # start codon 5 end
            ax1 = plt.subplot(2, 2, 1)
            sns.heatmap(data=tis_5end, annot=None, linewidths=0.5, ax=ax1, cmap=now_cmap)
            # label_y = ax1.get_yticklabels()
            # plt.setp(label_y, rotation=0, horizontalalignment='right')
            ax1.set_title('RPFs 5end')
            ax1.set_ylabel("RPFs length")
            ax1.set_xlabel("from start codon")
            # start codon 3 end
            ax2 = plt.subplot(2, 2, 2)
            sns.heatmap(data=tis_3end, annot=None, linewidths=0.5, ax=ax2, cmap=now_cmap)
            # label_y = ax1.get_yticklabels()
            # plt.setp(label_y, rotation=0, horizontalalignment='right')
            ax2.set_title('RPFs 3end')
            ax2.set_ylabel("RPFs length")
            ax2.set_xlabel("from start codon")
            # stop codon 5 end
            ax3 = plt.subplot(2, 2, 3)
            sns.heatmap(data=tts_5end, annot=None, linewidths=0.5, ax=ax3, cmap=now_cmap)
            # label_y = ax1.get_yticklabels()
            # plt.setp(label_y, rotation=0, horizontalalignment='right')
            ax3.set_title('RPFs 5end')
            ax3.set_ylabel("RPFs length")
            ax3.set_xlabel("from stop codon")
            # stop codon 3 end
            ax4 = plt.subplot(2, 2, 4)
            sns.heatmap(data=tts_3end, annot=None, linewidths=0.5, ax=ax4, cmap=now_cmap)
            # label_y = ax1.get_yticklabels()
            # plt.setp(label_y, rotation=0, horizontalalignment='right')
            ax4.set_title('RPFs 3end')
            ax4.set_ylabel("RPFs length")
            ax4.set_xlabel("from stop codon")

            plt.tight_layout()
            # plt.show()
            fig.savefig(fname=out_pdf)
            fig.savefig(fname=out_png)

        # draw the raw offset heatmap
        out_pdf = self.output_prefix + "_SSCBM_offset.pdf"
        out_png = self.output_prefix + "_SSCBM_offset.png"
        draw_figure(self.tis_5end, self.tis_3end, self.tts_5end, self.tts_3end, out_pdf, out_png)

        # draw the scaled offset heatmap
        out_pdf_s = self.output_prefix + "_SSCBM_offset_scale.pdf"
        out_png_s = self.output_prefix + "_SSCBM_offset_scale.png"
        tis_5end_s = self.offset_scale(self.tis_5end)
        tis_3end_s = self.offset_scale(self.tis_3end)
        tts_5end_s = self.offset_scale(self.tts_5end)
        tts_3end_s = self.offset_scale(self.tts_3end)
        draw_figure(tis_5end_s, tis_3end_s, tts_5end_s, tts_3end_s, out_pdf_s, out_png_s)


    def make_frame_offset(self):
        '''
        make the offset dict

        1. peak of ribosome protect fragments about 30nt, 
        the read length distribution also varies significantly due to different protocols and species.
        So shift 1 nt (depend on the reads length) for the both 5/3 reads end to solve this problem (5end + 3end = 2).
        Empirical value: 2nt for Eukaryotes, 1nt for prokaryotes
        
        2. set the minimum and maximum window for p-site location
        3. set the first window for p-site calculation

        '''

        # assign the empirical value for the offset start
        self.exp_offset = 11 + self.peak_length - 30
        
        if self.adj_tis_offset:
            for length in range(self.min_length, self.max_length + 1):
                self.exp_offset = self.adj_tis_offset[length][1]

                left_offset1, left_offset2, left_offset3 = self.exp_offset, self.exp_offset + 1, self.exp_offset + 2

                self.frame_offset_len[length] = [left_offset1, left_offset2, left_offset3]
                self.frame_offset[length] = [0, 0, 0]
                
        else:
            for length in range(self.min_length, self.max_length + 1):
                shift_nt = (length - self.peak_length) // self.shift_nt
                if shift_nt < -3:
                    shift_nt = -3
                elif shift_nt > 6:
                    shift_nt = 6
                else:
                    pass

                left_offset1 = self.exp_offset + shift_nt
                left_offset2 = self.exp_offset + 1 + shift_nt
                left_offset3 = self.exp_offset + 2 + shift_nt

                self.frame_offset_len[length] = [left_offset1, left_offset2, left_offset3]
                self.frame_offset[length] = [0, 0, 0]

    def calc_frame(self, map_start, cds_start, reads_length, num1, num2, num3):
        '''
        calculate the offset frame of each reads

        1. trim the offset_length dependent on the reads length
        2. fit the offset length of each reads in the codon frame
        3. arrange the [0, 1, 2] frame reads into frame_offset dictionary
        '''

        if (map_start + self.frame_offset_len[reads_length][num1] - cds_start) % 3 == 0:
            try:
                self.frame_offset[reads_length][num1] += 1
            except KeyError:
                pass
        elif (map_start + self.frame_offset_len[reads_length][num2] - cds_start) % 3 == 0:
            try:
                self.frame_offset[reads_length][num2] += 1
            except KeyError:
                pass
        elif (map_start + self.frame_offset_len[reads_length][num3] - cds_start) % 3 == 0:
            try:
                self.frame_offset[reads_length][num3] += 1
            except KeyError:
                pass

    def get_mono_frame(self, mrna_attr, cds_start):
        '''
        through the entire bam file,
        1. reterieve the reads map start site of gene body
        2. Remove any reads that are too long or too short
        3. arrange predicted psite for each reads
        '''

        for line in mrna_attr.bam:
            map_start, map_end = line.get_blocks()[0]
            reads_length = line.infer_read_length()
            if self.min_length <= reads_length <= self.max_length:
                self.calc_frame(map_start, cds_start, reads_length, 0, 1, 2)

    def get_frame_offset(self):
        '''
        through the entire mrna dict file,
        1. define the frame offset dict file
        2. eliminate the transcript without read mapped
        3. set the gene body range for offset calculation, trim the stop codon from gene CDS
        4. check the category of sequencing profile
        monosome, disome, trisome
        5. disome and trisome need to perform two and three times psite calculations respectively

        '''

        self.make_frame_offset()

        for _, mrna_attr in self.mrna_dict.items():
            
            if len(mrna_attr.bam) == 0:
                continue
            else:
                cds_start, cds_end = mrna_attr.utr5_length, mrna_attr.utr5_length + mrna_attr.cds_length - 3 - 2
                self.get_mono_frame(mrna_attr, cds_start)


    def format_frame_offset(self):
        columns = ["length", "frame0", "rpfs0", "frame1", "rpfs1", "frame2", "rpfs2", "p_site", "rpfs"]

        # different flag for the disome and trisome sequencing profiles
        # next version
        ribo = [["first"], ["second"], ["third"]]
        frame_offset = pd.DataFrame(self.frame_offset).T
        rows, cols = frame_offset.shape[0], frame_offset.shape[1]
        
        def get_offset_len(temp_offset):
            for k, v in dict(temp_offset.idxmax(1)).items():
                if v == 0:
                    temp_offset.loc[k, "p_site"] = self.frame_offset_len[k][v] + 1
                elif v == 1:
                    temp_offset.loc[k, "p_site"] = self.frame_offset_len[k][v] + 1
                elif v == 2:
                    temp_offset.loc[k, "p_site"] = self.frame_offset_len[k][v] + 1 - 3

        # calculate the frame offset for each ribosome profile
        for part in range(0, cols, 3):
            temp_offset = frame_offset.iloc[:, part:part + 3].copy()
            get_offset_len(temp_offset)
            temp_offset["sum"] = temp_offset.apply(lambda x: x.iloc[0:3].sum(), axis=1)
            temp_offset.columns = ["rpfs0", "rpfs1", "rpfs2", "p_site", "rpfs"]
            temp_offset[["frame0", "frame1", "frame2"]] = [i[part:part + 3] for i in self.frame_offset_len.values()]
            temp_offset["length"] = temp_offset.index
            self.merge_frame_offset = pd.concat([self.merge_frame_offset, temp_offset[columns]])

        self.merge_frame_offset["periodicity"] = self.merge_frame_offset[["rpfs0", "rpfs1", "rpfs2"]].max(axis = 1) / self.merge_frame_offset["rpfs"]
        self.merge_frame_offset["periodicity"] = self.merge_frame_offset["periodicity"] * 100
        self.merge_frame_offset["periodicity"] = self.merge_frame_offset["periodicity"].astype(float).round(2)
        self.merge_frame_offset["ribo"] = ribo[part // 3] * rows
        self.merge_frame_offset[columns] = self.merge_frame_offset[columns].astype(int)


    def adjust_frame_offset(self):
        '''
        The length of reads varies in different sequencing files,
        so offset should be normalise to suitable length.

        1. first we need to check the peak length from offset dict

        2. In general, the offset codon frame of shorter reads will not be larger than the offset of peak reads length commonly,
        so the abnormal offset will shift to smaller offset in the range of one codon frame.

        3. Similarly, an offset codon frame of longer reads will not be shorter than peak reads length.

        '''

        # adjust the frame shift of reads shorter than peak length
        for length in range(self.peak_length, self.min_length, -1):
            psite1 = self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length, 'p_site'].reset_index(drop=True)[0]
            psite2 = self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length - 1, 'p_site'].reset_index(drop=True)[0]

            if psite2 + 2 < psite1:
                self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length - 1, 'p_site'] = psite2 + 3
            elif psite2 - 2 > psite1:
                self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length - 1, 'p_site'] = psite2 - 3
            else:
                continue

        # adjust the frame shift of reads longer than peak length
        for length in range(self.peak_length, self.max_length, 1):
            psite1 = self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length, 'p_site'].reset_index(drop=True)[0]
            psite2 = self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length + 1, 'p_site'].reset_index(drop=True)[0]
            
            if psite2 - 2 > psite1:
                self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length + 1, 'p_site'] = psite2 - 3
            elif psite2 + 2 < psite1:
                self.merge_frame_offset.loc[self.merge_frame_offset['length'] == length + 1, 'p_site'] = psite2 + 3
            else:
                continue

    def write_frame_offset(self):
        '''
        arrange the frame offset with predicted psite
        '''
        for length in range(0, self.merge_frame_offset.shape[0]):
            offset_line = self.merge_frame_offset.iloc[length, :]
            offset_psite = offset_line.iloc[[1, 3, 5]]
            offset_rpfs = offset_line.iloc[[2, 4, 6]]
        
            offset_index = (offset_psite - offset_line.iloc[-4] + 1) % 3

            offset_psite.index = offset_index
            offset_psite = offset_psite.loc[[0, 1, 2]]

            offset_rpfs.index = offset_index
            offset_rpfs = offset_rpfs.loc[[0, 1, 2]]

            self.merge_frame_offset.iloc[length, [1, 3, 5]] = offset_psite
            self.merge_frame_offset.iloc[length, [2, 4, 6]] = offset_rpfs

        self.merge_frame_offset.to_csv(self.output_prefix + "_RSBM_offset.txt", sep='\t', index=False)

    def draw_frame_heatmap(self):
        out_pdf = self.output_prefix + "_RSBM_offset.pdf"
        out_png = self.output_prefix + "_RSBM_offset.png"

        # reterieve the psite and rpfs
        raw_frame_offset = self.merge_frame_offset.loc[:, ["length", "p_site", 'rpfs0', 'rpfs1', 'rpfs2']]

        # rename the ylabel
        offset_num = raw_frame_offset.shape[0]
        raw_frame_offset.index = [str(raw_frame_offset["length"].to_list()[i]) + '_' + str(raw_frame_offset["p_site"].to_list()[i]) for i in range(0, offset_num)]

        # format the offset table and rename the columns
        raw_frame_offset = raw_frame_offset.drop(columns=["length", "p_site"])
        raw_frame_offset.columns = ["frame0", "frame1", "frame2"]
        raw_frame_offset = raw_frame_offset.apply(pd.to_numeric)

        # scale the frame offset to range [0, 1]
        scale_frame_offset = self.offset_scale(raw_frame_offset)
        scale_frame_offset = scale_frame_offset.apply(pd.to_numeric)

        # draw the frame offset heatmap
        matplotlib.use('Agg')
        now_cmap = 'Blues'

        fig = plt.figure(figsize=(12, 6), dpi=300)

        # raw RPFs frame offset
        ax1 = plt.subplot(1, 2, 1)
        sns.heatmap(data=raw_frame_offset, annot=None, linewidths=0.5, ax=ax1, cmap=now_cmap)
        ax1.set_title('raw counts')
        ax1.set_ylabel("RPFs length")
        ax1.set_xlabel("open reading frame")

        # scale RPFs frame offset
        ax2 = plt.subplot(1, 2, 2)
        sns.heatmap(data=scale_frame_offset, annot=None, linewidths=0.5, ax=ax2, cmap=now_cmap)
        ax2.set_title('scaled counts')
        ax2.set_ylabel("RPFs length")
        ax2.set_xlabel("open reading frame")

        plt.tight_layout()
        # plt.show()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
