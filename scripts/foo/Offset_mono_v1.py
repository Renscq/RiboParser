#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Offset_mono_v1.py


import os.path
import sys
from collections import OrderedDict
from itertools import islice

import pandas as pd
import pysam


class Mrna(object):

    def __init__(self, record):
        self.chromosome = record[0]
        self.gene_id = record[1]
        self.gene_name = record[2]
        self.transcript_id = record[3]
        self.start = record[4]
        self.end = record[5]
        self.utr5_length = int(record[6])
        self.cds_length = int(record[7])
        self.utr3_length = int(record[8])
        self.strand = record[9]
        self.rep_transcript = record[10]
        self.modified = record[11]
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
        self.column_name = ["length", "from_tis", "tis_5end", "tis_5end_per",
                            "from_tts", "tts_3end", "tts_3end_per", "sum"]
        # self.tis_offset = {"start_codon": OrderedDict(), "stop_codon": OrderedDict()}
        self.tis_offset = {"tis_5end": OrderedDict(), "tis_3end": OrderedDict(),
                           "tts_5end": OrderedDict(), "tts_3end": OrderedDict()}
        self.adj_tis_offset = pd.DataFrame(columns=self.column_name)
        self.frame_offset = OrderedDict()
        self.adj_frame_offset = pd.DataFrame(columns=["length", "frame1", "rpf_1",
                                                      "frame2", "rpf_2", "frame3", "rpf_3", "sum", "p_site", "ribo"])
        self.frame_offset_len = {}
        self.mode = args.mode

        self.mrna_dict = OrderedDict()
        self.length_dict = {}

        self.sample_file = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]

        self.pysam_input = None
        self.pysam_output = None
        self.output_prefix = None

        self.bam_file_list = None
        self.sample_dict = None

        self.peak_length = args.peak_length
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
            sys.stdout.writelines("import file: {bam}.\n".format(bam=self.sample_file))
            read_format = 'rb'
        elif self.sample_format.lower() == ".sam":
            sys.stdout.writelines("import file: {bam}.\n".format(bam=self.sample_file))
            read_format = 'r'
        else:
            sys.stdout.writelines("Unknown file format, please input the correct bam or sam file.")
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

        peak_length, self.peak_reads = max(self.length_dict.items(), key=lambda length: length[1])
        if self.peak_length:
            pass
        else:
            self.peak_length = peak_length

    @staticmethod
    def check_zero_offset(temp_offset):
        # offsets without reads coverage will examined and corrected.
        if temp_offset['from_tis'] == 0:
            temp_offset['from_tis'] = temp_offset['length'] - temp_offset['from_tts'] - 1
            return temp_offset

        elif temp_offset['from_tts'] == 0:
            temp_offset['from_tts'] = temp_offset['length'] - temp_offset['from_tis'] - 1
            return temp_offset

        else:
            return temp_offset

    def get_tis_offset(self):
        # build the offset table
        for number in range(self.min_length, self.max_length + 1):
            self.tis_offset["start_codon"][number] = {"tis_5end": {}, "tis_3end": {}}
            self.tis_offset["stop_codon"][number] = {"tts_5end": {}, "tts_3end": {}}
        # open the output file flag
        self.pysam_output = pysam.AlignmentFile(self.output_prefix + "_mrna.bam", 'wb', template=self.pysam_input)

        for mrna, mrna_attr in self.mrna_dict.items():
            # eliminate the transcript without read mapped
            if len(mrna_attr.bam) == 0:
                continue
            else:
                # get the TIS and TTS from mRNA dict.
                cds_start, cds_end = mrna_attr.utr5_length + 1, mrna_attr.utr5_length + mrna_attr.cds_length
                for line in mrna_attr.bam:
                    self.pysam_output.write(line)
                    map_start, map_end = line.get_blocks()[0]
                    # map start is not included, so add 1 nt.
                    # map_start = map_start + 1
                    read_length = line.infer_read_length()

                    # reads that match a specific length are retained.
                    if self.min_length <= read_length <= self.max_length:
                        # specifies the possible range of offset lengths around the start codon.
                        # start codon position like this [8, 9, 10, 11, "12", 13, 14, 15, 16]
                        if map_start + 9 <= cds_start <= map_end - 14:
                            offset_5end = cds_start - map_start
                            offset_3end = map_end - cds_start
                            try:
                                self.tis_offset["start_codon"][read_length]["tis_5end"][offset_5end] += 1
                            except KeyError:
                                self.tis_offset["start_codon"][read_length]["tis_5end"][offset_5end] = 1
                            try:
                                self.tis_offset["start_codon"][read_length]["tis_3end"][offset_3end] += 1
                            except KeyError:
                                self.tis_offset["start_codon"][read_length]["tis_3end"][offset_3end] = 1

                        # specifies the possible range of offset lengths around the stop codon.
                        elif map_start + 14 <= cds_end <= map_end - 9:
                            # shift -5 nt to align the last AA.  like the list [-6, -5, -4, -3, -2, -1, 0, 1, 2]
                            offset_5end = (cds_end - 2 - 3) - map_start
                            offset_3end = map_end - (cds_end - 2 - 3)
                            try:
                                self.tis_offset["stop_codon"][read_length]["tts_5end"][offset_5end] += 1
                            except KeyError:
                                self.tis_offset["stop_codon"][read_length]["tts_5end"][offset_5end] = 1
                            try:
                                self.tis_offset["stop_codon"][read_length]["tts_3end"][offset_3end] += 1
                            except KeyError:
                                self.tis_offset["stop_codon"][read_length]["tts_3end"][offset_3end] = 1
                    else:
                        pass

        self.pysam_output.close()
        self.pysam_input.close()

    def get_tis_offset2(self):
        # build the offset table
        for number in range(self.min_length, self.max_length + 1):
            self.tis_offset["start_codon"][number] = {"tis_5end": {}, "tis_3end": {}}
            self.tis_offset["stop_codon"][number] = {"tts_5end": {}, "tts_3end": {}}

        for mrna, mrna_attr in self.mrna_dict.items():
            # eliminate the transcript without read mapped
            if len(mrna_attr.bam) == 0:
                continue
            else:
                # get the TIS and TTS from mRNA dict.
                cds_start, cds_end = mrna_attr.utr5_length + 1, mrna_attr.utr5_length + mrna_attr.cds_length
                for line in mrna_attr.bam:
                    map_start, map_end = line.get_blocks()[0]
                    # map start is not included, so add 1 nt.
                    # map_start = map_start + 1
                    read_length = line.infer_read_length()

                    # reads that match a specific length are retained.
                    if self.min_length <= read_length <= self.max_length:
                        # specifies the possible range of offset lengths around the start codon.
                        # start codon position like this [8, 9, 10, 11, "12", 13, 14, 15, 16]

                        if map_start <= cds_start <= map_end:
                            offset_5end = cds_start - map_start
                            offset_3end = map_end - cds_start
                            try:
                                self.tis_offset["start_codon"][read_length]["tis_5end"][offset_5end] += 1
                            except KeyError:
                                self.tis_offset["start_codon"][read_length]["tis_5end"][offset_5end] = 1
                            try:
                                self.tis_offset["start_codon"][read_length]["tis_3end"][offset_3end] += 1
                            except KeyError:
                                self.tis_offset["start_codon"][read_length]["tis_3end"][offset_3end] = 1

                        # specifies the possible range of offset lengths around the stop codon.
                        elif map_start + 14 <= cds_end <= map_end - 9:
                            # shift -5 nt to align the last AA.  like the list [-6, -5, -4, -3, -2, -1, 0, 1, 2]
                            offset_5end = (cds_end - 2 - 3) - map_start
                            offset_3end = map_end - (cds_end - 2 - 3)
                            try:
                                self.tis_offset["stop_codon"][read_length]["tts_5end"][offset_5end] += 1
                            except KeyError:
                                self.tis_offset["stop_codon"][read_length]["tts_5end"][offset_5end] = 1
                            try:
                                self.tis_offset["stop_codon"][read_length]["tts_3end"][offset_3end] += 1
                            except KeyError:
                                self.tis_offset["stop_codon"][read_length]["tts_3end"][offset_3end] = 1
                    else:
                        pass

        self.pysam_input.close()

    def align_to_start_codon(self, length, start_codon, stop_codon):

        start_reads = start_codon.loc[start_codon["length"] == length].copy()
        stop_reads = stop_codon.loc[stop_codon["length"] == length].copy()

        # calculate the complementary length of offset and align it to the start codon.
        del stop_reads["length"]
        stop_reads["from_tis"] = length - 1 - stop_reads["from_tts"]

        shift_nt = round((length - self.peak_length) / 2)
        # shift_nt = shift_nt if shift_nt <= 3 else 3
        candidate_offset = [i for i in range(9 + shift_nt, 15 + shift_nt)]
        start_reads = start_reads[start_reads["from_tis"].isin(candidate_offset)]
        stop_reads = stop_reads[stop_reads["from_tis"].isin(candidate_offset)]

        # get the percentage of each candidate offset.
        start_reads["tis_5end_per"] = start_reads["tis_5end"].div(start_reads["tis_5end"].sum()) * 100
        stop_reads["tts_3end_per"] = stop_reads["tts_3end"].div(stop_reads["tts_3end"].sum()) * 100

        # alter the merged dataframe.
        temp_offset = pd.merge(start_reads, stop_reads, how='outer', left_on="from_tis", right_on="from_tis")
        temp_offset["length"] = length
        temp_offset = temp_offset.fillna(0)
        temp_offset = temp_offset.apply(self.check_zero_offset, axis=1)
        temp_offset["sum"] = temp_offset["tis_5end"] + temp_offset["tts_3end"]
        temp_offset = temp_offset[self.column_name]

        return temp_offset

    def align_to_tis(self, length, start_codon, stop_codon):

        start_reads = start_codon.loc[start_codon["length"] == length].copy()
        stop_reads = stop_codon.loc[stop_codon["length"] == length].copy()

        # calculate the complementary length of offset and align it to the start codon.
        del stop_reads["length"]
        stop_reads["from_tis"] = length - 1 - stop_reads["from_tts"]

        shift_nt = round((length - self.peak_length) / 2)
        # shift_nt = shift_nt if shift_nt <= 3 else 3
        candidate_offset = [i for i in range(9 + shift_nt, 15 + shift_nt)]
        start_reads = start_reads[start_reads["from_tis"].isin(candidate_offset)]
        stop_reads = stop_reads[stop_reads["from_tis"].isin(candidate_offset)]

        # get the percentage of each candidate offset.
        start_reads["tis_5end_per"] = start_reads["tis_5end"].div(start_reads["tis_5end"].sum()) * 100
        stop_reads["tts_3end_per"] = stop_reads["tts_3end"].div(stop_reads["tts_3end"].sum()) * 100

        # alter the merged dataframe.
        temp_offset = pd.merge(start_reads, stop_reads, how='outer', left_on="from_tis", right_on="from_tis")
        temp_offset["length"] = length
        temp_offset = temp_offset.fillna(0)
        temp_offset = temp_offset.apply(self.check_zero_offset, axis=1)
        temp_offset["sum"] = temp_offset["tis_5end"] + temp_offset["tts_3end"]
        temp_offset = temp_offset[self.column_name]

        return temp_offset

    def align_to_stop_codon(self, length, start_codon, stop_codon):

        start_reads = start_codon.loc[start_codon["length"] == length].copy()
        stop_reads = stop_codon.loc[stop_codon["length"] == length].copy()

        # calculate the complementary length of offset and align it to the start codon.
        del start_reads["length"]
        start_reads["from_tts"] = length - 1 - start_reads["from_tis"]
        shift_nt = round((length - self.peak_length) / 2)
        # shift_nt = shift_nt if shift_nt <= 3 else 3
        candidate_offset = [i for i in range(14 + shift_nt, 20 + shift_nt)]
        start_reads = start_reads[start_reads["from_tts"].isin(candidate_offset)]
        stop_reads = stop_reads[stop_reads["from_tts"].isin(candidate_offset)]

        # get the percentage of each candidate offset.
        start_reads["tis_5end_per"] = start_reads["tis_5end"].div(start_reads["tis_5end"].sum()) * 100
        stop_reads["tts_3end_per"] = stop_reads["tts_3end"].div(stop_reads["tts_3end"].sum()) * 100

        # alter the merged dataframe.
        temp_offset = pd.merge(start_reads, stop_reads, how='outer', left_on="from_tts", right_on="from_tts")
        temp_offset["length"] = length
        temp_offset = temp_offset.fillna(0)
        temp_offset = temp_offset.apply(self.check_zero_offset, axis=1)
        temp_offset["sum"] = temp_offset["tis_5end"] + temp_offset["tts_3end"]
        temp_offset = temp_offset[self.column_name]

        return temp_offset

    def format_tis_offset(self):
        self.adj_tis_offset[["tis_5end_per", "tts_3end_per"]] = self.adj_tis_offset[
            ["tis_5end_per", "tts_3end_per"]].astype(float).round(2)
        self.adj_tis_offset[["from_tis", "tis_5end", "from_tts", "tts_3end", "sum"]] = self.adj_tis_offset[
            ["from_tis", "tis_5end", "from_tts", "tts_3end", "sum"]].astype(int)

    def adjust_tis_offset(self, start_codon, stop_codon):
        column_name = ["length", "from_tis", "tis_5end", "tis_5end_per", "from_tts", "tts_3end", "tts_3end_per", "sum"]
        # monosome data only needs to align p-site on one side.
        for length in range(self.min_length, self.max_length + 1):
            start_offset = self.align_to_start_codon(length, start_codon, stop_codon)
            start_offset_peak = pd.DataFrame(start_offset.sort_values("sum", ascending=False).iloc[0]).T
            self.adj_tis_offset = pd.concat([self.adj_tis_offset, start_offset_peak])
        self.format_tis_offset()
        self.adj_tis_offset["p_site"] = self.adj_tis_offset["from_tis"] + 1
        self.adj_tis_offset["ribo"] = ["first"] * self.nt_num

    def write_tis_offset(self):

        start_codon = pd.DataFrame(columns=("length", "from_tis", "tis_5end"))
        for length, data in self.tis_offset['start_codon'].items():
            df = pd.DataFrame(data)
            df = df.fillna(0)
            df["length"] = length
            df["from_tis"] = df.index
            start_codon = pd.concat([start_codon, df[["length", "from_tis", "tis_5end"]]])

        stop_codon = pd.DataFrame(columns=("length", "from_tts", "tts_3end"))
        for length, data in self.tis_offset['stop_codon'].items():
            df = pd.DataFrame(data)
            df = df.fillna(0).astype(int)
            df["length"] = length
            df["from_tts"] = df.index
            stop_codon = pd.concat([stop_codon, df[["length", "from_tts", "tts_3end"]]])

        self.adjust_tis_offset(start_codon, stop_codon)
        self.adj_tis_offset.sort_values(['ribo', 'length'], inplace=True)

        self.adj_tis_offset = self.adj_tis_offset.reset_index(drop=True)
        psite_num = len(self.adj_tis_offset.index)
        for rows in range(psite_num):
            if rows > psite_num - 2:
                break
            elif self.adj_tis_offset.loc[rows, 'p_site'] - self.adj_tis_offset.loc[rows + 1, 'p_site'] > 1:
                self.adj_tis_offset.loc[rows, 'p_site'] -= 3
        self.adj_tis_offset.to_csv(self.output_prefix + "_tis_offset.txt", sep='\t', index=False)

    def calc_frame(self, map_start, cds_start, reads_length, num1, num2, num3):
        if (map_start + self.frame_offset_len[reads_length][num1] - cds_start + 1) % 3 == 0:
            try:
                self.frame_offset[reads_length][num1] += 1
            except KeyError:
                pass
        elif (map_start + self.frame_offset_len[reads_length][num2] - cds_start + 1) % 3 == 0:
            try:
                self.frame_offset[reads_length][num2] += 1
            except KeyError:
                pass
        elif (map_start + self.frame_offset_len[reads_length][num3] - cds_start + 1) % 3 == 0:
            try:
                self.frame_offset[reads_length][num3] += 1
            except KeyError:
                pass

    def get_mono_frame(self, mrna_attr, cds_start):
        for line in mrna_attr.bam:
            self.pysam_output.write(line)
            map_start, map_end = line.get_blocks()[0]
            # map start is not included, so add 1 nt.
            # map_start = map_start + 1
            reads_length = line.infer_read_length()
            if self.min_length <= reads_length <= self.max_length:
                self.calc_frame(map_start, cds_start, reads_length, 0, 1, 2)

    def make_frame_offset(self):
        for number in range(self.min_length, self.max_length + 1):
            # shift 1 nt for the both 5/3 read end (5end + 3end = 2)
            length_normalised = (number - self.peak_length) // 2
            left_offset1, left_offset2, left_offset3 = 11 + length_normalised, 12 + length_normalised, 13 + length_normalised

            self.frame_offset_len[number] = [left_offset1, left_offset2, left_offset3]
            self.frame_offset[number] = [0, 0, 0]

    def get_frame_offset(self):
        self.make_frame_offset()
        self.pysam_output = pysam.AlignmentFile(self.output_prefix + "_mrna.bam", 'wb', template=self.pysam_input)
        for mrna, mrna_attr in self.mrna_dict.items():
            # eliminate the transcript without read mapped
            if len(mrna_attr.bam) == 0:
                continue
            else:
                cds_start, cds_end = mrna_attr.utr5_length + 1, mrna_attr.utr5_length + mrna_attr.cds_length - 3
                self.get_mono_frame(mrna_attr, cds_start)

        self.pysam_output.close()
        self.pysam_input.close()

    def write_frame_offset(self):
        columns = ["length", "frame1", "rpf_1", "frame2", "rpf_2", "frame3", "rpf_3", "sum", "p_site", "ribo"]

        ribo = [["first"], ["second"], ["third"]]
        frame_offset = pd.DataFrame(self.frame_offset).T
        rows, cols = frame_offset.shape[0], frame_offset.shape[1]
        for part in range(0, cols, 3):
            temp_offset = frame_offset.iloc[:, part:part + 3].copy()
            temp_offset["p_site"] = [self.frame_offset_len[k][v] + 1 for k, v in dict(temp_offset.idxmax(1)).items()]
            temp_offset["sum"] = temp_offset.apply(lambda x: x.iloc[0:3].sum(), axis=1)
            temp_offset.columns = ["rpf_1", "rpf_2", "rpf_3", "p_site", "sum"]
            temp_offset[["frame1", "frame2", "frame3"]] = [i[part:part + 3] for i in self.frame_offset_len.values()]
            temp_offset["length"] = temp_offset.index
            temp_offset["ribo"] = ribo[part // 3] * rows
            self.adj_frame_offset = pd.concat([self.adj_frame_offset, temp_offset[columns]])

        self.adj_frame_offset.to_csv(self.output_prefix + "_frame_offset.txt", sep='\t', index=False)
