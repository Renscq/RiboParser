#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Ensembl_Ref.py


from collections import OrderedDict


class Transcripts(object):

    def __init__(self, section):
        self.mess = [section["mess"]]

        self.chrom = section["chrom"]
        self.source = section["source"]
        self.feature = section["feature"]

        self.start = section["start"]
        self.end = section["end"]
        self.strand = section["strand"]
        self.attr = section["attr"]
        self.attr_dict = section["attr_dict"]

        self.gene_biotype = section["attr_dict"]["gene_biotype"]
        self.transcript_biotype = section["attr_dict"]["transcript_biotype"]
        self.gene_id = section["attr_dict"]["gene_id"]
        self.transcript_id = section["attr_dict"]["transcript_id"]
        self.transcript_name = None
        self.gene_name = None

        self.exon_feature = []
        self.cds_feature = []
        self.five_prime_utr = []
        self.three_prime_utr = []
        self.start_codon = []
        self.stop_codon = []
        self.exons = []
        self.cds = []
        self.utr5 = 0
        self.utr3 = 0
        self.cds_length = 0

        self.rep_transcript = False

        self.modified = False

        def get_ids():
            try:
                self.gene_name = section["attr_dict"]["gene_name"]
            except KeyError:
                self.gene_name = "-"
            try:
                self.transcript_name = section["attr_dict"]["transcript_name"]
            except KeyError:
                self.transcript_name = "-"

        get_ids()

    def add_feature(self, section, cds_shift):
        if section["feature"] == "exon":
            self.mess.append(section["mess"])
            self.exon_feature.append(section["mess"].split('\t'))
            self.exons.append([section["start"], section["end"]])

            if len(self.exons) >= 2 and section["start"] < self.exons[-2][1] and self.strand == "+":
                self.exon_feature.sort(key=lambda exon: exon[3], reverse=False)
                self.exons.sort(key=lambda exon: exon[0], reverse=False)

            elif len(self.exons) >= 2 and section["end"] > self.exons[-2][1] and self.strand == "-":
                self.exon_feature.sort(key=lambda exon: exon[3], reverse=True)
                self.exons.sort(key=lambda exon: exon[0], reverse=True)

        elif section["feature"] == "CDS":
            self.mess.append(section["mess"])
            self.cds_feature.append(section["mess"].split('\t'))
            self.cds.append([section["start"], section["end"]])
            self.cds_length += (section["end"] - section["start"] + cds_shift)

            if len(self.cds) >= 2 and section["start"] < self.cds[-2][1] and self.strand == "+":
                self.cds_feature.sort(key=lambda cds: cds[3], reverse=False)
                self.cds.sort(key=lambda cds: cds[0], reverse=False)

            elif len(self.cds) >= 2 and section["end"] > self.cds[-2][1] and self.strand == "-":
                self.cds_feature.sort(key=lambda cds: cds[3], reverse=True)
                self.cds.sort(key=lambda cds: cds[0], reverse=True)

        elif section["feature"] == "start_codon":
            self.start_codon = [section["start"], section["end"]]

        elif section["feature"] == "stop_codon":
            self.stop_codon = [section["start"], section["end"]]

    def check_codon(self):

        if self.strand == '+':
            if self.start_codon and self.start_codon[0] < self.cds[0][0]:
                if self.cds[0][0] - self.start_codon[0] > 3:
                    pass
                else:
                    self.cds_feature[0][3] = str(self.start_codon[0])
                self.cds[0][0] -= 3
                self.cds_length += 3
            if self.stop_codon and self.stop_codon[1] > self.cds[-1][1]:
                if self.stop_codon[1] - self.cds[-1][1] > 3:
                    pass
                else:
                    self.cds_feature[-1][4] = str(self.stop_codon[1])
                    self.cds[-1][1] += 3
                    self.cds_length += 3

        elif self.strand == '-':
            if self.start_codon and self.start_codon[1] > self.cds[0][1]:
                if self.start_codon[1] - self.cds[0][1] > 3:
                    pass
                else:
                    self.cds_feature[0][4] = str(self.start_codon[1])
                    self.cds[0][1] += 3
                    self.cds_length += 3
            if self.stop_codon and self.stop_codon[0] < self.cds[-1][0]:
                if self.cds[-1][0] - self.stop_codon[0] > 3:
                    pass
                else:
                    self.cds_feature[-1][3] = str(self.stop_codon[0])
                    self.cds[-1][0] -= 3
                    self.cds_length += 3

    def plus_utr5(self):
        self.utr5 = 0
        for exon in self.exons:
            if self.cds[0][0] > exon[1]:
                self.utr5 += exon[1] - exon[0] + 1
            elif exon[0] <= self.cds[0][0] <= exon[1]:
                self.utr5 += self.cds[0][0] - exon[0]

    def plus_utr3(self):
        self.utr3 = 0
        for exon in reversed(self.exons):
            if self.cds[-1][-1] < exon[0]:
                self.utr3 += exon[1] - exon[0] + 1
            elif exon[0] <= self.cds[-1][-1] <= exon[1]:
                self.utr3 += exon[1] - self.cds[-1][-1]

    def minus_utr5(self):
        self.utr5 = 0
        for exon in self.exons:
            if self.cds[0][-1] < exon[0]:
                self.utr5 += exon[1] - exon[0] + 1
            elif exon[0] <= self.cds[0][-1] <= exon[1]:
                self.utr5 += exon[1] - self.cds[0][-1]

    def minus_utr3(self):
        self.utr3 = 0
        for exon in reversed(self.exons):
            if self.cds[-1][0] > exon[1]:
                self.utr3 += exon[1] - exon[0] + 1
            elif exon[0] <= self.cds[-1][0] <= exon[1]:
                self.utr3 += self.cds[-1][0] - exon[0]

    def add_utr(self, utr_len, chroms_dict):
        self.check_codon()
        # check the plus strand
        if self.strand == '+':
            self.plus_utr5()
            self.plus_utr3()
            # weather the utr5 more than 30
            if self.utr5 < 90:
                if self.exons[0][0] - utr_len < 1:
                    self.exons[0][0] = 1
                    self.utr5 = self.cds[0][0] - 1
                    self.exon_feature[0][3] = 1
                    self.start = self.exons[0][0]
                else:
                    self.exons[0][0] = self.exons[0][0] - utr_len
                    self.utr5 = self.cds[0][0] - self.exons[0][0]
                    self.exon_feature[0][3] = self.exons[0][0]
                    self.start = self.exons[0][0]
                self.modified = utr_len
                self.plus_utr5()
            else:
                pass
            # weather the utr3 more than 30
            if self.utr3 < 90:
                if self.exons[-1][-1] + utr_len > chroms_dict[self.chrom].end:
                    self.exons[-1][-1] = chroms_dict[self.chrom].end
                    self.exon_feature[-1][4] = self.exons[-1][-1]
                    self.end = self.exons[-1][-1]
                else:
                    self.exons[-1][-1] = self.exons[-1][-1] + utr_len
                    self.exon_feature[-1][4] = self.exons[-1][-1]
                    self.end = self.exons[-1][-1]
                self.modified = utr_len
                self.plus_utr3()
            else:
                pass
        # 判断负链
        elif self.strand == '-':
            self.minus_utr5()
            self.minus_utr3()
            # weather the utr5 more than 30
            if self.utr5 < 90:
                if self.exons[0][-1] + utr_len > chroms_dict[self.chrom].end:
                    self.exons[0][-1] = chroms_dict[self.chrom].end
                    self.exon_feature[0][4] = self.exons[0][-1]
                    self.end = self.exons[0][-1]
                else:
                    self.exons[0][-1] = self.exons[0][-1] + utr_len
                    self.exon_feature[0][4] = self.exons[0][-1]
                    self.end = self.exons[0][-1]
                self.modified = utr_len
                self.minus_utr5()
            else:
                pass
            # weather the utr3 more than 30
            if self.utr3 < 90:
                if self.exons[-1][0] - utr_len < 1:
                    self.exons[-1][0] = 1
                    self.exon_feature[-1][3] = 1
                    self.start = 1
                else:
                    self.exons[-1][0] = self.exons[-1][0] - utr_len
                    self.exon_feature[-1][3] = self.exons[-1][0]
                    self.start = self.exons[-1][0]
                self.modified = utr_len
                self.minus_utr3()
            else:
                pass


class Gene(object):

    def __init__(self, section):
        self.section = section
        self.chrom = section["chrom"]
        self.source = section["source"]
        self.feature = section["feature"]
        self.gene_id = section["attr_dict"]["gene_id"]
        self.start = section["start"]
        self.end = section["end"]
        self.strand = section["strand"]
        self.attr = section["attr"]
        self.attr_dict = section["attr_dict"]

        self.gene_type = self.attr_dict["gene_biotype"]
        self.transcript = OrderedDict()
        self.rep_transcript = None
        self.rep_length = 0
        self.trans_num = 0
        self.modified = False

    def add_transcript(self, now_transcript):
        if now_transcript.cds_length > self.rep_length:
            self.rep_length = now_transcript.cds_length
            self.rep_transcript = now_transcript.transcript_id

        if now_transcript.modified:
            self.modified = True
            self.start = now_transcript.start
            self.end = now_transcript.end

        self.transcript[now_transcript.transcript_id] = now_transcript
        self.trans_num += 1


class Chrom(object):

    def __init__(self, line):
        self.attr = line.description
        self.chromosome = line.id
        self.start = 1
        self.end = len(line.seq)
        self.seq = line.seq
