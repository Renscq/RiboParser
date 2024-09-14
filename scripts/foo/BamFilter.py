#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : BamFilter.py


import pysam
import os
import sys


class BamFilter():

    def __init__(self, args):
        
        # filename
        self.ibam = args.ibam
        self.obam = args.obam

        # option for reads filter
        self.length = args.length
        self.quality = args.quality
        self.secondary = args.unique
        self.supplementary = args.supplementary

        # file store
        self.input_bam = None
        self.flt_bam = None


    def get_length_list(self):

        self.length = [int(i) for i in self.length.split(',')]


    def filter_reads(self):
        """
        filte the mapped reads with three condition
        
        secondary mapping
        supplymentary mapping
        high quality
        fitted length

        """
        for read in self.input_bam.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            if read.is_secondary != self.secondary:
                continue

            if read.is_supplementary != self.supplementary:
                continue

            if read.mapping_quality < self.quality:
                continue

            if read.infer_read_length() in self.length:
                self.flt_bam.write(read)

    def import_bam(self):
        """
        check the file format
        import the alignment file with bam or sam file

        """
        sample_suffix = os.path.splitext(self.ibam)[1]

        if sample_suffix.lower() == ".bam":
            print("import file: {bam}.\n".format(bam=self.ibam), flush=True)
            read_format = 'rb'
        elif sample_suffix.lower() == ".sam":
            print("import file: {bam}.\n".format(bam=self.ibam), flush=True)
            read_format = 'r'
        else:
            print("Unknown file format, please input the correct bam or sam file.", flush=True)
            sys.exit()

        self.input_bam = pysam.AlignmentFile(self.ibam, read_format)

        self.flt_bam = pysam.AlignmentFile(self.obam, "wb", template=self.input_bam)

        self.get_length_list()
        self.filter_reads()

        self.input_bam.close()
        self.flt_bam.close()

