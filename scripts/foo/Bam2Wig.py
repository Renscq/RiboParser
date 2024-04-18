#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Bam2Wig.py


import pysam
import os
import sys
import numpy as np
import pandas as pd
import polars as pl
from collections import OrderedDict


class Bam2Wig(object):

    def __init__(self, args):
        """
        arguments for these:
        bam file
        psite file
        multiple mapping
        bed file for density
        rpm file

        """
        self.bam = args.bam
        self.sample_suffix = os.path.splitext(args.bam)[1]
        self.read_format = None
        self.pysam_input = None

        self.psite = args.psite
        self.offset = {}

        self.tag_num = args.times
        self.secondary = args.secondary
        self.supplementary = args.supplementary
        self.format = args.format

        self.chrom = OrderedDict()
        self.bed_minus = {}
        self.bed_plus = {}
        self.bed_minus_df = None
        self.bed_plus_df = None

        self.norm = args.norm
        self.merge = args.merge
        self.output = args.output

    def read_offset(self):
        """
        import the offset file to locate the p-site of reads from bam file.

        1. import the offset file
        length psite
        27  12
        28  12
        29  13
        30  13
        ...

        2. make the offset dict for reads trim

        """
        offset = pd.read_csv(self.psite, sep='\t', header=0, names=None)
        offset.dropna(axis=0, inplace=True)
        offset["p_site"] = offset["p_site"].astype(int)
        offset["p_site"] = offset["p_site"] - 1
        self.offset = offset[["length", "p_site"]].groupby(offset["length"])["p_site"].apply(list).T.to_dict()
        print(self.offset, flush=True)

        # check the data type.
        if "third" in offset["ribo"].unique():
            self.profile = "trisome"
        elif "second" in offset["ribo"].unique():
            self.profile = "disome"
        elif "first" in offset["ribo"].unique():
            self.profile = "monosome"
        else:
            print("Unknown type of sequence profile, please re-check the offset file!", flush=True)
            sys.exit()

    def create_chroms_dict(self):
        """
        create a dict to save the chromosome and reads number
        """
        # self.chrom = dict(zip(self.pysam_input.header.references, self.pysam_input.header.lengths))
        for chrom in self.pysam_input.header.references:
            self.chrom[chrom] = 0
            self.bed_minus[chrom] = {}
            self.bed_plus[chrom] = {}

    def set_tag_num(self):
        # set the tag_num to discard the reads aligned too many loci
        if self.tag_num == 0:
            self.tag_num = len(self.mrna_dict)
        else:
            pass

    def convert_rpf_density(self):
        """
        import the bam file and convert the reads to density

        1. import the file
        2. remove the reads out of specified range
        3. remove the reads aligned to multiple loci
        3. retrieve the psite
        4. merge the psite density to chrom_dict

        """
        
        now_rows = 0
        for read in self.pysam_input.fetch(until_eof=True):

            # processing
            now_rows += 1
            if now_rows % 100000 == 0:
                print("now rows:" + str(now_rows), flush=True)
            else:
                pass

            # discard the unmapped reads
            if read.is_unmapped:
                continue

            # discard the secondary mapping
            if  read.is_secondary == self.secondary == True:
                continue
            
            # discard the supplementary mapping
            if  read.is_supplementary == self.secondary == True:
                continue

            # filter the unique reads
            if read.has_tag("NH") and read.get_tag("NH") > self.tag_num:
                continue
            
            # split positive and negative strand
            chrom = read.reference_name
            length = read.infer_read_length()
            
            # split positive and negative strand
            if read.is_reverse:
                try:
                    # start = read.aligned_pairs[self.offset[length][0]][1]
                    start = read.positions[-self.offset[length][0]]
                except KeyError:
                    continue

                self.chrom[chrom] += 1

                try:
                    self.bed_minus[chrom][start] += 1
                except KeyError:
                    self.bed_minus[chrom][start] = 1

            else:
                try:
                    # start = read.aligned_pairs[self.offset[length][0]][1]
                    start = read.positions[self.offset[length][0]]
                except KeyError:
                    continue

                self.chrom[chrom] += 1

                try:
                    self.bed_plus[chrom][start] += 1
                except KeyError:
                    self.bed_plus[chrom][start] = 1

        print("now rows:" + str(now_rows), flush=True)

    def import_bam(self):
        """
        check the file format
        import the alignment file with bam or sam file

        """
        if self.sample_suffix.lower() == ".bam":
            print("import file: {bam}.\n".format(bam=self.bam), flush=True)
            self.read_format = 'rb'
        elif self.sample_suffix.lower() == ".sam":
            print("import file: {bam}.\n".format(bam=self.bam), flush=True)
            self.read_format = 'r'
        else:
            print("Unknown file format, please input the correct bam or sam file.", flush=True)
            sys.exit()

        self.pysam_input = pysam.AlignmentFile(self.bam, self.read_format)
        
        self.create_chroms_dict()

        self.convert_rpf_density()
        
        self.pysam_input.close()
        

    def convert_dict_to_dataframe(self):
        """
        1. split minus and plus strand
        2. for each chrom, convert the dict to dataframe
        3. merge all chrom dataframe

        """

        minus_list = []
        for chrom, density in self.bed_minus.items():
            tmp_df = pd.DataFrame([density]).T.reset_index()
            tmp_df.columns = ["start", "density"]
            tmp_df.insert(loc=0, column="chr", value=chrom)

            minus_list.append(tmp_df)
        self.bed_minus_df = pd.concat(minus_list, axis=0)
        self.bed_minus_df.loc[:, "density"] = -self.bed_minus_df["density"]

        plus_list = []
        for chrom, density in self.bed_plus.items():
            tmp_df = pd.DataFrame([density]).T.reset_index()
            tmp_df.columns = ["start", "density"]
            tmp_df.insert(0, "chr", chrom)

            plus_list.append(tmp_df)
        self.bed_plus_df = pd.concat(plus_list, axis=0)

    def norm_rpm(self):
        """
        1. summary the bam reads count
        2. convert the density to rpm
        rpm = ( density * 1000000 ) / total reads
        """
        if self.norm:
            total_reads = pd.DataFrame([self.chrom]).sum()

            self.bed_minus_df.loc[:, "density"] = self.bed_minus_df["density"] * 1e6 / total_reads[0]
            self.bed_plus_df.loc[:, "density"] = self.bed_plus_df["density"] * 1e6 / total_reads[0]

            self.bed_minus_df.loc[:, "density"] = self.bed_minus_df["density"].round(3)
            self.bed_plus_df.loc[:, "density"] = self.bed_plus_df["density"].round(3)
        else:
            print('Skip!', flush=True)

    def output_bed(self):
        """
        1. drop all INF and NA values
        2. convert the ribosome density to bedgraph / wig format
        3. merge plus / minus strand density
        4. output the density file
        """
        self.bed_minus_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.bed_plus_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        
        self.bed_minus_df.dropna(inplace=True)
        self.bed_plus_df.dropna(inplace=True)

        if self.merge:
            bed_df = pd.concat([self.bed_plus_df, self.bed_minus_df], axis=0)
            bed_df.reset_index(drop=True).sort_values(["chr", "start"])

            if self.format == "bedgraph":
                bed_df.insert(loc=2, column="end", value=bed_df["start"] + 1)
                bed_df[["start", "end"]] = bed_df[["start", "end"]].astype(int)

                bed_pl = pl.DataFrame(bed_df)
                bed_pl.write_csv(self.output + ".bedgraph", separator = '\t', has_header = False)
            else:
                bed_df["start"] = bed_df["start"].astype(int)
                bed_pl = pl.DataFrame(bed_df)
                bed_pl.write_csv(self.output + ".wig", separator = '\t', has_header = False)

        else:
            if self.format == 'bedgraph':
                self.bed_minus_df.insert(loc=2, column="end", value=self.bed_minus_df["start"] + 1)
                self.bed_plus_df.insert(loc=2, column="end", value=self.bed_plus_df["start"] + 1)

                self.bed_minus_df[["start", "end"]] = self.bed_minus_df[["start", "end"]].astype(int)
                self.bed_plus_df[["start", "end"]] = self.bed_plus_df[["start", "end"]].astype(int)

                bed_minus_pl = pl.DataFrame(self.bed_minus_df)
                bed_plus_pl = pl.DataFrame(self.bed_plus_df)

                bed_minus_pl.write_csv(self.output + "_minus.bedgraph", separator = '\t', has_header = False)
                bed_plus_pl.write_csv(self.output + "_plus.bedgraph", separator = '\t', has_header = False)
            else:
                self.bed_minus_df["start"] = self.bed_minus_df["start"].astype(int)
                self.bed_plus_df["start"] = self.bed_plus_df["start"].astype(int)
            
                bed_minus_pl = pl.DataFrame(self.bed_minus_df)
                bed_plus_pl = pl.DataFrame(self.bed_plus_df)

                bed_minus_pl.write_csv(self.output + "_minus.wig", separator = '\t', has_header = False)
                bed_plus_pl.write_csv(self.output + "_plus.wig", separator = '\t', has_header = False)
