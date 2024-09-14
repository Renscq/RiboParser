#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Ribo.py


import os
import sys
from collections import OrderedDict
from itertools import islice
from Bio import SeqIO
from math import ceil
import numpy as np
import pandas as pd
import polars as pl
import pysam


class Mrna(object):

    def __init__(self, record):
        self.chromosome = record[0]
        self.gene_id = record[1]
        # self.gene_name = record[2]
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


class Ribo(object):

    def __init__(self, args):
        self.mrna_file = args.transcript
        self.mrna_seq = args.sequence
        self.mrna_dict = OrderedDict()
        self.longest = args.longest
        self.silence = args.silence

        self.bam = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]
        self.pysam_input = None

        self.psite = args.psite
        self.periodicity = args.periodicity
        self.min_length = args.min
        self.max_length = args.max
        self.offset = {}
        self.profile = None

        self.thread = args.thread

        self.total_rpf_df = None
        self.output = args.output

    def read_transcript(self):
        '''
        @Message  : import the gene sequence and message
        @Input    : mrna_seq --> the gene sequences
                    mrna_file --> the gene message
        @Return   : mrna_dict --> the gene message
        @Flow     : step1 --> import the mrna sequence generate from rpf_Reference script
                        contains: [5'-utr, cds, 3'-utr]
                    step2 --> import the gene message generate from rpf_Reference script
                        contains id same as the Mrna class
        '''

        # import the mrna sequences
        mrna_seq = SeqIO.parse(self.mrna_seq, "fasta")
        mrna_sequence = OrderedDict()
        for line in mrna_seq:
            # sys.stdout.writelines("import gene:  {gene}\n".format(gene=line.id))
            mrna_sequence[line.id] = line.seq

        # mrna_sequence = pysam.FastaFile(self.mrna_seq)

        # import the gene message
        with open(self.mrna_file, 'r') as trans_file_in:
            if self.longest:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    if int(record[6]) % 3 != 0 and not self.silence:
                        sys.stdout.write("{gene} CDS did not fit 3nt periodicity. \n".format(gene=record[2]))

                    if record[9] == "True":
                        now_mrna = Mrna(record)
                        # sys.stdout.writelines(now_mrna.transcript_id + '\n')
                        now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                        try:
                            now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                            now_mrna.rpf = np.zeros(len(now_mrna.seq))
                            # now_mrna.rpf = [0] * len(now_mrna.seq)
                            self.mrna_dict[now_mrna.transcript_id] = now_mrna
                        except KeyError:
                            continue
                    else:
                        continue
            else:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    if int(record[6]) % 3 != 0:
                        sys.stdout.write("{gene} CDS did not fit 3nt periodicity. \n".format(gene=record[2]))

                    now_mrna = Mrna(record)
                    # sys.stdout.writelines(now_mrna.transcript_id + '\n')
                    now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                    try:
                        now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                        now_mrna.rpf = np.zeros(len(now_mrna.seq))
                        self.mrna_dict[now_mrna.transcript_id] = now_mrna
                    except KeyError:
                        continue

    def read_offset(self):
        '''
        @Message  : import the offset and check the type of sequencing profile
        @Input    : psite --> the offset file
        @Return   : self.offset --> the offset dictionary
        @Flow     : step1 --> import the offset generate from the rpf_Offset script
                    step2 --> retrieve the psite form offset file
                    step3 --> check the type of sequencing profie from offset file
        '''

        offset = pd.read_csv(self.psite, sep='\t', header=0, names=None)
        offset.dropna(axis=0, inplace=True)
        offset["p_site"] = offset["p_site"].astype(int)
        offset["p_site"] = offset["p_site"] - 1

        raw_offset = offset[["length", "p_site"]].groupby(offset["length"])["p_site"].apply(list).T.to_dict()
        print("The import offset table: ", flush=True)
        print(raw_offset, flush=True)

        # filter the min_length, max_length and periodicity
        offset = offset[offset['periodicity'] >= self.periodicity]
        offset = offset[(offset["length"] >= self.min_length) & (offset["length"] <= self.max_length)]

        self.offset = offset[["length", "p_site"]].groupby(offset["length"])["p_site"].apply(list).T.to_dict()
        print("Filterd RPFs offset table: ", flush=True)
        print(self.offset, flush=True)

        # check the data type.
        if "third" in offset["ribo"].unique():
            self.profile = "trisome"
        elif "second" in offset["ribo"].unique():
            self.profile = "disome"
        elif "first" in offset["ribo"].unique():
            self.profile = "monosome"
        else:
            sys.stdout.writelines("Unknown type of sequence profile, please re-check the offset file!")
            sys.exit()

    def monosome_p_site(self):
        '''
        @Message  : retrieve the psite from each reads
        @Input    : pysam_input --> the mapping bam file interface
                    mrna_dict --> the gene message
        @Return   : output --> description
        @Flow     : step1 --> import the gene mapping bam file
                    step2 --> delete the reads not in the specified range
                    step3 --> retrieve the psite from each reads
        '''

        for line in self.pysam_input.fetch(until_eof=True):
            if line.reference_name in self.mrna_dict:
                map_start, map_end = line.get_blocks()[0]
                read_length = line.infer_read_length()

                if (self.min_length <= read_length <= self.max_length) and (read_length in self.offset):
                    if self.offset[read_length][0] == -1:
                        p_site = line.get_reference_positions()
                        self.mrna_dict[line.reference_name].rpf[p_site] += 1

                    else:
                        p_site = map_start + self.offset[read_length][0]
                        self.mrna_dict[line.reference_name].rpf[p_site] += 1

                else:
                    pass

    def disome_p_site(self):
        """
        like the mono_p_site,
        parser the reads with function twice,
        """

        for line in self.pysam_input.fetch(until_eof=True):
            if line.reference_name in self.mrna_dict:
                map_start, map_end = line.get_blocks()[0]
                read_length = line.infer_read_length()

                if (self.min_length <= read_length <= self.max_length) and (read_length in self.offset):

                    if self.offset[read_length][0] == -1:
                        p_site = map_start + np.arange(read_length)

                        for site in p_site:
                            try:
                                self.mrna_dict[line.reference_name].rpf[site] += 1
                            except KeyError:
                                continue
                            except IndexError:
                                continue
                    else:

                        try:
                            p_site = map_start + self.offset[read_length][0]
                            self.mrna_dict[line.reference_name].rpf[p_site] += 1
                            p_site = map_start + self.offset[read_length][1]
                            self.mrna_dict[line.reference_name].rpf[p_site] += 1
                        except KeyError:
                            # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                            continue
                        except IndexError:
                            continue
                else:
                    pass

    def trisome_p_site(self):
        """
        like the disome_p_site
        """

        for line in self.pysam_input.fetch(until_eof=True):
            if line.reference_name in self.mrna_dict:
                map_start, map_end = line.get_blocks()[0]
                read_length = line.infer_read_length()

                if (self.min_length <= read_length <= self.max_length) and (read_length in self.offset):
                    if self.offset[read_length][0] == -1:
                        p_site = map_start + np.arange(read_length)

                        for site in p_site:
                            try:
                                self.mrna_dict[line.reference_name].rpf[site] += 1
                            except KeyError:
                                continue
                            except IndexError:
                                continue
                    else:

                        try:
                            p_site = map_start + self.offset[read_length][0]
                            self.mrna_dict[line.reference_name].rpf[p_site] += 1
                            p_site = map_start + self.offset[read_length][1]
                            self.mrna_dict[line.reference_name].rpf[p_site] += 1
                            p_site = map_start + self.offset[read_length][2]
                            self.mrna_dict[line.reference_name].rpf[p_site] += 1
                        except KeyError:
                            # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                            continue
                        except IndexError:
                            continue
                else:
                    pass

    @staticmethod
    def flt_star_results(args):
        '''
        @Message  : filter the bam file generated by STAR
        @Input    : bam_in_file --> the path of input bam file
                    bam_out_file --> the path of output bam file
                    gene_list --> the list of gene name
                    tag --> the number of unique reads
                    sample_format --> the format of bam file
        @Return   : length_dict --> the dict of length distribution
                    bam_seq_dict --> the dict of sequence profile
        @Flow     : step1: open the bam file
                    step2: read the bam file and store the information into a dict
                    step3: calculate the length distribution
                    step4: output the filtered bam file
        '''
        
        bam_in_file, mrna_dict, min_length, max_length, offset, sample_format = args
        
        with pysam.AlignmentFile(bam_in_file, sample_format) as bam_in:
            for gene in mrna_dict.keys():
                try:
                    gene_reads = bam_in.fetch(gene)
                except ValueError:
                    continue

                # filter the mapped reads or not
                for reads in gene_reads:
                    if reads.is_unmapped:
                        continue

                    if reads.reference_name in mrna_dict:
                        map_start, map_end = reads.get_blocks()[0]
                        read_length = reads.infer_read_length()

                        if (min_length <= read_length <= max_length) and (read_length in offset):
                            if offset[read_length][0] == -1:
                                p_site = map_start + np.arange(read_length)
                                for site in p_site:
                                    try:
                                        mrna_dict[reads.reference_name].rpf[site] += 1
                                    except KeyError:
                                        continue
                                    except IndexError:
                                        continue
                            else:

                                p_site = map_start + offset[read_length][0]
                                try:
                                    mrna_dict[reads.reference_name].rpf[p_site] += 1

                                except KeyError:
                                    continue
                                    # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                                except IndexError:
                                    continue
                                    # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                        else:
                            pass

            bam_in.close()

        return mrna_dict
    
    def read_bam(self):
        '''
        @Message  : import the mapping bam file
        @Input    : bam --> the mapping bam file
        @Return   : 
        @Flow     : step1 --> run
        '''
        
        """
        1. import the mapping bam file
        2. retrieve the reads and add them to gene message dictionary

        """
        if self.sample_format.lower() == ".bam":
            sys.stdout.writelines("import file: {bam}.\n".format(bam=self.bam))
            read_format = 'rb'
        elif self.sample_format.lower() == ".sam":
            sys.stdout.writelines("import file: {bam}.\n".format(bam=self.bam))
            read_format = 'r'
        else:
            sys.stdout.writelines("Unknown file format, please input the correct bam or sam file.")
            sys.exit()

        self.pysam_input = pysam.AlignmentFile(self.bam, read_format)

        # calculate the length distribution and eliminate the reads aligned to the negative strand.
        if self.profile == "monosome":
            self.monosome_p_site()
            # # run the flt_star_results function with multi-thread
            # gene_list = list(self.mrna_dict.keys())
            # gene_list_split = np.array_split(gene_list, self.thread)
            # gene_list_split = [list(split) for split in gene_list_split]

            # args = [(self.bam, self.mrna_dict, self.min_length, self.max_length, self.offset, 'rb') for i in range(self.thread)]

            # # run the format_rpf_density function with multi-thread
            # from multiprocessing import Pool
            # pool = Pool(processes=self.thread)
            # rpf_list = pool.map(self.flt_star_results, args)
            # pool.close()
            # pool.join()

            # # merge the long data from each sub-list
            # self.mrna_dict = pd.concat(rpf_list, ignore_index=True)

        elif self.profile == "disome":
            self.disome_p_site()
        elif self.profile == "trisome":
            self.trisome_p_site()

        self.pysam_input.close()
    
    @staticmethod
    def format_rpf_density(args):
        '''
        @Message  : format the rpf density to long data
        @Input    : gene_list --> the gene list
                    mrna_dict --> the gene message
                    out_prefix --> the output prefix
        @Return   : the long data in pandas dataframe format
        @Flow     : step1 --> retrieve the gene message from the rpf density
                    step2 --> check the 3nt periodicity and trim the utr length
                    step3 --> format the gene message to the long data
                    step4 --> format the rpf density to three columns
                    step5 --> return the long data
        '''

        gene_list, mrna_dict, out_prefix = args
        rpf_list = []

        for name in gene_list:
            # retrieve the gene message from the rpf density
            isoform = mrna_dict[name]
            rpf_dict = OrderedDict()

            # trim the cds frame to fit the 3nt periodicity
            if isoform.cds_length % 3 != 0:
                shift_cds = isoform.cds_length % 3
                isoform.cds_length -= shift_cds
                isoform.utr3_length += shift_cds

            shift_5 = isoform.utr5_length % 3
            shift_3 = isoform.utr3_length % 3

            utr5_length = isoform.utr5_length - shift_5
            utr3_length = isoform.utr3_length - shift_3
            trim_length = isoform.length - shift_5 - shift_3
            
            # format the gene message to the long data
            rpf_dict['name'] = [name] * int(trim_length / 3)
            rpf_dict['now_nt'] = list(range(1 + shift_5, isoform.length - shift_3, 3))
            rpf_dict['from_tis'] = list(range(-utr5_length // 3, (trim_length - utr5_length) // 3))
            rpf_dict['from_tts'] = list(range((utr3_length - trim_length) // 3 + 1, utr3_length // 3 + 1))
            rpf_dict['region'] = ['5utr'] * (utr5_length // 3) + ['cds'] * (isoform.cds_length // 3) + ['3utr'] * (utr3_length // 3)

            if shift_3 != 0:
                trim_seq = isoform.seq[shift_5: -shift_3]
                rpf_dict['codon'] = [str(trim_seq[i:i + 3]) for i in range(0, trim_length, 3)]

                trim_rpf = isoform.rpf[shift_5: -shift_3]
                rpf_dict[out_prefix + '_f0'] = trim_rpf[0::3]
                rpf_dict[out_prefix + '_f1'] = trim_rpf[1::3]
                rpf_dict[out_prefix + '_f2'] = trim_rpf[2::3]
            elif shift_3 == 0:
                trim_seq = isoform.seq[shift_5::]
                rpf_dict['codon'] = [str(trim_seq[i:i + 3]) for i in range(0, trim_length, 3)]

                trim_rpf = isoform.rpf[shift_5::]
                rpf_dict[out_prefix + '_f0'] = trim_rpf[0::3]
                rpf_dict[out_prefix + '_f1'] = trim_rpf[1::3]
                rpf_dict[out_prefix + '_f2'] = trim_rpf[2::3]
            try:
                rpf = pd.DataFrame(rpf_dict)
                rpf_list.append(rpf)
            except ValueError:
                print('Checking the input genes: ' + str(name), flush=True)

        total_rpf_df = pd.concat(rpf_list, ignore_index=True)
        return total_rpf_df

    def run_format_rpf_with_multi_thread(self):
        '''
        @Message  : run the format_rpf_density function with multi-thread
        @Input    : self.mrna_dict --> the gene message
                    self.thread --> the number of thread
        @Return   : self.total_rpf_df --> the long data in pandas dataframe format
        @Flow     : step1 --> split the gene list into several sub-list
                    step2 --> run the format_rpf_density function with multi-thread
                    step3 --> merge the long data from each sub-list
        '''

        gene_list = list(self.mrna_dict.keys())

        gene_list_split = np.array_split(gene_list, self.thread)
        gene_list_split = [list(split) for split in gene_list_split]

        args = [(gene_list_split[i], self.mrna_dict, self.output) for i in range(self.thread)]

        # run the format_rpf_density function with multi-thread
        from multiprocessing import Pool
        pool = Pool(processes=self.thread)
        rpf_list = pool.map(self.format_rpf_density, args)
        pool.close()
        pool.join()

        # merge the long data from each sub-list
        self.total_rpf_df = pd.concat(rpf_list, ignore_index=True)
    
        numeric_cols = self.total_rpf_df.select_dtypes(include='number').columns
        self.total_rpf_df[numeric_cols] = self.total_rpf_df[numeric_cols].astype(int)

    def output_density(self):
        
        total_reads_pl = pl.DataFrame(self.total_rpf_df)
        total_reads_pl.write_csv(self.output + "_rpf.txt", separator = '\t', include_header = True)
