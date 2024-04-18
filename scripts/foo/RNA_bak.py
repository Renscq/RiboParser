#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : RNA.py


import os
import sys
from collections import OrderedDict
from itertools import islice
from Bio import SeqIO
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
        self.reads = []
        self.seq = None


class RNA(object):

    def __init__(self, args):
        self.mrna_file = args.transcript
        self.mrna_seq = args.sequence
        self.mrna_dict = OrderedDict()
        self.longest = args.longest

        self.bam = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]
        self.pysam_input = None

        self.psite = args.psite
        self.min_length = args.min
        self.max_length = args.max

        self.offset = {}
        self.profile = None
        self.rolling = args.rolling
        self.pair_end = args.pair_end

        self.thread = args.thread

        self.output = args.output

    def read_transcript(self):
        '''
        @Message  : read the transcript file
        @Input    : self.mrna_file --> mrna message table generated by riboParser
                    self.mrna_seq --> mrna sequence file
        @Return   : self.mrna_dict --> mrna dict contain the mrna message
        @Flow     : step1 --> read the mrna sequence file
                    step2 --> read the mrna message file
                    step3 --> filter the longest transcript
                    
        '''
        
        # read the mrna sequence file
        record = SeqIO.parse(self.mrna_seq, "fasta")
        mrna_sequence = OrderedDict()
        for line in record:
            # sys.stdout.writelines("import gene:  {gene}\n".format(gene=line.id))
            mrna_sequence[line.id] = line.seq

        # read the mrna message file
        with open(self.mrna_file, 'r') as trans_file_in:
            if self.longest:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    # if int(record[7]) % 3 != 0:
                    #     sys.stdout.write("Delete {gene}. CDS length doesn't fit 3nt periodicity. \n".format(gene=record[3]))

                    if record[9] == "True":
                        now_mrna = Mrna(record)
                        # sys.stdout.writelines(now_mrna.transcript_id + '\n')
                        now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                        try:
                            now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                            now_mrna.reads = [0] * len(now_mrna.seq)
                            self.mrna_dict[now_mrna.transcript_id] = now_mrna
                        except KeyError:
                            continue
                    else:
                        continue
            else:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    # if int(record[7]) % 3 != 0:
                    #     sys.stdout.write("Delete {gene}. CDS length doesn't fit 3nt periodicity. \n".format(gene=record[3]))
                    
                    now_mrna = Mrna(record)
                    # sys.stdout.writelines(now_mrna.transcript_id + '\n')
                    now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                    try:
                        now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                        now_mrna.reads = [0] * len(now_mrna.seq)
                        self.mrna_dict[now_mrna.transcript_id] = now_mrna
                    except KeyError:
                        continue

    def read_offset(self):
        '''
        @Message  : read the offset file
        @Input    : self.psite --> input offset file name
        @Return   : self.offset --> offset dictionary
        @Flow     : step1 --> read the offset file
                    step2 --> drop the null value
                    step3 --> change the offset type to int
        '''
        
        offset = pd.read_csv(self.psite, sep='\t', header=0, names=None)
        offset.dropna(axis=0, inplace=True)
        offset["p_site"] = offset["p_site"].astype(int)
        offset["p_site"] = offset["p_site"] - 1
        self.offset = offset[["length", "p_site"]].groupby(offset["length"])["p_site"].apply(list).T.to_dict()
        print(self.offset)

    def read_bam(self):
        '''
        @Message  : read the bam file
        @Input    : self.bam --> input bam file name
        @Return   : self.pysam_input --> pysam.AlignmentFile
        @Flow     : step1 --> check the input file format
                    step2 --> read the bam file
                    step3 --> check the bam index
                    step4 --> build the bam index
        '''
        
        # check the input file format
        if self.sample_format.lower() == ".bam":
            sys.stdout.writelines("import file: {bam}.\n".format(bam=self.bam))
            self.read_format = 'rb'
        elif self.sample_format.lower() == ".sam":
            sys.stdout.writelines("import file: {bam}.\n".format(bam=self.bam))
            self.read_format = 'r'
        else:
            sys.stdout.writelines("Unknown file format, please input the correct bam or sam file.")
            sys.exit()

        self.pysam_input = pysam.AlignmentFile(self.bam, self.read_format)
        
        # check the bam index
        if self.pysam_input.has_index():
            pass
        else:
            import subprocess
            print("build the bam index for : {bam}".format(bam=self.bam), flush=True)

            subprocess.call("samtools index {bam}".format(bam=self.bam), shell=True)

            self.pysam_input = pysam.AlignmentFile(self.bam, self.read_format)

    def mono_density(self):
        '''
        @Message  : calculate the reads monosome density
        @Input    : self.mrna_dict --> mrna dict contain the mrna message
                    self.pysam_input --> pysam.AlignmentFile
        @Return   : self.mrna_dict --> mrna dict contain the mrna message and reads density
        @Flow     : step1 --> delete the unmapped reads
                    step2 --> delete the reverse reads
                    step3 --> filter the reads length
                    step4 --> count the reads density
        '''

        # fetch the reads from bam file
        for line in self.pysam_input.fetch(until_eof=True):
            # delete the unmapped reads
            if line.is_unmapped:
                continue

            # delete the reverse reads
            if line.is_reverse and not self.pair_end:
                continue

            # count the reads density
            if line.reference_name in self.mrna_dict:
                map_start, map_end = line.get_blocks()[0]
                read_length = line.infer_read_length()
                if self.min_length <= read_length <= self.max_length:
                    if self.offset[read_length][0] == -1:
                        p_site = map_start + np.arange(read_length)

                        for site in p_site:
                            try:
                                self.mrna_dict[line.reference_name].reads[site] += 1
                            except KeyError:
                                continue
                            except IndexError:
                                continue
                    else:
                        try:
                            p_site = map_start + self.offset[read_length][0]
                            self.mrna_dict[line.reference_name].reads[p_site] += 1
                        except KeyError:
                            continue
                            # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                        except IndexError:
                            continue
                            # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                else:
                    pass
        self.pysam_input.close()

    def roll_density(self):
        '''
        @Message  : calculate the reads rolling polysome density with multi-thread
        @Input    : self.mrna_dict --> mrna dict contain the mrna message
                    self.pysam_input --> pysam.AlignmentFile
        @Return   : self.mrna_dict --> mrna dict contain the mrna message and reads density
        @Flow     : step1 --> delete the unmapped reads
                    step2 --> delete the reverse reads
                    step3 --> filter the reads length
                    step4 --> count the reads density with rolling polysome
        '''
        
        # fetch the reads from bam file
        for line in self.pysam_input.fetch(until_eof=True):
            # delete the unmapped reads
            if line.is_unmapped:
                continue

            # delete the reverse reads
            if line.is_reverse and not self.pair_end:
                continue

            map_start, map_end = line.get_blocks()[0]
            read_length = line.infer_read_length()
            
            if self.min_length <= read_length <= self.max_length:

                if self.offset[read_length][0] == -1:
                    p_site = map_start + np.arange(read_length)

                    for site in p_site:
                        try:
                            self.mrna_dict[line.reference_name].reads[site] += 1
                        except KeyError:
                            continue
                        except IndexError:
                            continue
                else:

                    try:
                        offset = self.offset[read_length][0]
                        p_site = map_start + self.offset[read_length][0]
                        while offset <= read_length:
                            self.mrna_dict[line.reference_name].reads[p_site] += 1
                            offset += 32
                            p_site += 32
                    except KeyError:
                        continue
                        # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                    except IndexError:
                        continue
                        # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
            else:
                pass
        
        self.pysam_input.close()

    def calculate_density(self):
        '''
        @Message  : function for calculate the reads density
        @Input    : self.rolling --> bool, True for rolling polysome, False for monosome
        @Return   : self.mrna_dict --> mrna dict contain the mrna message and reads density
        @Flow     : 
        '''

        # run the mono_density or roll_density function
        if self.rolling:
            self.roll_density()
        else:
            self.mono_density()

    @staticmethod
    def format_reads(args):
        '''
        @Message  : function for format the reads density
        @Input    : self.mrna_dict --> mrna dict contain the mrna message and reads density
        @Return   : self.total_reads_df --> dataframe contain the reads density
        @Flow     : step1 --> trim the cds frame to fit the 3nt periodicity
                    step2 --> trim the reads density
                    step3 --> format the reads density
                    step4 --> merge the reads density
        '''
        
        gene_list, mrna_dict, out_prefix = args
        reads_list = []

        for name in gene_list:
            # retrieve the gene message from the reads density
            isoform = mrna_dict[name]
            reads_dict = OrderedDict()

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

            # trim the reads density
            reads_dict['name'] = [name] * int(trim_length / 3)
            reads_dict['now_nt'] = list(range(1 + shift_5, isoform.length - shift_3, 3))
            reads_dict['from_tis'] = list(range(-utr5_length // 3, (trim_length - utr5_length) // 3))
            reads_dict['from_tts'] = list(range((utr3_length - trim_length) // 3 + 1, utr3_length // 3 + 1))
            reads_dict['region'] = ['5utr'] * (utr5_length // 3) + ['cds'] * (isoform.cds_length // 3) + ['3utr'] * (utr3_length // 3)

            # trim the split gene nucleotide sequence to codon
            if shift_3 != 0:
                trim_seq = isoform.seq[shift_5: -shift_3]
                reads_dict['codon'] = [str(trim_seq[i:i + 3]) for i in range(0, trim_length, 3)]

                trim_reads = isoform.reads[shift_5: -shift_3]
                reads_dict[out_prefix + '_f0'] = trim_reads[0::3]
                reads_dict[out_prefix + '_f1'] = trim_reads[1::3]
                reads_dict[out_prefix + '_f2'] = trim_reads[2::3]
            elif shift_3 == 0:
                trim_seq = isoform.seq[shift_5::]
                reads_dict['codon'] = [str(trim_seq[i:i + 3]) for i in range(0, trim_length, 3)]

                trim_reads = isoform.reads[shift_5::]
                reads_dict[out_prefix + '_f0'] = trim_reads[0::3]
                reads_dict[out_prefix + '_f1'] = trim_reads[1::3]
                reads_dict[out_prefix + '_f2'] = trim_reads[2::3]
            try:
                reads = pd.DataFrame(reads_dict)
                reads_list.append(reads)
            except ValueError:
                print('Skip the genes: ' + str(name), flush=True)

        total_reads_df = pd.concat(reads_list, ignore_index=True)
        return total_reads_df

    def run_format_reads_with_multi_thread(self):
        """
        1. split the gene list into several sub-list
        2. run the format_rreads_density function with multi-thread
        3. merge the long data from each sub-list

        """
        
        # split the gene list into several sub-list
        gene_list = list(self.mrna_dict.keys())
        gene_list_len = len(gene_list)

        # gene_list_split = [gene_list[i:i + gene_list_len // self.thread] for i in range(0, gene_list_len, gene_list_len // self.thread)]
        gene_list_split = np.array_split(gene_list, self.thread)
        gene_list_split = [list(split) for split in gene_list_split]

        args = [(gene_list_split[i], self.mrna_dict, self.output) for i in range(self.thread)]

        # run the format_rreads_density function with multi-thread
        from multiprocessing import Pool
        pool = Pool(processes=self.thread)
        results = pool.map(self.format_reads, args)
        pool.close()
        pool.join()

        # merge the long data from each sub-list
        self.total_reads_df = pd.concat(results, ignore_index=True)

    def output_density(self):
        
        # total_reads_df["sum"] = total_readsf_df[['frame_1', 'frame_2', 'frame_3']].sum(axis=1)
        # total_reads_df.to_csv(rna_attr.output + "_rna.txt", sep = '\t', index = False)
        total_reads_pl = pl.DataFrame(self.total_reads_df)
        total_reads_pl.write_csv(self.output + "_rna.txt", separator = '\t', has_header = True)