#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Quality.py


import os.path

import random
import sys
from collections import Counter
from collections import OrderedDict
from itertools import chain
from itertools import islice

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns

from multiprocessing import Pool

class Quality(object):

    def __init__(self, args):
        # opts for mrna file import
        self.longest = args.longest
        self.mrna_file = args.transcript
        self.mrna_dict = OrderedDict()
        self.length_dict = {}

        # opts for bam file import
        self.thread = args.thread
        self.pysam_input = None
        self.pysam_output = None
        self.sample_file = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]

        self.bam_seq_dict = {}
        self.sample_dict = None
        self.tag = args.tag
        self.reverse = args.reverse
        self.align = args.align

        # opts for bam length summary
        self.peak_length = 0
        self.peak_reads = 0
        self.mono = 1
        self.profile = None

        # opts for rpf saturation
        self.saturation_flag = args.saturation
        self.saturation = []
        self.x_ticks = np.array([i * 10 for i in range(1, 10)])

        # opts for file output
        self.output_prefix = args.output

    def read_transcript(self):
        '''
        @Message  : read the transcript file and store the information into a dict.
        @Input    : self.mrna_file --> the path of transcript file
        @Return   : self.mrna_dict --> the dict of transcript information
        @Flow     : step1: open the transcript file
                    step2: read the transcript file and store the information into a dict
        '''
        
        trans_file = self.mrna_file
        with open(trans_file, 'r') as trans_file_in:
            for line in islice(trans_file_in, 1, None):
                record = line.strip().split('\t')
                self.mrna_dict[record[2]] = [0] * 9
    
    def sort_index_bam(self):
        # check the data type, to support the SAM/BAM file format both.
        if self.sample_format.lower() == ".bam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            self.sample_format = 'rb'
        elif self.sample_format.lower() == ".sam":
            print("import file: {bam}.\n".format(bam=self.sample_file), flush=True)
            self.sample_format = 'r'
        else:
            print(
                "Unknown file format, please input the correct bam or sam file.", flush=True)
            sys.exit()
        
        # check the bam index file
        if os.path.exists(self.sample_file + ".bai"):
            pass
        else:
            print("index file: {bai} doesn't exist, create the index file.\n".format(bai=self.sample_file + ".bai"), flush=True)
            pysam.sort("-o", self.output_prefix + ".temp.sorted.bam", self.sample_file, "-@", str(self.thread))
            pysam.index(self.output_prefix + ".temp.sorted.bam", "-@", str(self.thread))
            self.sample_file = self.output_prefix + ".temp.sorted.bam"

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
        
        bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format = args

        tag_num = 0
        if tag == 0:
            tag_num = 99
        elif tag == 1:
            tag_num = 1
        
        bam_seq_dict = {}
        length_dict = {}
        
        with pysam.AlignmentFile(bam_in_file, sample_format) as bam_in:
            with pysam.AlignmentFile(bam_out_file, 'wb', template=bam_in) as bam_out:
                for gene in gene_list:
                    try:
                        gene_reads = bam_in.fetch(gene)
                    except ValueError:
                        continue

                    # filter the mapped reads or not
                    for reads in gene_reads:

                        # if reads.query_name == "A00202:1499:HWGTNDSX7:3:2518:22914:34194":
                        #     print(reads)
                        #     pass

                        if reads.is_unmapped:
                            continue

                        # filter the unique reads
                        if reads.has_tag('NH') and reads.get_tag('NH') <= tag_num:
                            read_length = reads.infer_read_length()
                            if reads.is_reverse:
                                try:
                                    length_dict[read_length][1] += 1
                                except KeyError:
                                    length_dict[read_length] = [0, 1]

                                if reverse_flag:
                                    try:
                                        bam_seq_dict[reads.query_name].add(reads.reference_name)
                                        bam_out.write(reads)
                                    except KeyError:
                                        bam_seq_dict[reads.query_name] = set()
                                        bam_seq_dict[reads.query_name].add(reads.reference_name)
                                        bam_out.write(reads)

                            else:
                                try:
                                    length_dict[read_length][0] += 1
                                except KeyError:
                                    length_dict[read_length] = [1, 0]
                                # mrna_dict[reads.reference_name].bam.append(reads)
                                try:
                                    bam_seq_dict[reads.query_name].add(reads.reference_name)
                                    bam_out.write(reads)
                                except KeyError:
                                    bam_seq_dict[reads.query_name] = set()
                                    bam_seq_dict[reads.query_name].add(reads.reference_name)
                                    bam_out.write(reads)
                bam_out.close()
            bam_in.close()

        return length_dict, bam_seq_dict
    
    @staticmethod
    def flt_hisat2_results(args):
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

        bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format = args

        tag_num = 0
        if tag == 0:
            zs_flg = True
            tag_num = len(gene_list)
        elif tag == 1:
            zs_flg = False
            tag_num = 1

        bam_seq_dict = {}
        length_dict = {}

        with pysam.AlignmentFile(bam_in_file, sample_format) as bam_in:
            with pysam.AlignmentFile(bam_out_file, 'wb', template=bam_in) as bam_out:
                for gene in gene_list:
                    try:
                        gene_reads = bam_in.fetch(gene)
                    except ValueError:
                        continue

                    # filter the mapped reads or not
                    for reads in gene_reads:
                        if reads.is_unmapped:
                            continue
                        # filter the unique reads
                        if reads.has_tag('NH') and reads.get_tag('NH') <= tag_num and reads.get_tag('ZS') == zs_flg:
                            read_length = reads.infer_read_length()
                            if reads.is_reverse:
                                try:
                                    length_dict[read_length][1] += 1
                                except KeyError:
                                    length_dict[read_length] = [0, 1]

                                if reverse_flag:
                                    try:
                                        bam_seq_dict[reads.query_name].add(reads.reference_name)
                                        bam_out.write(reads)
                                    except KeyError:
                                        bam_seq_dict[reads.query_name] = set()
                                        bam_seq_dict[reads.query_name].add(reads.reference_name)
                                        bam_out.write(reads)

                            else:
                                try:
                                    length_dict[read_length][0] += 1
                                except KeyError:
                                    length_dict[read_length] = [1, 0]
                                # mrna_dict[reads.reference_name].bam.append(reads)
                                try:
                                    bam_seq_dict[reads.query_name].add(reads.reference_name)
                                    bam_out.write(reads)
                                except KeyError:
                                    bam_seq_dict[reads.query_name] = set()
                                    bam_seq_dict[reads.query_name].add(reads.reference_name)
                                    bam_out.write(reads)

        return length_dict, bam_seq_dict
    
    @staticmethod
    def flt_bowtie2_results(args):
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

        bam_in_file, bam_out_file, gene_list, tag, reverse_flag, sample_format = args

        tag_num = 0
        if tag == 0:
            zs_flg = True
            tag_num = len(gene_list)
        elif tag == 1:
            zs_flg = False
            tag_num = 1

        bam_seq_dict = {}
        length_dict = {}

        with pysam.AlignmentFile(bam_in_file, sample_format) as bam_in:
            with pysam.AlignmentFile(bam_out_file, 'wb', template=bam_in) as bam_out:
                for gene in gene_list:
                    try:
                        gene_reads = bam_in.fetch(gene)
                    except ValueError:
                        continue
                    
                    # filter the mapped reads or not
                    for reads in gene_reads:
                        if reads.is_unmapped:
                            continue
                        # filter the unique reads
                        if reads.has_tag('AS') and reads.get_tag('AS') <= tag_num and reads.get_tag('XS') == zs_flg:
                            read_length = reads.infer_read_length()
                            if reads.is_reverse:
                                try:
                                    length_dict[read_length][1] += 1
                                except KeyError:
                                    length_dict[read_length] = [0, 1]

                                if reverse_flag:
                                    try:
                                        bam_seq_dict[reads.query_name].add(reads.reference_name)
                                        bam_out.write(reads)
                                    except KeyError:
                                        bam_seq_dict[reads.query_name] = set()
                                        bam_seq_dict[reads.query_name].add(reads.reference_name)
                                        bam_out.write(reads)

                            else:
                                try:
                                    length_dict[read_length][0] += 1
                                except KeyError:
                                    length_dict[read_length] = [1, 0]
                                # mrna_dict[reads.reference_name].bam.append(reads)
                                try:
                                    bam_seq_dict[reads.query_name].add(reads.reference_name)
                                    bam_out.write(reads)
                                except KeyError:
                                    bam_seq_dict[reads.query_name] = set()
                                    bam_seq_dict[reads.query_name].add(reads.reference_name)
                                    bam_out.write(reads)
        
        return length_dict, bam_seq_dict

    def fliter_mrna_reads(self):
        '''
        @Message  : retrieve the corrected rpf reads from bam file
        @Input    : self.sample_file --> the input bam file
                    self.sample_format --> the input bam file format
        @Return   : self.pysam_input --> the input bam file
                    self.pysam_output --> the output bam file
        @Flow     : step1 --> check the bam file format
                    step2 --> check the reads aligner algorithm of bam file
                    step3 --> import the bam file
                    step4 --> return the input and output bam file
        '''

        # flt by different algorithm
        pool = Pool(processes=self.thread)
        # splits = [gene_list[i:i + gene_list_len // self.thread] for i in range(0, gene_list_len, ceil(gene_list_len / self.thread))]

        splits = np.array_split(list(self.mrna_dict.keys()), self.thread)
        splits = [list(split) for split in splits]

        if self.align.upper() == 'HISAT2':
            args = [(self.sample_file, self.output_prefix + '_split_' + str(i) + '.bam', split, self.tag, self.reverse, self.sample_format) for i, split in enumerate(splits)]
            results = pool.map(self.flt_hisat2_results, args)

        elif self.align.upper() == 'STAR':
            args = [(self.sample_file, self.output_prefix + '_split_' + str(i) + '.bam', split, self.tag, self.reverse, self.sample_format) for i, split in enumerate(splits)]
            results = pool.map(self.flt_star_results, args)

        elif self.align.upper() == 'BOWTIE2' or self.align.upper() == 'BOWTIE':
            args = [(self.sample_file, self.output_prefix + '_split_' + str(i) + '.bam', split, self.tag, self.reverse, self.sample_format) for i, split in enumerate(splits)]
            results = pool.map(self.flt_bowtie2_results, args)
        else:
            print("Unknown align algorithm, please input the correct align algorithm.", flush=True)
            sys.exit()

        pool.close()
        pool.join()

        # merge length_dict with pandas
        # merge bam_seq_dict by dict keys and values
        first = True

        for length_dict, bam_seq_dict in results:
            if first:
                self.length_dict.update(length_dict)
                self.bam_seq_dict.update(bam_seq_dict)
                first = False
            else:
                for key, value in length_dict.items():
                    try:
                        self.length_dict[key][0] += value[0]
                        self.length_dict[key][1] += value[1]
                    except KeyError:
                        self.length_dict[key] = value
                        
                if self.saturation_flag:
                    self.bam_seq_dict = {**self.bam_seq_dict, **bam_seq_dict}
                    # for key, value in bam_seq_dict.items():
                    #     self.bam_seq_dict.setdefault(key, set()).update(value)
                else:
                    pass

        del results

    def merge_sort_index_bam(self):
        # merge the bam file
        bam_in_list = [self.output_prefix + "_split_" + str(i) + ".bam" for i in range(self.thread)]
        merge_parameters = ["-f", "-@", str(self.thread), self.output_prefix + ".bam"] + bam_in_list
        pysam.merge(*merge_parameters)

        # remove the temp bam file
        if os.path.exists(self.output_prefix + ".temp.sorted.bam"):
            os.remove(self.output_prefix + ".temp.sorted.bam")

        if os.path.exists(self.output_prefix + ".temp.sorted.bam.bai"):
            os.remove(self.output_prefix + ".temp.sorted.bam.bai")
        
        for i in range(self.thread):
            os.remove(self.output_prefix + "_split_" + str(i) + ".bam")

        # sort and index the bam file
        pysam.sort("-o", self.output_prefix + ".bam", self.output_prefix + ".bam", "-@", str(self.thread))
        pysam.index(self.output_prefix + ".bam", "-@", str(self.thread))

    def detect_seq_type(self):
        # auto-detect the type of sequence profile

        self.peak_length, self.peak_reads = max(self.length_dict.items(), key=lambda length: length[1])

        if 23 < self.peak_length < 35:
            print("{bam} is detected to be monosome-seq.\n".format(bam=self.sample_file), flush=True)
            self.profile = "monosome"
            self.mono = [19.5, 40.5]
        elif 53 < self.peak_length < 65:
            print("{bam} is detected to be disome-seq.\n".format(bam=self.sample_file), flush=True)
            self.profile = "disome"
            self.mono = [49.5, 70.5]
        elif 83 < self.peak_length < 95:
            print("{bam} is detected to be disome-seq.\n".format(bam=self.sample_file), flush=True)
            self.profile = "trisome"
            self.mono = [79.5, 100.5]
        elif 35 < self.peak_length < 53 or 65 < self.peak_length < 83 or 95 < self.peak_length:
            print("Warning! {bam} doesn't fit the empirical length distribution.!\n".format(bam=self.sample_file),
                  flush=True)
            print(''' 
            Monosome RPFs peak length is usually ~ 30 nt.
            Please check the files and run the detect_offset.py with specified peak_length!''',
                  flush=True)
            # sys.stdout.writelines('''
            # [25, 34] for the monosome, RPFs peak is usually ~ 30 nt.
            # [55, 65] for the disome, RPFs peak is usually ~ 60 nt.
            # [85, 95] for the trisome, RPFs peak is usually ~ 90 nt.
            # Please check the files and run the detect_offset.py with specified peak_length!''')

    def draw_the_length_distr(self, sorted_length):
        # output figure name
        out_pdf = self.output_prefix + "_length_distribution.pdf"
        out_png = self.output_prefix + "_length_distribution.png"

        matplotlib.use('AGG')

        # draw the figure
        length_df = pd.DataFrame(sorted_length).T
        fig = plt.figure(figsize=(6, 6), dpi=300)

        ax1 = fig.add_subplot(2, 1, 1)
        # ax1.bar(length_df.index, length_df[0], width=0.8, color="#FF9900")
        # sns.barplot(x=length_df.index, y=length_df[0], color="#f86934")
        sns.lineplot(x=length_df.index, y=length_df[0], color="#FF9900")
        ax1.set_xticks(length_df.index)
        ax1.set_xticklabels(length_df.index)
        ax1.set_xlim(self.mono)
        # ax1.set_ylim([0, round(length_df[0].max() + abs(length_df[0].max()) * 0.1)])
        ax1.set_ylabel('RPFs number')
        ax1.set_xlabel('RPFs length (nt)')
        ax1.set_title('plus strand')
        plt.xticks(rotation=90)

        # draw the minus strand RPFs length distribution
        ax2 = fig.add_subplot(2, 1, 2)
        # ax2.bar(length_df.index, length_df[1], width=0.8, color="#0099FF")
        # sns.barplot(x=length_df.index, y=length_df[0], color="#f86934")
        sns.lineplot(x=length_df.index, y=length_df[1], color="#0099FF")
        ax2.set_xticks(length_df.index)
        ax2.set_xticklabels(length_df.index)
        ax2.set_xlim(self.mono)
        # ax2.set_ylim([0, round(length_df[1].max() * 1.1)])
        ax2.set_ylabel('RPFs number')
        ax2.set_xlabel('RPFs length (nt)')
        ax2.set_title('Minus strand')

        plt.xticks(rotation=90)
        plt.suptitle('RPFs length distribution')
        plt.tight_layout()
        # plt.show()

        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()

    def write_length_distr(self):
        # write the length distribution.
        sorted_length = dict(sorted(self.length_dict.items(), key=lambda length: length[0], reverse=False))

        with open(self.output_prefix + "_length_distribution.txt", 'w') as length_out:
            length_out.writelines('\t'.join(["Length", "Plus", "Minus"]) + '\n')
            for reads_length, reads_num in sorted_length.items():
                length_out.writelines('\t'.join([str(reads_length), str(reads_num[0]), str(reads_num[1])]) + '\n')

        # draw the RPFs length distribution
        self.draw_the_length_distr(sorted_length)

    @staticmethod
    def sample_keys(step, keys):
        keys_list = list(keys)
        return random.sample(keys, step)
    
    def rpf_saturation(self):
        '''
        @Message  : calculate the RPFs saturation.
        @Input    : self.bam_seq_dict --> dict, {read_name: [gene1, gene2, ...]}
        @Return   : self.mrna_dict --> dict, {gene1: [site1, site2, ...]}
        @Flow     : step1 --> create the empty dict for resample steps and mrna count.
                    step2 --> create the function to resample the reads.
                    step3 --> runt the resample reads function.
                    step4 --> create the function to fill the empty mrna count dict.
                    step5 --> fill the empty mrna count dict.
        '''
        
        mapped_reads = len(self.bam_seq_dict)
        steps = [int(i * 0.05 * mapped_reads) for i in range(1, 10)]
        site = 0

        for step in steps:
            
            temp = random.sample(list(self.bam_seq_dict.keys()), step)
            temp_gene = list(chain.from_iterable(self.bam_seq_dict[i] for i in temp))

            temp_gene_dict = dict(Counter(temp_gene))

            gene_num = len(temp_gene_dict)
            self.saturation.append(gene_num)

            for gene, num in temp_gene_dict.items():
                self.mrna_dict[gene][site] = num

            site += 1

    @staticmethod
    def process_step(step_data):
        step, bam_seq_dict = step_data
        temp = random.sample(list(bam_seq_dict.keys()), step)
        temp_gene = []

        for i in temp:
            temp_gene.extend(bam_seq_dict[i])

        temp_gene_dict = dict(Counter(temp_gene))
        gene_num = len(temp_gene_dict)

        result = (gene_num, temp_gene_dict)
        return result

    def rpf_saturation_thread(self):
        from multiprocessing import Pool

        mapped_reads = len(self.bam_seq_dict)
        steps = [int(i * 0.05 * mapped_reads) for i in range(1, 10)]
        site = 0

        # Create a process pool
        pool = Pool(processes=self.thread)

        # Run the process_step function in parallel for each step
        step_data = [(step, self.bam_seq_dict) for step in steps]
        results = pool.map(self.process_step, step_data)

        for gene_num, temp_gene_dict in results:
            self.saturation.append(gene_num)

            for gene, num in temp_gene_dict.items():
                self.mrna_dict[gene][site] = num

            site += 1
            
    def draw_gene_saturation(self):
        out_pdf = self.output_prefix + "_gene_saturation.pdf"
        out_png = self.output_prefix + "_gene_saturation.png"
        out_gene = self.output_prefix + "_gene_saturation.txt"

        total_gene_num = len(self.mrna_dict)

        # convert the list to dataframe ,and set the title with 'gene'
        gene_df = pd.DataFrame(self.saturation + [total_gene_num], columns=['Count'])
        gene_df['Part'] = self.x_ticks.tolist() + [0]
        gene_df = gene_df[['Part', 'Count']]
        gene_df.to_csv(out_gene, sep='\t', index=False)

        matplotlib.use('AGG')

        fig = plt.figure(figsize=(8, 4), dpi=300)
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.bar(self.x_ticks, self.saturation, width=3, color="#FF9900")
        ax1.plot(self.x_ticks, self.saturation, linewidth=0.8)
        plt.xticks(rotation=90)
        ax1.set_xlabel('reads proportion (%)')
        ax1.set_ylabel('gene number')
        plt.title('covered gene saturation')

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.bar(self.x_ticks, [total_gene_num - i for i in self.saturation], width=3, color="#FF9900")
        ax2.plot(self.x_ticks, [total_gene_num - i for i in self.saturation], linewidth=0.8)
        plt.xticks(rotation=90)
        ax2.set_xlabel('reads proportion (%)')
        ax2.set_ylabel('gene number')
        plt.title('uncovered gene saturation')

        plt.tight_layout()
        # plt.show()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()

    def draw_rpf_saturation(self):
        # output figure name
        out_pdf = self.output_prefix + "_reads_saturation.pdf"
        out_png = self.output_prefix + "_reads_saturation.png"
        out_rpf = self.output_prefix + "_reads_saturation.txt"
        

        mrna_df = pd.DataFrame(self.mrna_dict).T
        mrna_df.loc[:, 'mean'] = mrna_df.mean(axis=1)
        mrna_df = mrna_df.sort_values(['mean'], ascending=True)
        mrna_df.to_csv(out_rpf, sep='\t')

        quantile_list = mrna_df['mean'].quantile([0.25, 0.5, 0.75])
        mrna_0_25 = mrna_df[mrna_df['mean'] < quantile_list.iloc[0]]
        mrna_25_50 = mrna_df[(quantile_list.iloc[0] <= mrna_df['mean']) & (mrna_df['mean'] < quantile_list.iloc[1])]
        mrna_50_75 = mrna_df[(quantile_list.iloc[1] <= mrna_df['mean']) & (mrna_df['mean'] < quantile_list.iloc[2])]
        mrna_75_100 = mrna_df[quantile_list.iloc[2] <= mrna_df['mean']]

        # draw the figure
        matplotlib.use('AGG')
        colors = ['#ffc7ce', '#c6efce', '#ffeb9c', '#9cc3e5']
        flierprops = dict(marker='o', markersize=2, markerfacecolor='#c95859', markeredgecolor='none')
        fig = plt.figure(figsize=(17, 4), dpi=300)

        fig.add_subplot(1, 4, 1)
        ax1 = sns.boxplot(data=mrna_0_25.iloc[:, 0:9], color='#30adfe', flierprops=flierprops)
        ax1.set(xlabel='reads proportion (%)', ylabel='Reads count')
        ax1.set_yscale('log')
        ax1.set_xticklabels(self.x_ticks, rotation=90)

        fig.add_subplot(1, 4, 2)
        ax2 = sns.boxplot(data=mrna_25_50.iloc[:, 0:9], color='#30adfe', flierprops=flierprops)
        ax2.set(xlabel='reads proportion (%)', ylabel='Reads count')
        ax2.set_yscale('log')
        ax2.set_xticklabels(self.x_ticks, rotation=90)

        fig.add_subplot(1, 4, 3)
        ax3 = sns.boxplot(data=mrna_50_75.iloc[:, 0:9], color='#30adfe', flierprops=flierprops)
        ax3.set(xlabel='reads proportion (%)', ylabel='Reads count')
        ax3.set_yscale('log')
        ax3.set_xticklabels(self.x_ticks, rotation=90)

        fig.add_subplot(1, 4, 4)
        ax4 = sns.boxplot(data=mrna_75_100.iloc[:, 0:9], color='#30adfe', flierprops=flierprops)
        ax4.set(xlabel='reads proportion (%)', ylabel='Reads count')
        ax4.set_yscale('log')
        ax4.set_xticklabels(self.x_ticks, rotation=90)

        plt.suptitle('Reads saturation')
        plt.tight_layout()
        # plt.show()
        fig.savefig(fname=out_pdf)
        fig.savefig(fname=out_png)
        plt.close()