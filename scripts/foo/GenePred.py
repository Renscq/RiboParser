#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : test01.py
# @Script  : GenePred.py


import os
import stat
import subprocess
import sys
from bisect import bisect_left
from bisect import bisect_right
from collections import OrderedDict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Seq


class Chrom(object):

    def __init__(self, line):
        self.attr = line.description
        self.chromosome = str(line.id)
        self.start = 1
        self.end = len(line.seq)
        self.seq = line.seq


class GenePred(object):
    def __init__(self, args):
        self.gene = OrderedDict()
        self.gtf = args.gtf
        self.gtf_format = os.path.splitext(self.gtf)[-1]
        self.gp_df = None
        self.columns = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
                        'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat',
                        'IstringexonFrames']
        self.coding = args.coding
        self.longest = args.longest
        self.utr = args.utr

        # genome file
        self.fasta = args.genome
        self.chroms_dict = OrderedDict()
        self.mrna_seq = OrderedDict()
        self.cds_seq = OrderedDict()

        # output file
        self.out_prefix = args.output
        self.whole = args.whole
        self.gp_file = self.out_prefix + '.genepred'
        self.gtf_new = self.out_prefix + '.norm.gtf'
        self.whole_file = self.out_prefix + '.whole.txt'
        self.txt_file = self.out_prefix + '.norm.txt'
        self.rna_file = self.out_prefix + '.norm.rna.fa'
        self.cds_file = self.out_prefix + '.norm.cds.fa'

    def gtf2gp(self):
        '''
        @Message  : function to convert the gtf to genepred.
        @Input    : self.gtf_format --> format of the gtf file
                    self.gtf --> gtf file
                    self.gp_file --> genepred file

        @Return   : output --> description
        @Flow     : step1 --> check the gff3ToGenePred and gtfToGenePred in the current environment
                    step2 --> convert the gtf or gff to genepred
        '''
        
        script_path = os.path.split((os.path.abspath(sys.argv[0])))[0]
        nowpath = os.path.dirname(script_path)
        if self.gtf_format in ['.gff', '.gff3']:
            gff3ToGenePred = subprocess.run(['which', 'gff3ToGenePred'], capture_output=True, text=True)
            if gff3ToGenePred.returncode == 0:
                run_tmp = subprocess.run(['gff3ToGenePred', '-warnAndContinue', '-rnaNameAttr=attr',
                                        '-geneNameAttr=attr', self.gtf, self.gp_file])
            else:
                print('{0} not found in current environment.'.format('gff3ToGenePred'), flush=True)
                sys.exit(1)

        elif self.gtf_format == '.gtf':
            gtfToGenePred = subprocess.run(['which', 'gtfToGenePred'], capture_output=True, text=True)
            if gtfToGenePred.returncode == 0:
                run_tmp = subprocess.run(['gtfToGenePred', '-allErrors', '-genePredExt',
                                        '-ignoreGroupsWithoutExons', self.gtf, self.gp_file])
            else:
                print('{0} not found in current environment.'.format('gtfToGenePred'), flush=True)
                sys.exit(1)
        else:
            print('Unknown file format {0}.'.format(self.gtf_format),)
            print('Please input annotation file in GTF or GFF format.', flush=True)
            sys.exit(1)

    def gp2gtf(self):
        '''
        @Message  : function to convert the genepred to gtf file.
        @Input    : self.gp_file --> genepred file
                    self.gtf_new --> new gtf file
        @Return   : output --> genePred file
        @Flow     : step1 --> check the genePredToGtf in the current environment
                    step2 --> convert the genepred to gtf
        '''
        
        # script_path = os.path.split((os.path.abspath(sys.argv[0])))[0]
        # nowpath = os.path.dirname(script_path)

        genePredToGtf = subprocess.run(['which', 'genePredToGtf'], capture_output=True, text=True)

        if genePredToGtf.returncode == 0:
            run_tmp = subprocess.run(['genePredToGtf', 'file', '-utr', '-honorCdsStat',
                                      '-source=ribo', self.gp_file, self.gtf_new])
        else:
            print('{0} not found in current environment.'.format('genePredToGtf'), flush=True)
            sys.exit(1)
        # os.remove(self.gp_file)

    def read_genepred(self):
        '''
        @Message  : function import the parser the genepred file .
        @Input    : self.gp_file --> gene pred file derived from gtf file
        @Return   : self.gp_df --> dataframe of the genepred file
        @Flow     : step1 --> read the genepred file and rename the columns
                    step2 --> calculate the utr5, cds and utr3 length
        '''

        self.gp_df = pd.read_csv(self.gp_file, sep='\t', header=None)
        os.remove(self.gp_file)
        self.gp_df.columns = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
                              'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat',
                              'IstringexonFrames']

        utr5_list = []
        cds_list = []
        utr3_list = []
        print('Import gtf annotation.', flush=True)

        for idx, rows in self.gp_df.iterrows():

            if idx % 1000 == 0:
                print('rows: ' + str(idx), flush=True)

            exonStarts = np.array(rows['exonStarts'].split(',')[:-1], dtype='int')
            exonEnds = np.array(rows['exonEnds'].split(',')[:-1], dtype='int')

            cdsStart = rows['cdsStart']
            cdsEnd = rows['cdsEnd']

            if rows['cdsStartStat'] == 'none' or rows['cdsEndStat'] == 'none':
                utr5_list.append(0)
                cds_list.append(0)
                utr3_list.append(0)

            elif rows['exonCount'] == 1:
                if rows.strand == '+':
                    utr5_list.append(int(abs(cdsStart - exonStarts)))
                    cds_list.append(int(abs(cdsEnd - cdsStart)))
                    utr3_list.append(int(abs(exonEnds - cdsEnd)))

                elif rows.strand == '-':
                    utr5_list.append(int(abs(exonEnds - cdsEnd)))
                    cds_list.append(int(abs(cdsEnd - cdsStart)))
                    utr3_list.append(int(abs(cdsStart - exonStarts)))

            elif rows['exonCount'] > 1:
                if rows.strand == '+':
                    intronStarts = exonStarts[1:]
                    intronEnds = exonEnds[:-1]
                    introns = intronStarts - intronEnds
                    intron_length = abs(sum(introns))

                    # get utr5 length
                    now_right = bisect_right(exonStarts, cdsStart) - 1
                    utr5_intron = sum(introns[0:now_right])
                    utr5_length = cdsStart - exonStarts[0] - utr5_intron

                    # get utr3 length
                    now_left = bisect_left(exonEnds, cdsEnd)
                    utr3_intron = sum(introns[now_left:])
                    utr3_length = exonEnds[-1] - cdsEnd - utr3_intron

                    # get cds length
                    cds_length = abs(cdsEnd - cdsStart) - intron_length + utr5_intron + utr3_intron

                    utr5_list.append(utr5_length)
                    cds_list.append(cds_length)
                    utr3_list.append(utr3_length)

                elif rows.strand == '-':
                    intronStarts = exonEnds[:-1]
                    intronEnds = exonStarts[1:]
                    introns = intronEnds - intronStarts
                    intron_length = abs(sum(introns))

                    # get utr5 length
                    now_right = bisect_right(exonStarts, cdsStart) - 1
                    utr3_intron = sum(introns[0:now_right])
                    utr3_length = cdsStart - exonStarts[0] - utr3_intron

                    # get utr3 length
                    now_left = bisect_left(exonEnds, cdsEnd)
                    utr5_intron = sum(introns[now_left:])
                    utr5_length = exonEnds[-1] - cdsEnd - utr5_intron
                    
                    # get cds length
                    cds_length = abs(cdsEnd - cdsStart) - intron_length + utr5_intron + utr3_intron

                    utr5_list.append(utr5_length)
                    cds_list.append(cds_length)
                    utr3_list.append(utr3_length)

            else:
                print('Error: gene range was wrong at :' + rows.name + '\n', flush=True)
                continue
        
        print('rows: ' + str(idx), flush=True)

        self.gp_df['utr5_length'] = utr5_list
        self.gp_df['cds_length'] = cds_list
        self.gp_df['utr3_length'] = utr3_list

    def get_rep_transcript_bak(self):
        '''
        @Message  : function to filter the representative transcripts.
        @Input    : self.gp_df --> dataframe of the genepred file
        @Return   : self.gp_df --> dataframe of the genepred file with annotated column [rep_transcript]
        @Flow     : step1 --> group the dataframe by gene name
                    step2 --> label the max cds length
                    step3 --> label the max utr5 + cds + utr3 length
                    step4 --> label the first transcript
                    step5 --> add the flag of representative transcript
        @Note     : low efficiency of this mothod to label the representative transcript.
        '''
        
        print('Filter representative transcripts.', flush=True)

        counter = 0
        def add_true_flag(group):

            nonlocal counter
            counter += 1
            if counter % 500 == 0:
                print('Representative gene: ' + str(counter), flush=True)

            if len(group) == 1:
                group['rep_transcript'] = True
                return group
            
            else:
                max_cds_length = group['cds_length'].max()
                max_utr5_cds_utr3 = group.loc[group['cds_length'] == max_cds_length, ['utr5_length', 'cds_length', 'utr3_length']].sum(axis=1).max()

                max_cds_index = group.loc[group['cds_length'] == max_cds_length].index
                max_utr5_cds_utr3_index = group.loc[group['utr5_length'] + group['cds_length'] + group['utr3_length'] == max_utr5_cds_utr3].index

                if len(max_cds_index) == 1:
                    group.loc[max_cds_index, 'rep_transcript'] = True
                elif len(max_utr5_cds_utr3_index) == 1:
                    group.loc[max_utr5_cds_utr3_index, 'rep_transcript'] = True
                else:
                    group.loc[max_cds_index[0], 'rep_transcript'] = True
                return group

        self.gp_df['rep_transcript'] = False
        self.gp_df = self.gp_df.groupby('name2').apply(add_true_flag)
        print('Representative gene: ' + str(counter), flush=True)

        # print(self.gp_df)

    def get_rep_transcript(self):
        '''
        @Message  : function to filter the representative transcripts.
        @Input    : self.gp_df --> dataframe of the genepred file
        @Return   : self.gp_df --> dataframe of the genepred file with annotated column [rep_transcript]
        @Flow     : step1 --> calculate the gene length, length = cds_length + utr5_length + utr3_length
                    step2 --> sort the dataframe by gene name, cds length, gene length, utr5 length and utr3 length
                    step3 --> reset the index of the dataframe
                    step4 --> group the dataframe by gene name and label the max cds length and max gene length
                    step5 --> label the first transcript
                    step6 --> reset the index of the dataframe
        '''

        print('Filter representative transcripts.', flush=True)

        self.gp_df['gene_length'] = self.gp_df['cds_length'] + self.gp_df['utr5_length'] + self.gp_df['utr3_length']
        self.gp_df = self.gp_df.sort_values(by=['name2', 'cds_length', 'gene_length', 'utr5_length', 'utr3_length'], 
                                            ascending=[False, False, False, False, False])

        self.gp_df.reset_index(drop=True, inplace=True)
        self.gp_df['rep_transcript'] = False
        self.gp_df.loc[self.gp_df.groupby('name2').head(1).index, 'rep_transcript'] = True
        self.gp_df.sort_values(by=['chrom', 'txStart', 'name2', 'rep_transcript'], 
                               ascending=[True, True, True, False],
                               inplace=True)
        self.gp_df.reset_index(drop=True, inplace=True)

    def add_utr(self):
        '''
        @Message  : function to add the pseudo utr.
        @Input    : self.gp_df --> dataframe of the genepred file
        @Return   : self.gp_df --> dataframe of the genepred file with annotated column [modified]
        @Flow     : step1 --> create the list to store the new exon starts and ends
                    step2 --> skip the transcript without cds
                    step3 --> calculate the new exon start site to extend the utr5
                    step4 --> calculate the new exon end site to extend the utr3
                    step5 --> add the pseudo utr to the transcript
        '''
        
        exonStarts_list, exonEnds_list = [], []
        txStart_list, txEnd_list = [], []
        utr5_list, utr3_list = [], []

        self.gp_df.loc[:, 'modified'] = 'False'

        # add pseudo utr
        if self.utr == 0:
            pass

        else:
            print('Add pseudo UTR.', flush=True)

            for idx, rows in self.gp_df.iterrows():
                if idx % 1000 == 0:
                    print('rows: ' + str(idx), flush=True)

                # skip the transcript without cds
                if rows['cdsStartStat'] == 'none' or rows['cdsEndStat'] == 'none':
                    exonStarts_list.append(rows.exonStarts)
                    exonEnds_list.append(rows.exonEnds)
                    txStart_list.append(rows.txStart)
                    txEnd_list.append(rows.txEnd)
                    utr5_list.append(rows.utr5_length)
                    utr3_list.append(rows.utr3_length)
                    continue

                # check exon starts for each genes
                exonStarts = rows['exonStarts'].split(',')[:-1]
                
                if int(rows['utr5_length']) >= self.utr:
                    # skip the transcript with utr5 length >= self.utr
                    exonStarts_list.append(rows.exonStarts)
                    txStart_list.append(rows.txStart)
                    utr5_list.append(rows.utr5_length)

                else:
                    # calculate the new exon start site to extend the utr5
                    new_start = int(exonStarts[0]) - self.utr

                    if new_start < 0:
                        # set the negative new start site to 0, since the chromosome start site can not be negative
                        exonStarts[0] = '0'
                        exonStarts_list.append(','.join(exonStarts) + ',')
                        txStart_list.append(0)
                        utr5_list.append(rows['cdsStart'])
                        self.gp_df.loc[idx, 'modified'] = 'True'
                    else:
                        # set the new start site to new_start
                        exonStarts[0] = str(new_start)
                        exonStarts_list.append(','.join(exonStarts) + ',')
                        txStart_list.append(rows.txStart - self.utr)
                        utr5_list.append(rows.utr5_length + self.utr)
                        self.gp_df.loc[idx, 'modified'] = 'True'

                # check exon ends for each genes
                exonEnds = rows['exonEnds'].split(',')[:-1]

                if int(rows['utr3_length']) >= self.utr:
                    # skip the transcript with utr3 length >= self.utr
                    exonEnds_list.append(rows.exonEnds)
                    txEnd_list.append(rows.txEnd)
                    utr3_list.append(rows.utr3_length)
                else:
                    # calculate the new exon end site to extend the utr3
                    new_end = int(exonEnds[-1]) + self.utr

                    if new_end > self.chroms_dict[str(rows['chrom'])].end:
                        # set the new end site to chromosome end site, since the chromosome end site can not be larger than the chromosome end site
                        exonEnds[-1] = str(self.chroms_dict[str(rows['chrom'])].end)
                        exonEnds_list.append(','.join(exonEnds) + ',')
                        txEnd_list.append(self.chroms_dict[str(rows['chrom'])].end)
                        utr3_list.append(rows.utr3_length + self.chroms_dict[str(rows['chrom'])].end - rows['cdsEnd'])
                        self.gp_df.loc[idx, 'modified'] = 'True'
                    else:
                        # set the new end site to new_end
                        exonEnds[-1] = str(new_end)
                        exonEnds_list.append(','.join(exonEnds) + ',')
                        txEnd_list.append(new_end)
                        utr3_list.append(rows.utr3_length + self.utr)
                        self.gp_df.loc[idx, 'modified'] = 'True'

            self.gp_df['exonStarts'] = exonStarts_list
            self.gp_df['exonEnds'] = exonEnds_list
            self.gp_df['txStart'] = txStart_list
            self.gp_df['txEnd'] = txEnd_list
            self.gp_df['utr5_length'] = utr5_list
            self.gp_df['utr3_length'] = utr3_list

    def read_genome(self):
        '''
        @Message  : function import the genome sequences.
        @Input    : self.fasta --> genome file
        @Return   : self.chroms_dict --> dictionary of the chromosome sequences
        @Flow     : step1 --> read the genome file with SeqIO
                    step2 --> store the chromosome sequences in the dictionary
        '''
        
        record = SeqIO.parse(self.fasta, "fasta")

        for line in record:
            print("import chromosome: {chrom}".format(chrom=line.id), flush=True)
            self.chroms_dict[str(line.id)] = Chrom(line)


    def get_seq(self):
        '''
        @Message  : function to retrieve the mRNA sequences.
        @Input    : self.gp_df --> dataframe of the genepred file
        @Return   : self.mrna_seq --> dictionary of the mRNA sequences
        @Flow     : step1 --> iterate the dataframe to get the mRNA sequences
                    step2 --> store the mRNA sequences in the dictionary
                    step3 --> store the cds sequences in the dictionary
        '''

        for idx, rows in self.gp_df.iterrows():

            # skip the transcript without cds
            if rows['cds_length'] == 0:
                continue

            chromesome = str(rows['chrom'])
            rna = list(zip(rows.exonStarts.split(',')[:-1], rows.exonEnds.split(',')[:-1]))
            tmp_rna_seq = Seq.Seq('')

            if rows['strand'] == "-":
                for exon in rna:
                    exon_start, exon_end = int(exon[0]), int(exon[1])
                    tmp_rna_seq += self.chroms_dict[chromesome].seq[exon_start: exon_end]
                tmp_rna_seq = tmp_rna_seq.reverse_complement()

            else:
                for exon in rna:
                    exon_start, exon_end = int(exon[0]), int(exon[1])
                    tmp_rna_seq += self.chroms_dict[chromesome].seq[exon_start: exon_end]

            self.mrna_seq[rows['name']] = tmp_rna_seq
        
            self.cds_seq[rows['name']] = tmp_rna_seq[rows['utr5_length']: rows['utr5_length'] + rows['cds_length']]

    def write_txt(self):
        '''
        @Message  : function to output the txt file.
        @Input    : self.gp_df --> dataframe of the genepred file
        @Return   : gp_df_new --> dataframe of the genepred file
                    txt_df --> dataframe of the txt file
        @Flow     : step1 --> output the whole genepred file
                    step2 --> filter the coding transcript
                    step3 --> filter the longest transcript
                    step4 --> output the genePred file
                    step5 --> format the genePred file to txt file
        '''
        
        if self.whole:
            self.gp_df.to_csv(self.whole_file, sep='\t', header=True, index=None)

        if self.coding:
            self.gp_df = self.gp_df.loc[self.gp_df['cds_length'] != 0, ]

        if self.longest:
            self.gp_df = self.gp_df.loc[self.gp_df['rep_transcript'] == 'True', ]

        self.gp_df.loc[:, 'name'] = self.gp_df.loc[:, 'name'].str.replace('cds-', '').str.replace('rna-', '')
        self.gp_df.loc[:, 'name2'] = self.gp_df.loc[:, 'name2'].str.replace('gene-', '')

        # output the genePred file
        gp_df_new = self.gp_df.loc[:, self.columns].copy()
        gp_df_new.to_csv(self.gp_file, sep='\t', header=None, index=None)

        # output the txt file
        txt_columns = ['chrom', 'name2', 'name', 'txStart', 'txEnd', 'utr5_length', 'cds_length', 'utr3_length',
                       'strand', 'rep_transcript', 'modified']
        txt_df = self.gp_df.loc[:, txt_columns].copy()

        txt_df.columns = ['chromosome', 'gene_id', 'transcript_id', 'start', 'end', 'utr5_length',
                          'cds_length', 'utr3_length', 'strand', 'rep_transcript', 'modified']
        
        txt_df.loc[:, ['utr5_length', 'cds_length', 'utr3_length']] = txt_df.loc[:, ['utr5_length', 'cds_length', 'utr3_length']].astype(int)
        txt_df = txt_df.loc[txt_df['cds_length'] != 0, ]
        txt_df.to_csv(self.txt_file, sep='\t', header=True, index=None)

    def write_seq(self):
        '''
        @Message  : function to output the mRNA and cds sequences.
        @Input    : self.rna_file --> mrna sequences file
                    self.cds_file --> cds sequences file
        '''
        
        with open(self.rna_file, 'w') as seq_out:
            for rna, mess in self.mrna_seq.items():
                seq_out.writelines('\n'.join(['>' + rna, str(mess)]) + '\n')

        with open(self.cds_file, 'w') as cds_out:
            for cds, mess in self.cds_seq.items():
                cds_out.writelines('\n'.join(['>' + cds, str(mess)]) + '\n')
