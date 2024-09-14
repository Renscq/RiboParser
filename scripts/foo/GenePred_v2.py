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
        self.chromosome = line.id
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
        # output file
        self.out_prefix = args.output
        self.gp_file = self.out_prefix + '.genepred'
        self.gtf_new = self.out_prefix + '.norm.gtf'
        self.txt_file = self.out_prefix + '.norm.txt'
        self.seq_file = self.out_prefix + '.norm.fa'

    def gtf2gp(self):
        script_path = os.path.split((os.path.abspath(sys.argv[0])))[0]
        nowpath = os.path.dirname(script_path)
        if self.gtf_format in ['.gff', '.gff3']:
            if not os.path.exists(nowpath + '/ucsctools/gff3ToGenePred'):
                print('{0} not found in package.'.format('gff3ToGenePred'), flush=True)
                sys.exit(1)
            run_tmp = subprocess.run([nowpath + '/ucsctools/gff3ToGenePred', '-warnAndContinue', '-rnaNameAttr=attr',
                                      '-geneNameAttr=attr', self.gtf, self.gp_file])
        elif self.gtf_format == '.gtf':
            if not os.path.exists(nowpath + '/ucsctools/gtfToGenePred'):
                print('{0} not found in package.'.format('gtfToGenePred'), flush=True)
                sys.exit(1)
            run_tmp = subprocess.run([nowpath + '/ucsctools/gtfToGenePred', '-allErrors', '-genePredExt',
                                      '-ignoreGroupsWithoutExons', self.gtf, self.gp_file])
        else:
            print('Unknown file format {0}.'.format(self.gtf_format),)
            print('Please input annotation file in GTF or GFF format.', flush=True)
            sys.exit(1)

    def gp2gtf(self):
        script_path = os.path.split((os.path.abspath(sys.argv[0])))[0]
        nowpath = os.path.dirname(script_path)

        if not os.path.exists(nowpath + '/ucsctools/genePredToGtf'):
            print('{0} not found in package.'.format('genePredToGtf'), flush=True)
            sys.exit(1)
        else:
            # os.chmod(nowpath + '../ucsctools/genePredToGtf', stat.S_IXOTH)
            run_tmp = subprocess.run([nowpath + '/ucsctools/genePredToGtf', 'file', '-utr', '-honorCdsStat',
                                      '-source=ribo', self.gp_file, self.gtf_new])
        # os.remove(self.gp_file)

    def read_genepred(self):
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
            if rows['name'] == 'rna-XM_016990516.1':
                print('')
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

                    utr5_list.append(abs(cdsStart - exonStarts))
                    cds_list.append(abs(cdsEnd - cdsStart))
                    utr3_list.append(abs(exonEnds - cdsEnd))
                elif rows.strand == '-':
                    utr5_list.append(abs(exonEnds - cdsEnd))
                    cds_list.append(abs(cdsEnd - cdsStart))
                    utr3_list.append(abs(cdsStart - exonStarts))

            elif rows['exonCount'] > 1:
                if rows.strand == '+':
                    # summary 5utr intron, 3utr intron, cds intron
                    intronStarts = exonStarts[1:]
                    intronEnds = exonEnds[:-1]

                    utr5_Intron_Start = intronStarts[intronStarts <= cdsStart]
                    utr5_Intron_End = intronEnds[intronEnds <= cdsStart]

                    utr3_Intron_Start = intronStarts[intronStarts >= cdsEnd]
                    utr3_Intron_End = intronEnds[intronEnds >= cdsEnd]

                    cds_Intron_Start = intronStarts[intronStarts >= cdsStart]
                    cds_Intron_Start = cds_Intron_Start[cds_Intron_Start <= cdsEnd]
                    cds_Intron_End = intronEnds[intronEnds >= cdsStart]
                    cds_Intron_End = cds_Intron_End[cds_Intron_End <= cdsEnd]

                    if utr5_Intron_Start.any() and utr5_Intron_End.any():
                        utr5_Intron = abs(sum(utr5_Intron_Start - utr5_Intron_End))
                    else:
                        utr5_Intron = 0
                    if utr3_Intron_Start.any() and utr3_Intron_End.any():
                        utr3_Intron = abs(sum(utr5_Intron_Start - utr5_Intron_End))
                    else:
                        utr3_Intron = 0
                    if cds_Intron_Start.any() and cds_Intron_End.any():
                        cds_Intron = abs(sum(cds_Intron_Start - cds_Intron_End))
                    else:
                        cds_Intron = 0

                    introns = intronStarts - intronEnds
                    intron_length = utr5_Intron + utr3_Intron + cds_Intron

                    # get utr5 length
                    now_right = bisect_right(exonStarts, cdsStart) - 1
                    utr5_length = cdsStart - exonStarts[0] - sum(introns[0:now_right])
                    # get utr3 length
                    now_left = bisect_left(exonEnds, cdsEnd)
                    utr3_length = exonEnds[-1] - cdsEnd - sum(introns[now_left:])
                    # get cds length
                    cds_length = abs(cdsEnd - cdsStart) - intron_length

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
                    utr3_length = cdsStart - exonStarts[0] - sum(introns[0:now_right])
                    # get utr3 length
                    now_left = bisect_left(exonEnds, cdsEnd)
                    utr5_length = exonEnds[-1] - cdsEnd - sum(introns[now_left:])
                    # get cds length
                    cds_length = abs(cdsEnd - cdsStart) - intron_length

                    utr5_list.append(utr5_length)
                    cds_list.append(cds_length)
                    utr3_list.append(utr3_length)

            else:
                print('Error: gene range was wrong at :' + rows.name + '\n', flush=True)
                continue
        self.gp_df['utr5_length'] = utr5_list
        self.gp_df['cds_length'] = cds_list
        self.gp_df['utr3_length'] = utr3_list

    def get_rep_transcript(self):

        print('Filter representative transcripts.', flush=True)
        gp_df_max = self.gp_df.groupby(['name2'])[['name', 'cds_length']].max()
        id_dict = dict(zip(gp_df_max.name.values, gp_df_max.name.values))
        flags = self.gp_df.loc[:, 'name'].apply(lambda x: 'True' if x in id_dict else 'False')
        self.gp_df.loc[:, 'rep_transcript'] = flags

    def add_utr(self):
        exonStarts_list, exonEnds_list = [], []
        txStart_list, txEnd_list = [], []
        utr5_list, utr3_list = [], []

        self.gp_df.loc[:, 'modified'] = 'False'
        if self.utr == 0:
            pass
        else:
            print('Add pseudo utr.', flush=True)

            for idx, rows in self.gp_df.iterrows():
                if idx % 1000 == 0:
                    print('rows: ' + str(idx), flush=True)

                if rows['cdsStartStat'] == 'none' or rows['cdsEndStat'] == 'none':
                    exonStarts_list.append(rows.exonStarts)
                    exonEnds_list.append(rows.exonEnds)
                    txStart_list.append(rows.txStart)
                    txEnd_list.append(rows.txEnd)
                    utr5_list.append(rows.utr5_length)
                    utr3_list.append(rows.utr3_length)
                    continue

                # exon starts
                exonStarts = rows['exonStarts'].split(',')[:-1]
                if int(rows['utr5_length']) >= self.utr:
                    exonStarts_list.append(rows.exonStarts)
                    txStart_list.append(rows.txStart)
                    utr5_list.append(rows.utr5_length)

                else:
                    new_start = int(exonStarts[0]) - self.utr
                    if new_start < 0:
                        exonStarts[0] = '0'
                        exonStarts_list.append(','.join(exonStarts) + ',')
                        txStart_list.append(0)
                        utr5_list.append(rows['cdsStart'])
                        self.gp_df.loc[idx, 'modified'] = 'True'
                    else:
                        exonStarts[0] = str(new_start)
                        exonStarts_list.append(','.join(exonStarts) + ',')
                        txStart_list.append(rows.txStart - self.utr)
                        utr5_list.append(rows.utr5_length + self.utr)
                        self.gp_df.loc[idx, 'modified'] = 'True'

                # exon ends
                exonEnds = rows['exonEnds'].split(',')[:-1]
                if int(rows['utr3_length']) >= self.utr:
                    exonEnds_list.append(rows.exonEnds)
                    txEnd_list.append(rows.txEnd)
                    utr3_list.append(rows.utr3_length)
                else:
                    new_end = int(exonEnds[-1]) + self.utr
                    if new_end > self.chroms_dict[rows['chrom']].end:
                        exonEnds[-1] = str(self.chroms_dict[rows['chrom']].end)
                        exonEnds_list.append(','.join(exonEnds) + ',')
                        txEnd_list.append(self.chroms_dict[rows['chrom']].end)
                        utr3_list.append(rows.utr3_length + self.chroms_dict[rows['chrom']].end - rows['cdsEnd'])
                        self.gp_df.loc[idx, 'modified'] = 'True'
                    else:
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

        record = SeqIO.parse(self.fasta, "fasta")
        for line in record:
            print("import chromosome: {chrom}".format(chrom=line.id), flush=True)
            self.chroms_dict[str(line.id)] = Chrom(line)

    def get_seq(self):
        self.mrna_seq = OrderedDict()

        for idx, rows in self.gp_df.iterrows():
            if rows['cds_length'] == 0:
                continue

            chromesome = str(rows['chrom'])
            rna = list(zip(rows.exonStarts.split(',')[:-1], rows.exonEnds.split(',')[:-1]))
            tmp_seq = Seq.Seq('')

            if rows['strand'] == "-":
                for exon in rna:
                    exon_start, exon_end = int(exon[0]), int(exon[1])
                    tmp_seq += self.chroms_dict[chromesome].seq[exon_start: exon_end]
                tmp_seq = tmp_seq.reverse_complement()

            else:
                for exon in rna:
                    exon_start, exon_end = int(exon[0]), int(exon[1])
                    tmp_seq += self.chroms_dict[chromesome].seq[exon_start: exon_end]

            self.mrna_seq[rows['name']] = tmp_seq

    def write_txt(self):
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
        txt_df.insert(loc=2, column='gene_name', value='')
        txt_df.columns = ['chromosome', 'gene_id', 'gene_name', 'transcript_id', 'start', 'end', 'utr5_length',
                          'cds_length', 'utr3_length', 'strand', 'rep_transcript', 'modified']
        txt_df.loc[:, ['utr5_length', 'cds_length', 'utr3_length']] = txt_df.loc[:, ['utr5_length', 'cds_length',
                                                                                     'utr3_length']].astype(int)
        txt_df = txt_df.loc[txt_df['cds_length'] != 0, ]
        txt_df.to_csv(self.txt_file, sep='\t', header=True, index=None)

    def write_seq(self):
        with open(self.seq_file, 'w') as seq_out:
            for rna, mess in self.mrna_seq.items():
                seq_out.writelines('\n'.join(['>' + rna, str(mess)]) + '\n')
