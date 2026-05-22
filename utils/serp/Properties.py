#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/11/3 9:50
# @Project : riboParser
# @Script  : Properties.py


import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import ProtParam


class SeqProperties(object):

    def __init__(self, seq):
        '''
        @Message  : function for sequence properties calculation
        @Input    : param --> properties
        @Return   : output --> sequence propeerties
        @Flow     : step1 --> run ProteinAnalysis
        '''
        
        self.prpt = ProtParam.ProteinAnalysis(seq)
        self.gravy = str(round(self.prpt.gravy(), 3))
        if self.prpt.flexibility():
            self.flexibility = str(round(sum(self.prpt.flexibility()) / len(self.prpt.flexibility()), 3))
        else:
            self.flexibility = 'nan'
        self.instability = str(round(self.prpt.instability_index(), 3))
        self.isoelectric_point = str(round(self.prpt.isoelectric_point(), 3))
        self.structure = self.prpt.secondary_structure_fraction()
        self.helix = str(round(self.structure[0], 3))
        self.turn = str(round(self.structure[1], 3))
        self.sheet = str(round(self.structure[2], 3))


class Sequence(object):
    def __init__(self, args):
        self.fasta = args.fasta

        self.codon_dict = {
            "AAA": ["Lys", "K"],
            "AAC": ["Asn", "N"],
            "AAG": ["Lys", "K"],
            "AAT": ["Asn", "N"],
            "ACA": ["Thr", "T"],
            "ACC": ["Thr", "T"],
            "ACG": ["Thr", "T"],
            "ACT": ["Thr", "T"],
            "AGA": ["Arg", "R"],
            "AGC": ["Ser", "S"],
            "AGG": ["Arg", "R"],
            "AGT": ["Ser", "S"],
            "ATA": ["Ile", "I"],
            "ATC": ["Ile", "I"],
            "ATG": ["Met", "M"],
            "ATT": ["Ile", "I"],
            "CAA": ["Gln", "Q"],
            "CAC": ["HIS", "H"],
            "CAG": ["Gln", "Q"],
            "CAT": ["HIS", "H"],
            "CCA": ["Pro", "P"],
            "CCC": ["Pro", "P"],
            "CCG": ["Pro", "P"],
            "CCT": ["Pro", "P"],
            "CGA": ["Arg", "R"],
            "CGC": ["Arg", "R"],
            "CGG": ["Arg", "R"],
            "CGT": ["Arg", "R"],
            "CTA": ["Leu", "L"],
            "CTC": ["Leu", "L"],
            "CTG": ["Leu", "L"],
            "CTT": ["Leu", "L"],
            "GAA": ["Glu", "E"],
            "GAC": ["Asp", "D"],
            "GAG": ["Glu", "E"],
            "GAT": ["Asp", "D"],
            "GCA": ["Ala", "A"],
            "GCC": ["Ala", "A"],
            "GCG": ["Ala", "A"],
            "GCT": ["Ala", "A"],
            "GGA": ["Gly", "G"],
            "GGC": ["Gly", "G"],
            "GGG": ["Gly", "G"],
            "GGT": ["Gly", "G"],
            "GTA": ["Val", "V"],
            "GTC": ["Val", "V"],
            "GTG": ["Val", "V"],
            "GTT": ["Val", "V"],
            "TAA": ["Stop", "*"],
            "TAC": ["Tyr", "Y"],
            "TAG": ["Stop", "*"],
            "TAT": ["Tyr", "Y"],
            "TCA": ["Ser", "S"],
            "TCC": ["Ser", "S"],
            "TCG": ["Ser", "S"],
            "TCT": ["Ser", "S"],
            "TGA": ["Stop", "*"],
            "TGC": ["Cys", "C"],
            "TGG": ["Trp", "W"],
            "TGT": ["Cys", "C"],
            "TTA": ["Leu", "L"],
            "TTC": ["Phe", "F"],
            "TTG": ["Leu", "L"],
            "TTT": ["Phe", "F"],
        }
        
        self.codon_anno = pd.DataFrame.from_dict(self.codon_dict).T
        self.codon_anno = self.codon_anno.reset_index()
        self.codon_anno.columns = ['Codon', 'AA', 'Abbr.']

        self.cds_data_frame = ''

        # opts for file output
        self.output_prefix = args.output
        self.output_gene_cu = self.output_prefix + '.Gene.CodonUsage.xlsx'
        self.output_whole_cu = self.output_prefix + '.Whole.CodonUsage.txt'
        self.output_seq = self.output_prefix + '.Properties.txt'


    def count_codons_numpy(self, cds_sequence):
        '''
        @Message  : function for codon count summary.
        @Input    : param --> cds_sequence
        @Return   : output --> codon and count
        @Flow     : step1 --> split the cds sequence to codon sequence
                    step2 --> summary the codon count
        '''
        
        # ensure the length of the sequence is a multiple of 3
        if len(cds_sequence) % 3 != 0:
            raise ValueError("CDS length is not a multiple of 3")
        
        # split the sequence into codons
        codons = np.array([cds_sequence[i:i+3] for i in range(0, len(cds_sequence), 3)])
        unique, counts = np.unique(codons, return_counts=True)
        return dict(zip(unique, counts))
    
    def create_codon_table(self):
        '''
        @Message  : create the codon and amino acid table.
        @Input    : param --> self.fasta
        @Return   : output --> self.cds_data_frame
        @Flow     : step1 --> summary the codon and count of cds sequence
                    step2 --> merge all gene codon count to one data table
        '''
        
        codon_dict = OrderedDict()

        for record in SeqIO.parse(self.fasta, "fasta"):
            # convert the lower case to upper case
            record.seq = record.seq.upper()
            codon_count = self.count_codons_numpy(str(record.seq))
            codon_count['Length'] = len(record.seq) / 3
            codon_dict[record.id] = codon_count

        # convert codon count to codon table
        codon_table = pd.DataFrame(codon_dict).T.fillna(0)
        codon_table = codon_table.reset_index(names=['Gene'])

        # convert codon table to long format
        self.cds_data_frame = pd.melt(codon_table, id_vars=['Gene', 'Length'], 
                              value_vars=codon_table.columns.difference(['Gene', 'Length']),
                              var_name='Codon', value_name='Count')
        
        # merge codon annotation and sort the table
        self.cds_data_frame = self.cds_data_frame.merge(self.codon_anno, on='Codon', how='left')
        self.cds_data_frame = self.cds_data_frame.sort_values(by=['Gene', 'Abbr.']).reset_index(drop=True)


    def calc_gene_codon_usage(self):
        '''
        @Message  : function for gene codon usage calculation
        @Input    : param --> self.cds_data_frame
        @Return   : output --> Frequency, RSCU, CAI
        @Flow     : step1 --> calculate the codon frequency
                    step2 --> calculate the RSCU and CAI
                    step3 --> fill NA with 0 and round the value
                    step4 --> convert the data frame to wide format
        '''
        
        # calculate codon frequency
        self.cds_data_frame['Frequency'] = self.cds_data_frame['Count'] / self.cds_data_frame['Length'] * 1000

        # calculate RSCU and CAI
        self.cds_data_frame['RSCU'] = self.cds_data_frame.groupby(['Gene', 'AA'])['Frequency'].transform(lambda x: x / x.mean())
        self.cds_data_frame['CAI'] = self.cds_data_frame.groupby(['Gene', 'AA'])['Frequency'].transform(lambda x: x / x.max())

        # fill NA with 0
        self.cds_data_frame.fillna({'Frequency': 0, 'RSCU': 0, 'CAI': 0}, inplace=True)

        # round the value
        self.cds_data_frame['Frequency'] = self.cds_data_frame['Frequency'].round(4)
        self.cds_data_frame['RSCU'] = self.cds_data_frame['RSCU'].round(4)
        self.cds_data_frame['CAI'] = self.cds_data_frame['CAI'].round(4)

        # convert the data frame to wide format
        cds_freq_data_frame_wide = self.cds_data_frame.pivot_table(index=['Gene', 'Length'], columns='Codon', values='Frequency').reset_index()
        cds_rscu_data_frame_wide = self.cds_data_frame.pivot_table(index=['Gene', 'Length'], columns='Codon', values='RSCU').reset_index()
        cds_cai_data_frame_wide = self.cds_data_frame.pivot_table(index=['Gene', 'Length'], columns='Codon', values='CAI').reset_index()

        # output the three data frame
        cds_freq_data_frame_wide.to_csv(self.output_prefix + '_frequency.txt', sep='\t', index=False)
        cds_rscu_data_frame_wide.to_csv(self.output_prefix + '_rscu.txt', sep='\t', index=False)
        cds_cai_data_frame_wide.to_csv(self.output_prefix + '_cai.txt', sep='\t', index=False)

        # with pd.ExcelWriter(self.output_gene_cu) as writer:
        #     cds_freq_data_frame_wide.to_excel(writer, sheet_name='Frequency', index=False)
        #     cds_rscu_data_frame_wide.to_excel(writer, sheet_name='RSCU', index=False)
        #     cds_cai_data_frame_wide.to_excel(writer, sheet_name='CAI', index=False)

    def calc_whole_codon_usage(self):
        '''
        @Message  : function for whole codon usage calculation
        @Input    : param --> self.cds_data_frame
        @Return   : output --> Frequency, RSCU, CAI
        @Flow     : step1 --> calculate the codon frequency
                    step2 --> calculate the RSCU and CAI
                    step3 --> fill NA with 0 and round the value
                    step4 --> save the data frame to file
        '''
        
        cds_data_frame_sum = self.cds_data_frame.groupby(['Codon', 'AA', 'Abbr.'])['Count'].sum().reset_index()
        cds_data_frame_sum['Frequency'] = cds_data_frame_sum['Count'] / cds_data_frame_sum['Count'].sum() * 1000
        cds_data_frame_sum['RSCU'] = cds_data_frame_sum['Frequency'] / cds_data_frame_sum.groupby('AA')['Frequency'].transform('mean')
        cds_data_frame_sum['CAI'] = cds_data_frame_sum['Frequency'] / cds_data_frame_sum.groupby('AA')['Frequency'].transform('max')
        
        cds_data_frame_sum.fillna({'Frequency': 0, 'RSCU': 0, 'CAI': 0}, inplace=True)

        cds_data_frame_sum['Frequency'] = cds_data_frame_sum['Frequency'].round(4)
        cds_data_frame_sum['RSCU'] = cds_data_frame_sum['RSCU'].round(4)
        cds_data_frame_sum['CAI'] = cds_data_frame_sum['CAI'].round(4)

        cds_data_frame_sum.to_csv(self.output_prefix + '_whole_codon_usage.txt', sep='\t', index=False)

    def protein_analysis(self):
        '''
        @Message  : function for protein properties calculation
        @Input    : param --> self.fasta
        @Return   : output --> sequence properties
        @Flow     : step1 --> read the fasta file
                    step2 --> write the sequence properties to file
        '''
        
        with open(self.output_seq, 'w') as out_seq:
            out_seq.writelines('\t'.join(['ID', 'Seq', 'Length', 'Gravy', 'Flexibility', 'Instability',
                                          'Isoelectric_Point', 'Helix', 'Turn', 'Sheet']) + '\n')
            fa_seq = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))

            for now_idx, now_seq in fa_seq.items():
                now_aa = str(now_seq.seq.translate())
                now_len = str(len(now_aa))

                if now_aa[-1] == '*':
                    now_aa = now_aa[:-1]
                    now_len = str(len(now_aa))

                if '*' in now_aa:
                    sys.stdout.writelines(
                        'Skip sequence contain stop codon: {id}\t{aa}'.format(id=now_idx, aa=now_aa) + '\n')
                    continue

                prpt = SeqProperties(now_aa)
                out_seq.writelines('\t'.join([now_idx, now_aa, now_len, prpt.gravy, prpt.flexibility, prpt.instability,
                                              prpt.isoelectric_point, prpt.helix, prpt.turn, prpt.sheet]) + '\n')
