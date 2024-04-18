#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @time    : 2021/11/3 9:50
# @Project : riboParser
# @Script  : Properties.py


import sys
from collections import OrderedDict

import CAI
from Bio import SeqIO
from Bio.SeqUtils import CodonUsage
from Bio.SeqUtils import ProtParam


class SeqProperties(object):
    def __init__(self, seq):
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
        self.codon_anno = {
            'AAA': ['Lys', 'K', 2], 'AAC': ['Asn', 'N', 2], 'AAG': ['Lys', 'K', 2], 'AAT': ['Asn', 'N', 2],
            'ACA': ['Thr', 'T', 4], 'ACC': ['Thr', 'T', 4], 'ACG': ['Thr', 'T', 4], 'ACT': ['Thr', 'T', 4],
            'AGA': ['Arg', 'R', 6], 'AGC': ['Ser', 'S', 6], 'AGG': ['Arg', 'R', 6], 'AGT': ['Ser', 'S', 6],
            'ATA': ['Ile', 'I', 3], 'ATC': ['Ile', 'I', 3], 'ATG': ['Met', 'M', 1], 'ATT': ['Ile', 'I', 3],
            'CAA': ['Gln', 'Q', 2], 'CAC': ['HIS', 'H', 2], 'CAG': ['Gln', 'Q', 2], 'CAT': ['HIS', 'H', 2],
            'CCA': ['Pro', 'P', 4], 'CCC': ['Pro', 'P', 4], 'CCG': ['Pro', 'P', 4], 'CCT': ['Pro', 'P', 4],
            'CGA': ['Arg', 'R', 6], 'CGC': ['Arg', 'R', 6], 'CGG': ['Arg', 'R', 6], 'CGT': ['Arg', 'R', 6],
            'CTA': ['Leu', 'L', 6], 'CTC': ['Leu', 'L', 6], 'CTG': ['Leu', 'L', 6], 'CTT': ['Leu', 'L', 6],
            'GAA': ['Glu', 'E', 2], 'GAC': ['Asp', 'D', 2], 'GAG': ['Glu', 'E', 2], 'GAT': ['Asp', 'D', 2],
            'GCA': ['Ala', 'A', 4], 'GCC': ['Ala', 'A', 4], 'GCG': ['Ala', 'A', 4], 'GCT': ['Ala', 'A', 4],
            'GGA': ['Gly', 'G', 4], 'GGC': ['Gly', 'G', 4], 'GGG': ['Gly', 'G', 4], 'GGT': ['Gly', 'G', 4],
            'GTA': ['Val', 'V', 4], 'GTC': ['Val', 'V', 4], 'GTG': ['Val', 'V', 4], 'GTT': ['Val', 'V', 4],
            'TAA': ['*', '*', 3], 'TAC': ['Tyr', 'Y', 2], 'TAG': ['*', '*', 3], 'TAT': ['Tyr', 'Y', 2],
            'TCA': ['Ser', 'S', 6], 'TCC': ['Ser', 'S', 6], 'TCG': ['Ser', 'S', 6], 'TCT': ['Ser', 'S', 6],
            'TGA': ['*', '*', 3], 'TGC': ['Cys', 'C', 2], 'TGG': ['Trp', 'W'], 'TGT': ['Cys', 'C', 2],
            'TTA': ['Leu', 'L', 6], 'TTC': ['Phe', 'F', 2], 'TTG': ['Leu', 'L', 6], 'TTT': ['Phe', 'F', 2]
        }
        self.codon_anno_div = {}
        self.seq_dict = OrderedDict()

        # opts for file output
        self.output_prefix = args.output
        self.output_aa = self.output_prefix + '.CodonUsage.txt'
        self.output_seq = self.output_prefix + '.Properties.txt'

    def calc_cai(self):
        cu = CodonUsage.CodonAdaptationIndex()
        try:
            cu.generate_index(self.fasta)
        except ZeroDivisionError:

            print('ZeroDivisionError', flush=True)
            print('Sequence file is not suitable for CAI calculation.', flush=True)
            print('Some codon is empty', flush=True)
            sys.exit()
        cai_dict = cu.index
        frequency_dict = cu.codon_count
        total_aa = sum(cu.codon_count.values())

        fa_seq = [rec.seq for rec in SeqIO.parse(self.fasta, 'fasta')]
        rscu_dict = CAI.RSCU(fa_seq)

        for codon, aa in self.codon_anno.items():
            # CAI
            if codon in cai_dict.keys():
                self.codon_anno[codon].append(str(round(cai_dict[codon], 3)))
            else:
                self.codon_anno[codon].append('nan')
            # RSCU
            if codon in rscu_dict.keys():
                self.codon_anno[codon].append(str(round(rscu_dict[codon] / aa[2], 3)))
            else:
                self.codon_anno[codon].append('nan')
            # frequency
            if codon in frequency_dict.keys():
                self.codon_anno[codon].append(str(round(frequency_dict[codon] * 1000 / total_aa, 3)))
            else:
                self.codon_anno[codon].append('nan')
            # AA count
            if codon in rscu_dict.keys():
                self.codon_anno[codon].append(str(frequency_dict[codon]))
            else:
                self.codon_anno[codon].append('nan')

        with open(self.output_aa, 'w') as out_aa:
            out_aa.writelines('\t'.join(['Codon', 'AA', 'Abbr.', 'synonym', 'CAI', 'RSCU', 'frequency', 'Count']) + '\n')
            for codon, usage in self.codon_anno.items():
                out_aa.writelines('\t'.join([codon, '\t'.join(usage)]) + '\n')

    def protein_analysis(self):
        with open(self.output_seq, 'w') as out_seq:
            out_seq.writelines('\t'.join(['ID', 'Seq', 'Length', 'Gravy', 'Flexibility', 'Instability',
                                          'Isoelectric_Point', 'Helix', 'Turn', 'Sheet']) + '\n')
            fa_seq = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))
            for now_idx, now_seq in fa_seq.items():
                now_aa = str(now_seq.seq.translate())
                now_len = str(len(now_aa))
                if '*' in now_aa:
                    sys.stdout.writelines(
                        'Skip sequence contain stop codon: {id}\t{aa}'.format(id=now_idx, aa=now_aa) + '\n')
                    continue
                prpt = SeqProperties(now_aa)
                out_seq.writelines('\t'.join([now_idx, now_aa, now_len, prpt.gravy, prpt.flexibility, prpt.instability,
                                              prpt.isoelectric_point, prpt.helix, prpt.turn, prpt.sheet]) + '\n')
