#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_end.py


import os
import sys
from collections import OrderedDict
from itertools import islice
from Bio import SeqIO
import argparse
import pandas as pd
import pysam


class Mrna(object):

    def __init__(self, record):
        self.chromosome = record[0]
        self.gene_id = record[1]
        self.gene_name = record[2]
        self.transcript_id = record[3]
        self.start = record[4]
        self.end = record[5]
        self.utr5_length = int(record[6])
        self.cds_length = int(record[7])
        self.utr3_length = int(record[8])
        self.strand = record[9]
        self.rep_transcript = record[10]
        self.modified = record[11]
        self.bam = []
        self.rpf = []
        self.seq = None


class Ribo(object):

    def __init__(self, args):
        self.mrna_file = args.transcript
        self.mrna_seq = args.sequence
        self.mrna_dict = OrderedDict()
        self.longest = str(args.longest)

        self.bam = args.bam
        self.sample_format = os.path.splitext(args.bam)[1]
        self.pysam_input = None
        self.pysam_output = None

        # self.psite = args.psite
        self.min_length = args.min
        self.max_length = args.max
        self.offset = {}
        self.profile = None

        self.output_prefix = args.output

    def read_transcript(self):

        record = SeqIO.parse(self.mrna_seq, "fasta")
        mrna_sequence = OrderedDict()
        for line in record:
            # sys.stdout.writelines("import gene:  {gene}\n".format(gene=line.id))
            mrna_sequence[line.id] = line.seq

        # mrna_sequence = pysam.FastaFile(self.mrna_seq)

        with open(self.mrna_file, 'r') as trans_file_in:
            if self.longest:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    # if int(record[7]) % 3 != 0:
                    #     sys.stdout.write("Delete {gene}. CDS length doesn't fit 3nt periodicity. \n".format(gene=record[3]))
                    #     continue

                    if record[10] == "True":
                        now_mrna = Mrna(record)
                        # sys.stdout.writelines(now_mrna.transcript_id + '\n')
                        now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                        try:
                            now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                            now_mrna.rpf = [0] * len(now_mrna.seq)
                            self.mrna_dict[now_mrna.transcript_id] = now_mrna
                        except KeyError:
                            continue
                    else:
                        continue
            else:
                for line in islice(trans_file_in, 1, None):
                    record = line.strip().split('\t')
                    now_mrna = Mrna(record)
                    # sys.stdout.writelines(now_mrna.transcript_id + '\n')
                    now_mrna.length = now_mrna.utr5_length + now_mrna.utr3_length + now_mrna.cds_length
                    try:
                        now_mrna.seq = mrna_sequence[now_mrna.transcript_id]
                        now_mrna.rpf = [0] * len(now_mrna.seq)
                        self.mrna_dict[now_mrna.transcript_id] = now_mrna
                    except KeyError:
                        continue

    @staticmethod
    def get_tes_site(utr5, cds, utr3):
        start_site = (utr5 + cds) - 100
        end_site = (utr5 + cds) + 100
        if start_site < 0:
            start_site = 0
        if end_site > utr5 + cds + utr3:
            end_site = utr5 + cds + utr3

        return start_site, end_site

    @staticmethod
    def get_tis_site(utr5, cds, utr3):
        start_site = utr5 - 100
        end_site = utr5 + 100
        if start_site < 0:
            start_site = 0
        if end_site > utr5 + cds:
            end_site = utr5 + cds

        return start_site, end_site

    def rpf_end_sum(self):
        rpf_list = []
        self.pysam_input = pysam.AlignmentFile(self.bam, "rb")
        self.pysam_output = pysam.AlignmentFile(self.output_prefix + "_end.bam", 'wb', template=self.pysam_input)
        for name, isoform in self.mrna_dict.items():
            start_site, end_site = self.get_tis_site(isoform.utr5_length, isoform.cds_length, isoform.utr3_length)
            for read in self.pysam_input.fetch(name, start_site, end_site):
                self.pysam_output.write(read)

            start_site, end_site = self.get_tes_site(isoform.utr5_length, isoform.cds_length, isoform.utr3_length)
            for read in self.pysam_input.fetch(name, start_site, end_site):
                self.pysam_output.write(read)

                # if self.min_length <= read_length <= self.max_length:
                #     try:
                #         p_site = map_start + self.offset[read_length][0]
                #         self.mrna_dict[line.reference_name].rpf[p_site] += 1
                #         p_site = map_start + self.offset[read_length][1]
                #         self.mrna_dict[line.reference_name].rpf[p_site] += 1
                #     except KeyError:
                #         # sys.stdout.writelines("skip reads {reads}!".format(reads=line))
                #         continue
                #     except IndexError:
                #         continue
                # else:
                #     pass


def rpf_end_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to parsing the ribo-seq.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument(
        '-t', dest="transcript", required=True, type=str,
        help="the name of input transcript file in TXT format."
    )
    input_group.add_argument(
        '-s', dest="sequence", required=True, type=str,
        help="the name of input transcript sequence file in FA format."
    )
    input_group.add_argument(
        '-b', dest="bam", required=True, type=str,
        help="the name of mapping file in BAM format."
    )
    input_group.add_argument(
        '-o', dest="output", required=True, type=str,
        help="the prefix of output file. (prefix + _rpf.txt)"
    )

    # arguments for the ribo-seq parsing
    parser.add_argument(
        '-l', dest="longest", action='store_true', required=False, default=False,
        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
             "Recommended : True"
    )
    parser.add_argument(
        '-m', dest="min", required=False, type=int, default=27,
        help="the minimum reads length to keep (default: %(default)s nt)."
    )
    parser.add_argument(
        '-M', dest="max", required=False, type=int, default=33,
        help="the maximum reads length to keep (default: %(default)s nt)."
    )

    args = parser.parse_args()
    return args


def main():

    sys.stdout.writelines('\nDetect the p-site offset.\n')
    sys.stdout.writelines('Step1: Checking the input Arguments.\n')
    args = rpf_end_args_parser()
    ribo_attr = Ribo(args)

    sys.stdout.writelines('\nStep2: Import the transcripts annotation.\n')
    ribo_attr.read_transcript()

    sys.stdout.writelines('\nStep3: Import the bam file.\n')
    ribo_attr.rpf_end_sum()

    sys.stdout.writelines('\nAll done.\n')


if __name__ == '__main__':
    main()

