#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : make_ensb_ref.py


import sys

import numpy as np
from Bio import SeqIO

from utils.ribo import ArgsParser
from utils.ribo.Ensembl_Ref import *


def readline(record):
    chrom, source, feature, start, end, score, strand, phase, attr = record.split('\t')

    attr_dict = OrderedDict()
    for num, mess in enumerate(attr.strip(';').split(';')):
        idx, values = mess.strip().split(' "')
        attr_dict[idx] = values.strip('"')

    section = {"chrom": chrom, "source": source, "feature": feature,
               "start": int(start), "end": int(end), "strand": strand, "attr": attr,
               "attr_dict": attr_dict, "mess": record}

    return section


def detect_cds_periodicity(gtf_filename):
    # CDS frame in gtf file could be fit the 3nt periodicity, commonly, the frame of CDS are Closed-interval
    # but some times the frame is Right-open-interval
    transcripts_dict = OrderedDict()
    now_mrna = 0
    with open(gtf_filename, 'r') as gtf_in:
        for line in gtf_in:
            record = line.strip()
            if not record or line.startswith('#'):
                continue

            section = readline(record)
            if section["attr_dict"]["gene_biotype"] == "protein_coding":
                if section["feature"] == "transcript":
                    transcripts_dict[section["attr_dict"]["transcript_id"]] = np.array([0, 0, 0])
                    now_mrna += 1
                    if now_mrna >= 2000:
                        sys.stdout.writelines("The first 1000 genes was used to detect the format of CDS position.\n")
                        break
                elif section["feature"] == "CDS":
                    cds_frame0 = section["end"] - section["start"]
                    cds_frame1 = section["end"] - section["start"] + 1
                    cds_frame2 = section["end"] - section["start"] - 1
                    transcripts_dict[section["attr_dict"]["transcript_id"]] += [cds_frame0, cds_frame1, cds_frame2]

    cds_shift = 1
    cds_type = ["Right-open-interval", "Closed-interval", "Open-interval"]
    for frame in [0, 1, -1]:
        mrna_length_array = np.array(list(map(lambda length: length[frame], transcripts_dict.values())))

        cds_frame = mrna_length_array % 3
        out_frame_num = np.count_nonzero(cds_frame)
        if out_frame_num == 0:
            # print(transcripts_dict)
            sys.stdout.writelines("The frame type of CDS is {now_type}.\n".format(now_type=cds_type[frame]))
            cds_shift = frame
            break
        elif frame == -1:
            sys.stdout.writelines("Some CDS in GTF file does not fit to 3 nt periodicity.\n")
            cds_shift = 1
        else:
            continue

    return cds_shift


def read_gtf(gtf_filename):
    title_list = []
    genes_dict = OrderedDict()
    transcripts_dict = OrderedDict()
    now_row = 0

    cds_shift = detect_cds_periodicity(gtf_filename)

    with open(gtf_filename, 'r') as gtf_in:
        for line in gtf_in:
            now_row += 1
            if now_row % 10000 == 0:
                sys.stdout.writelines("Rows:  {number}\n".format(number=now_row))

            record = line.strip()

            # skip the '#' lines and blank lines
            if not record:
                continue
            if line.startswith('#'):
                title_list.append(record)
                continue

            section = readline(record)

            # merge the gene lines
            if section["feature"] == "gene":
                genes_dict[section["attr_dict"]["gene_id"]] = Gene(section)

            # merge the mRNA
            elif section["feature"] == "transcript":
                now_rna = Transcripts(section)
                transcripts_dict[now_rna.transcript_id] = now_rna

            # merge the exon, cds, start_codon, stop_codon
            elif section["feature"] in ["exon", "CDS", "start_codon", "stop_codon"]:
                trans_id = section["attr_dict"]["transcript_id"]
                transcripts_dict[trans_id].add_feature(section, cds_shift)

            # skip the other genes
            else:
                # sys.stdout.write("Skip: {row} ".format(row=record))
                continue

    sys.stdout.writelines("Rows: {number}\n".format(number=now_row))

    return title_list, genes_dict, transcripts_dict


def gene_tree(utr_len, chroms_dict, transcripts_dict, genes_dict):
    for trans_id, trans_info in transcripts_dict.items():

        if trans_info.gene_id in genes_dict:
            if trans_info.transcript_biotype == "protein_coding":
                try:
                    trans_info.add_utr(utr_len, chroms_dict)
                except IndexError:
                    sys.stdout.write("IndexError: {gene}\n".format(gene=trans_info.gene_id))
                    continue

            genes_dict[trans_info.gene_id].add_transcript(trans_info)

        else:
            raise KeyError("Error: {trans} not found in {gene}!".format(trans=trans_id, gene=trans_info.gene_id))

    return genes_dict


def format_results(gene_mess):
    mrna_gtf = []
    gene_gtf = []
    mrna_txt = []
    mrna_region = []

    gene_gtf.append('\t'.join([gene_mess.chrom, gene_mess.source, gene_mess.feature, str(gene_mess.start),
                               str(gene_mess.end), '.', gene_mess.strand, '.', gene_mess.attr]))

    if gene_mess.gene_type == "protein_coding":
        mrna_gtf.append('\t'.join([gene_mess.chrom, gene_mess.source, gene_mess.feature, str(gene_mess.start),
                                   str(gene_mess.end), '.', gene_mess.strand, '.', gene_mess.attr]))

    for trans_ids, trans_info in gene_mess.transcript.items():
        if trans_info.transcript_biotype == "protein_coding":
            # Determine whether the CDS length is an integer multiple of 3.
            if trans_info.cds_length % 3 != 0:
                sys.stdout.write("Warning! {gene} CDS length doesn't fit the 3nt periodicity. \n".format(gene=trans_ids))
                continue

            if gene_mess.rep_transcript == trans_info.transcript_id:
                rep_transcript = True
            else:
                rep_transcript = False

            mrna_txt.append('\t'.join([trans_info.chrom, gene_mess.gene_id, trans_info.gene_name, trans_info.transcript_id,
                                       str(trans_info.start), str(trans_info.end), str(trans_info.utr5), str(trans_info.cds_length),
                                       str(trans_info.utr3), trans_info.strand, str(rep_transcript), str(trans_info.modified)]))
            mrna_region.append([trans_ids, trans_info.chrom, trans_info.exons, trans_info.strand])

            if not trans_info.modified:
                mrna_gtf.extend(trans_info.mess)
                gene_gtf.extend(trans_info.mess)

            else:
                mrna_gtf.append('\t'.join([trans_info.chrom, trans_info.source, trans_info.feature, str(trans_info.start),
                                           str(trans_info.end), '.', trans_info.strand, '.', trans_info.attr]))
                gene_gtf.append('\t'.join([trans_info.chrom, trans_info.source, trans_info.feature, str(trans_info.start),
                                           str(trans_info.end), '.', trans_info.strand, '.', trans_info.attr]))

                for exon in trans_info.exon_feature:
                    mrna_gtf.append('\t'.join([str(i) for i in exon]))
                    gene_gtf.append('\t'.join([str(i) for i in exon]))

                for cds in trans_info.cds_feature:
                    mrna_gtf.append('\t'.join([str(i) for i in cds]))
                    gene_gtf.append('\t'.join([str(i) for i in cds]))
        else:
            gene_gtf.extend(trans_info.mess)

    return mrna_gtf, gene_gtf, mrna_txt, mrna_region


def filter_genes(genes_dict, coding):
    filtered_gtf = []
    filtered_txt = []
    filtered_region = []

    if coding:
        for gene_name, gene_mess in genes_dict.items():
            mrna_gtf, gene_gtf, mrna_txt, mrna_region = format_results(gene_mess)
            filtered_gtf.extend(mrna_gtf)
            filtered_txt.extend(mrna_txt)
            filtered_region.extend(mrna_region)

    elif not coding:
        for gene_name, gene_mess in genes_dict.items():
            if gene_mess.gene_type == 'protein_coding':
                mrna_gtf, gene_gtf, mrna_txt, mrna_region = format_results(gene_mess)
                filtered_gtf.extend(gene_gtf)
                filtered_txt.extend(mrna_txt)
                filtered_region.extend(mrna_region)
            else:
                mrna_gtf, gene_gtf, mrna_txt, mrna_region = format_results(gene_mess)
                filtered_gtf.extend(gene_gtf)

    return filtered_gtf, filtered_txt, filtered_region


def read_genome(genome):
    chroms_dict = OrderedDict()

    record = SeqIO.parse(genome, "fasta")
    for line in record:
        sys.stdout.writelines("import chromosome: {chrom}\n".format(chrom=line.id))
        chroms_dict[line.id] = Chrom(line)

    return chroms_dict


def get_seq(chroms_dict, mrna_region):
    mrna_seq = OrderedDict()

    for transcript in mrna_region:
        transcript_seq = ''

        if transcript[-1] == "-":
            for exon in reversed(transcript[2]):
                exon_start, exon_end = exon[0], exon[1]
                transcript_seq += chroms_dict[transcript[1]].seq[exon_start - 1: exon_end]
            transcript_seq = transcript_seq.reverse_complement()

        else:
            for exon in transcript[2]:
                exon_start, exon_end = exon[0], exon[1]
                transcript_seq += chroms_dict[transcript[1]].seq[exon_start - 1: exon_end]

        mrna_seq[transcript[0]] = transcript_seq

    return mrna_seq


def output_results(output_prefix, title_list, filtered_gtf, mrna_txt, mrna_seq):
    gtf_out_file = output_prefix + '.norm.gtf'
    with open(gtf_out_file, 'w') as gtf_out:
        for line in title_list:
            gtf_out.writelines(''.join(line) + '\n')

        for line in filtered_gtf:
            gtf_out.writelines(''.join(line) + '\n')

    txt_out_file = output_prefix + '.norm.txt'
    with open(txt_out_file, 'w') as txt_out:
        txt_out.writelines('\t'.join(["chromosome", "gene_id", "gene_name", "transcript_id", "start", "end", "utr5_length",
                                      "cds_length", "utr3_length", "strand", "rep_transcript", "modified"]) + '\n')
        for line in mrna_txt:
            txt_out.writelines(''.join(line) + '\n')

    seq_out_file = output_prefix + '.norm.fa'
    with open(seq_out_file, 'w') as seq_out:
        for mrna, sequence in mrna_seq.items():
            seq_out.writelines('\n'.join([">" + mrna, str(sequence)]) + '\n')


def main():
    ArgsParser.now_time()
    sys.stdout.writelines('\nMake the gene annotation files.\n')
    sys.stdout.writelines('Step1: Checking the input Arguments.\n')
    args = ArgsParser.gtf_args_parser()

    sys.stdout.writelines('\nStep2: Import the gtf file.\n')
    title_list, genes_dict, transcripts_dict = read_gtf(args.transcript)

    sys.stdout.writelines('\nStep3: Import the genome file.\n')
    chroms_dict = read_genome(args.sequence)

    sys.stdout.writelines('\nStep4: Make the gene tree.\n')
    genes_dict = gene_tree(args.utr, chroms_dict, transcripts_dict, genes_dict)

    sys.stdout.writelines('\nStep5: Screening genes.\n')
    filtered_gtf, mrna_txt, mrna_region = filter_genes(genes_dict, args.coding)

    sys.stdout.writelines('\nStep6: Retrieve the mRNA sequence from genome.\n')
    mrna_seq = get_seq(chroms_dict, mrna_region)

    sys.stdout.writelines('\nStep7: Output the results.\n')
    output_results(args.output, title_list, filtered_gtf, mrna_txt, mrna_seq)

    sys.stdout.writelines('\nALL DONE!\n\n')
    ArgsParser.now_time()


if __name__ == "__main__":
    main()
