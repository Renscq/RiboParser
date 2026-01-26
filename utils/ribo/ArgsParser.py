#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : ArgsParser.py


import argparse
import os
import sys
import textwrap
import time
import warnings
warnings.filterwarnings('ignore')


def now_time():
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)


def args_print(args):
    args_dict = vars(args)

    for k, v in args_dict.items():
        print("{:<12}:  {:<}".format(k, str(v)), flush=True)

    sys.stdout.flush()


def file_check(*files):
    for my_file in files:
        if os.path.exists(my_file):
            continue
        else:
            print("\nFile {file_name} is not exists!\n".format(file_name=my_file), flush=True)
            sys.exit()


def gtf_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to build the references.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the input file name of genome sequence")
    input_group.add_argument('-t', dest="transcript", required=True, type=str, help="the input file name of gtf file")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _norm.gtf)")

    # arguments for the modification
    parser.add_argument('-u', dest="utr", required=False, type=int, default=99,
                        help="add the pseudo UTR to the leaderless transcripts (default: %(default)s nt).")
    parser.add_argument('-c', dest="coding", required=False, action="store_true", default=False,
                        help="only retain the protein coding transcripts (default: %(default)s).")

    args = parser.parse_args()
    file_check(args.sequence, args.transcript)
    args_print(args)

    return args


def make_ribo_ref():
    parser = argparse.ArgumentParser(description="This script is used to build the references.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-g', dest="genome", required=True, type=str,
                             help="the input file name of genome sequence")
    input_group.add_argument('-t', dest="gtf", required=True, type=str, help="the input file name of gtf file")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _norm.gtf)")

    # arguments for the modification
    parser.add_argument('-u', dest="utr", required=False, type=int, default=0,
                        help="add the pseudo UTR to the leaderless transcripts (default: %(default)s nt).")
    parser.add_argument('-c', dest="coding", required=False, action="store_true", default=False,
                        help="only retain the protein coding transcripts (default: %(default)s).")
    parser.add_argument('-l', dest="longest", required=False, action="store_true", default=False,
                        help="only retain the longest protein coding transcripts (default: %(default)s).")
    parser.add_argument('-w', dest="whole", required=False, action="store_true", default=False,
                        help="output whole message (default: %(default)s).")
    
    args = parser.parse_args()
    file_check(args.genome, args.gtf)
    args_print(args)

    return args


############################################################
def rpf_bam_check_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to summary the BAM condition.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the input file name of gene annotation")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the input file name of bam")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the modification
    parser.add_argument('--thread', dest="thread", required=False, type=int, default=1,
                        help='''the number of threads (default: %(default)s). 
                        Suitable for large bam files > 1G.
                        It will take a lot of memory.''')
    parser.add_argument('-g', dest="tag", choices=[0, 1], type=int, required=False, default=0,
                        help='''filter the number of reads mapped loci (default: %(default)s).
                        [0]: all reads will be used; [1]: reads with unique mapped loci will be used
                        ''')
    parser.add_argument('-a', dest="align", choices=['star', 'hisat2', 'bowtie2'], type=str, required=False,
                        default='star', help='''filter the number of reads mapped loci (default: %(default)s).
                        ''')
    parser.add_argument('-r', dest="reverse", action='store_true', required=False, default=False,
                        help="reads aligned to negative strand will also be counted. (default: %(default)s).")
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only keep the longest transcripts (default: %(default)s).")
    parser.add_argument('-s', dest="saturation", action='store_true', required=False, default=False,
                        help='''whether to calculate RPF saturation. (default: %(default)s).
                        This step will take a lot of time and memory.
                        ''')
    args = parser.parse_args()
    file_check(args.transcript, args.bam)
    args_print(args)

    return args


def digestion_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to Detect the digestion sites.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript file in TXT format.")
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the name of input transcript sequence file in FA format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the name of mapping file in BAM format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the name of output file. (prefix + _digestion_sites.txt)")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
                             "Recommended : True")
    parser.add_argument('--scale', dest="scale", action='store_true', required=False, default=False,
                        help="scale the motif matrix (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=20,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=100,
                        help="the maximum reads length to keep (default: %(default)s nt).")

    args = parser.parse_args()
    file_check(args.transcript, args.sequence, args.bam)
    args_print(args)

    return args


def offset_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to detect the P-site offset.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript filein TXT format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the name of mapping file in BAM format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _offset.txt)")

    # arguments for the offset detection
    # parser.add_argument('--mode', dest="mode", required=False, type=str, default='SSCBM', choices=['SSCBM', 'RSBM'],
    #                     help="specify the mode of offset detect [SSCBM, RSBM]. (default: %(default)s).")
    parser.add_argument('-a', dest="align", required=False, type=str, default='both', choices=['both', 'tis', 'tts'],
                        help="specify the alignment of reads for offset detect [both, tis, tts]. (default: %(default)s).")
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=27,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=33,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-p', dest="exp_peak", required=False, type=int, default=30,
                        help="Expected RPFs length fitted to ribosome structure [~30 nt] (default: %(default)s nt).")
    parser.add_argument('-s', dest="shift", required=False, type=int, default=2,
                        help="psite shift for different RPFs length. (default: %(default)s nt).")
    parser.add_argument('--silence', dest="silence", required=False, action='store_true', default=True,
                        help="discard the warning information. (default: %(default)s).")
    parser.add_argument('-d', dest="detail", action='store_true', required=False, default=False,
                        help="output the details of offset (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.transcript, args.bam)
    args_print(args)

    return args


def offset_rsbm_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to detect the P-site offset.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript filein TXT format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the name of mapping file in BAM format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _offset.txt)")

    # arguments for the offset detection
    # parser.add_argument('--mode', dest="mode", required=False, type=str, default='SSCBM',
    #                     help="specify the mode of offset detect [RSBM, SSCBM]. (default: %(default)s).")
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=27,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=33,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-p', dest="exp_peak", required=False, type=int, default=29,
                        help="RPFs peak length [~30 nt] (default: %(default)s nt).")
    parser.add_argument('-e', dest="exp_offset", required=False, type=int, default=11,
                        help="expected offset length (default: %(default)s).")
    parser.add_argument('-s', dest="shift", required=False, type=int, default=2,
                        help="psite shift for different RPFs length."
                             "Empirical value: 2nt for Eukaryotes, 1nt for prokaryotes. (default: %(default)s nt).")
    parser.add_argument('--silence', dest="silence", required=False, action='store_true', default=True,
                        help="discard the warning information. (default: %(default)s).")
    parser.add_argument('-d', dest="detail", action='store_true', required=False, default=False,
                        help="output the details of offset (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.transcript, args.bam)
    args_print(args)

    return args

def rna_offset_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to create the offset table of RNA-seq.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _offset.txt)")

    # arguments for the offset table creation
    parser.add_argument('-m', dest="min", required=False, type=int, default=25,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=151,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-e', dest="exp_offset", required=False, type=int, default=12,
                        help="Expected offset (default: %(default)s nt).")

    args = parser.parse_args()
    args_print(args)

    return args

def rpf_bam2bw_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to convert the bam to bedgraph.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-b', dest="bam", required=True, type=str,
                             help="name of input genome alignment file (bam/sam format).")
    input_group.add_argument('-p', dest="psite", required=True, type=str,
                            help="p-site offset file in TXT format. (default: %(default)s).")

    # arguments for the RPFs calculation
    parser.add_argument('-t', dest="times", required=False, type=int, default=3,
                        help="set the reads multiple aligned times. (default: %(default)s).")
    parser.add_argument('--second', dest="secondary", required=False, action='store_true', default=False,
                        help="discard the reads aligned secondary loci. (default: %(default)s).")
    parser.add_argument('--supply', dest="supplementary", required=False, action='store_true', default=False,
                        help="discard the supplementary reads (default: %(default)s ).")    
    parser.add_argument('-f', dest="format", required=False, choices=['bedgraph', 'wig'], type=str, default='bedgraph',
                        help="output the bed file. (default: %(default)s).")
    parser.add_argument('-n', dest="norm", required=False, action='store_true', default=False,
                        help="normalise the RPFs to RPM. (default: %(default)s).")
    parser.add_argument('-m', dest="merge", required=False, action='store_true', default=False,
                        help="merge minus and plus strand to one file. (default: %(default)s).")
    parser.add_argument('-o', dest="output", required=False, type=str,
                        help="output the bed file. (default: prefix + .bedgraph).")

    args = parser.parse_args()
    file_check(args.bam, args.psite)
    args_print(args)

    return args


def bam_filter_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to filter the reads length from bam file.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-i', dest="ibam", required=True, type=str,
                             help="the name of input bam file.")
    input_group.add_argument('-o', dest="obam", required=True, type=str,
                             help="the name of output bam file.")
    
    # arguments for the RPFs calculation
    parser.add_argument('-l', dest="length", required=False, type=str, default="27,28,29,30,31,32",
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-u', dest="unique", required=False, action='store_true', default=False,
                        help="discard the secondary reads (default: %(default)s ).")
    parser.add_argument('-s', dest="supplementary", required=False, action='store_true', default=False,
                        help="discard the supplementary reads (default: %(default)s ).")
    parser.add_argument('-q', dest="quality", required=False, type=int, default=28,
                        help="discard the low mapping quality reads (default: %(default)s ).")

    args = parser.parse_args()
    file_check(args.ibam)
    args_print(args)

    return args


def ribo_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to convert Ribo-seq bam to p-site density.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript file in TXT format.")
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the name of input transcript sequence file in FA format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, 
                             help="the name of mapping file in BAM format.")
    input_group.add_argument('-p', dest="psite", required=True, type=str,
                             help="the name of p-site offset file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (output = prefix + _rpf.txt)")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
                             "Recommended : True")
    parser.add_argument('-m', dest="min", required=False, type=int, default=27,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=33,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    input_group.add_argument('--period', dest='periodicity', required=False, type=float, default=40,
                             help="the minimum 3nt periodicity to keep. (default: %(default)s).")
    parser.add_argument('--silence', dest="silence", required=False, action='store_true', default=True,
                        help="discard the warning information. (default: %(default)s).")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads (default: %(default)s). It will take a lot of memory.")
    args = parser.parse_args()
    file_check(args.transcript, args.sequence, args.bam, args.psite)
    args_print(args)

    return args


def rna_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to convert RNA-seq bam to p-site density.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript file in TXT format.")
    input_group.add_argument('-s', dest="sequence", required=True, type=str,
                             help="the name of input transcript sequence file in FA format.")
    input_group.add_argument('-b', dest="bam", required=True, type=str, help="the name of mapping file in BAM format.")
    input_group.add_argument('-p', dest="psite", required=True, type=str,
                             help="the name of p-site offset file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _rna.txt)")

    # arguments for the RNA-seq parsing
    parser.add_argument('-l', dest="longest", action='store_true', required=False, default=False,
                        help="only retain the transcript with longest CDS of each gene (default: %(default)s)."
                             " Recommended : True")
    parser.add_argument('-m', dest="min", required=False, type=int, default=25,
                        help="the minimum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-M', dest="max", required=False, type=int, default=150,
                        help="the maximum reads length to keep (default: %(default)s nt).")
    parser.add_argument('-r', dest="rolling", action='store_true', required=False, default=False,
                        help="density is calculated once per ribosome width. (default: %(default)s)."
                             "only suitable for long RNA-seq reads.")
    parser.add_argument('--pe', dest="pair_end", action='store_true', required=False, default=False,
                        help="reads aligned to negative strand will also be counted. (default: %(default)s)."
                             "only suitable for pair-end RNA-seq reads.")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads (default: %(default)s). It will take a lot of memory.")
    args = parser.parse_args()
    file_check(args.transcript, args.sequence, args.bam, args.psite)
    args_print(args)

    return args


def shuffle_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to shuffle the RPFs data.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs density file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _shuffle.txt)")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="list", required=False, type=str,
                        help="the list of input genes for transcript id.")
    parser.add_argument('-s', dest="seed", required=False, type=int, default=0,
                        help="the random seed for shuffle. (default: %(default)s).")
    parser.add_argument('-i', dest="individual", action='store_true', required=False, default=False,
                        help="shuffle the RPFs data for each samples. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def retrieve_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to retrieve density files.")

    # arguments for the required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs density file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, 
                             help="prefix of output file name (default: filename + '_retrieve.txt'.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="list", required=False, type=str,
                        help="the list of input genes for transcript id.")
    parser.add_argument('-m', dest="min", required=False, type=int, default=0,
                        help="retain transcript with more than minimum RPFs (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM (default: %(default)s).")
    parser.add_argument('-f', dest="format", action='store_true', required=False, default=False,
                        help="melt three column data of each sample to one column (default: %(default)s).")
    parser.add_argument('-s', dest="split", action='store_true', required=False, default=False,
                        help="split gene rpf to each TXT file (default: %(default)s).")

    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args

def frame_shift_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the frame shift plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument(
        '-o', dest="output", required=True, type=str,
        help="the prefix of output file."
    )

    # arguments for the ribo-seq parsing
    input_group.add_argument('-t', dest="transcript", required=False, type=str,
                             help="the name of input transcript filein TXT format.")
    parser.add_argument('-p', dest="period", required=False, type=int, default=45,
                        help="the minimum in-frame value for frame shifting screen, range [0 - 1]. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="the number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="the number of codons before TTS will be discarded.. (default: %(default)s AA).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_merge_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to merge the density file.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-l', dest="list", required=True, type=str, help="the sample list in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _rpf_merged.txt)")

    args = parser.parse_args()
    args_print(args)

    return args


def periodicity_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the periodicity plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument(
        '-o', dest="output", required=True, type=str,
        help="the prefix of output file."
    )

    # arguments for the ribo-seq parsing
    input_group.add_argument('-t', dest="transcript", required=False, type=str,
                             help="the name of input transcript filein TXT format.")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="the number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="the number of codons before TTS will be discarded.. (default: %(default)s AA).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def metaplot_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the meta plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript filein TXT format.")
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix name of output file.")
    
    # arguments for the ribo-seq parsing
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="delete transcript with less than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--utr5', dest="utr5", type=int, required=False, default=20,
                        help="the codon number in 5-utr region (default: %(default)s AA).")
    parser.add_argument('--cds', dest="cds", type=int, required=False, default=50,
                        help="the codon number in cds region (default: %(default)s AA).")
    parser.add_argument('--utr3', dest="utr3", type=int, required=False, default=20,
                        help="the codon number in 3-utr region (default: %(default)s AA).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--mode', dest="mode", required=False, choices=["line", "bar"], type=str, default='bar',
                        help="specify the mode of metaplot. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def coverage_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the coverage meta plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript filein TXT format.")
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, help="the prefix of output file.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-b', dest="bin", required=False, type=str, default='30,100,30',
                        help="adjust the transcript to specified bins. 30 for 5'-UTR"
                             "and 3'-UTR, 100 for CDS. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads. (default: %(default)s).")
    parser.add_argument('--outlier', dest="outlier", action='store_true', required=False, default=False,
                        help="filter the outliers (default: %(default)s).")
    parser.add_argument('--set', dest="set", choices=['intersect', 'union'], required=False, type=str, default='union',
                        help="filter the gene list with 5-UTR / CDS / 3-UTR. (default: %(default)s).")
    parser.add_argument('--heat', dest="heatmap", action='store_true', required=False, default=False,
                        help="draw the coverage heatmap of whole gene. (default: %(default)s).")
    parser.add_argument('--bar', dest="barplot", action='store_true', required=False, default=False,
                        help="draw the coverage barplot of whole gene. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def percentage_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the plot of rpf coverage percentage.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-t', dest="transcript", required=True, type=str,
                             help="the name of input transcript filein TXT format.")
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, help="the prefix of output file.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_corr_args_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to draw the correlation of rpf density.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _rpf_merged.txt)")

    args = parser.parse_args()
    args_print(args)

    return args


def meta_codon_plot_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the meta codon plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-c', dest="codon", required=False, type=str, help="the codon list in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-a', dest="around", required=False, type=int, default=20,
                        help="retrieve length of codon upstream and downstream. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="The number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="The number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('-u', dest="unique", action='store_true', required=False, default=False,
                        help="delete the cross repetition codon in different window. (default: %(default)s).")
    parser.add_argument('-s', dest="scale", action='store_true', required=False, default=False,
                        help="scale the window density with gene density. (default: %(default)s).")
    parser.add_argument('--smooth', dest="smooth", required=False, default=None, type=str,
                        help="smooth the window density [eg, 3,1]. (default: %(default)s).")
    parser.add_argument('--thread', dest="thread", type=int, required=False, default=1,
                        help="the number of threads (default: %(default)s).")
    parser.add_argument('--fig', dest="fig", required=False, action='store_true', default=False,
                        help="output the figure. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def codon_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the codon usage plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=30,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=15,
                        help="The number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=5,
                        help="The number of codons before TTS will be discarded.. (default: %(default)s AA).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_gene_plot_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the gene rpf plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    gene_name = parser.add_mutually_exclusive_group()
    gene_name.add_argument('-l', dest="list", required=False, type=str, help="the gene name list in TXT format.")
    gene_name.add_argument('-g', dest="gene", required=False, type=str, help="the gene name.")

    # arguments for the RPFs calculation
    parser.add_argument('--log', dest="log", required=False, choices=['2', '10'], type=str, default='Not',
                        help="set the y-axis to log scaling. (default: %(default)s).")

    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_pausing_args_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to calculate the relative pausing score.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the pausing calculation
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for pausing calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for pausing calculation. (default: %(default)s).")
    parser.add_argument('-b', dest="background", required=False, type=int, default=2,
                        help="set the codon number before and after p-site as the background. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=10,
                        help="The number of codons after TIS will be discarded. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=5,
                        help="The number of codons before TTS will be discarded. (default: %(default)s AA).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the pausing score. (default: %(default)s).")
    parser.add_argument('--stop', dest="stop",  action='store_true', required=False, default=False,
                        help="remove the stop codon. (default: %(default)s).")    
    parser.add_argument('--fig', dest="figure", choices=['none', 'png', 'pdf'], required=False, default='none',
                        help="draw the rpf pausing score of each gene (it will takes a lot of time). (default: %(default)s).")
    parser.add_argument('--all', dest="all", action='store_true', required=False, default=False,
                        help="output all pausing score of each gene. (default: %(default)s).")

    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_occupancy_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the codon occupancy plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the RPFs calculation
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=30,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=15,
                        help="The number of codons after TIS will be discarded. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=5,
                        help="The number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the occupancy. (default: %(default)s).")
    parser.add_argument('--stop', dest="stop",  action='store_true', required=False, default=False,
                        help="rmove the stop codon. (default: %(default)s).")
    parser.add_argument('--all', dest="all", action='store_true', required=False, default=False,
                        help="output all RPFs density. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_cdt_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the codon decoding time plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('--rpf', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('--rna', dest="rna", required=True, type=str,
                             help="the name of input reads file in TXT format.")
    input_group.add_argument('-l', dest="list", required=True, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the RPFs calculation
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for codon decoding time calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for codon decoding time calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=30,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    # parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
    #                     help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=15,
                        help="The number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=5,
                        help="The number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the codon decoding time. (default: %(default)s).")
    parser.add_argument('--stop', dest="stop",  action='store_true', required=False, default=False,
                        help="rmove the stop codon. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_cst_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the codon decoding time plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('--rpf', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('--rna', dest="rna", required=True, type=str,
                             help="the name of input reads file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the RPFs calculation
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for codon decoding time calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for codon decoding time calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=30,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-t', dest="times", required=False, type=int, default=10,
                        help="Specify the number of iteration times required for computation. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="The number of codons after TIS will be discarded. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="The number of codons before TTS will be discarded. (default: %(default)s AA).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the codon selection time. (default: %(default)s).")
    parser.add_argument('--stop', dest="stop",  action='store_true', required=False, default=False,
                        help="rmove the stop codon. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    file_check(args.rna)
    args_print(args)

    return args


def rpf_cov_args_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to calculate CoV in the CDS region.",
                                     epilog=textwrap.dedent('''\
        The group file needs to contain at least two column:
        +----------+---------+
        | name     | group  |
        +==========+=========+
        | wt1      | wt      |
        | wt2      | wt      |
        | treat1   | treat   |
        | treat2   | treat   |
        | ko1      | ko      |
        | ko2      | ko      |
        +----------+---------+
    '''))

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-g', dest="group", required=False, type=str, default=None,
                            help="specify the list of sample group. (default: None)")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (prefix + _Cov.txt)")

    # arguments for the CoV calculation
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=5,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=15,
                        help="The number of codons after TIS will be discarded. (default: %(default)s).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=5,
                        help="The number of codons before TES will be discarded. (default: %(default)s).")
    parser.add_argument('--fig', dest="figure", action='store_true', required=False, default=False,
                        help="show the figure. (default: %(default)s).")
    
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)
    return args


def cumulative_cov_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to calculate the cumulative CoV.")

    # arguments for the required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs density file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, 
                             help="prefix of output file name (default: filename + '_cumulative_CoV.txt'.")

    # arguments for the ribo-seq parsing
    parser.add_argument('-l', dest="list", required=False, type=str,
                        help="the list of input genes for transcript id.")
    parser.add_argument('-m', dest="min", required=False, type=int, default=0,
                        help="retain transcript with more than minimum RPFs (default: %(default)s).")
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM (default: %(default)s).")
    
    parser.add_argument('-t', dest="trim", required=False, type=int, default=50,
                        help="trim transcript with specific length (default: %(default)s nt).")
    
    parser.add_argument('-s', dest="split", action='store_true', required=False, default=False,
                        help="split gene rpf to each TXT file (default: %(default)s).")
    parser.add_argument('-z', dest="zero", action='store_true', required=False, default=False,
                        help="set the start site to zero (default: %(default)s).")
    
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_odd_ratio_args_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to calculate the odd ratio of stalling codon.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-l', dest="list", required=False, type=str, default=None,
                             help="the gene name list in TXT format. (default: whole).")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file.")

    # arguments for the odd ratio calculation
    parser.add_argument('--thread', dest="thread", required=False, type=int, default=1,
                        help="set the thread number. (default: %(default)s).")
    
    parser.add_argument('-s', dest="site", choices=['E', 'P', 'A'], required=False, type=str, default='P',
                        help="set the E/P/A-site for odd ratio calculation. (default: %(default)s).")
    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for odd ratio calculation. (default: %(default)s).")
    parser.add_argument('-c', dest="control", required=True, type=str,
                        help="specify the name of control samples, separated by commas.")
    parser.add_argument('-t', dest="treat", required=True, type=str,
                        help="specify the name of treat samples, separated by commas.")
    
    parser.add_argument('-n', dest="normal", action='store_true', required=False, default=False,
                        help="normalize the RPFs count to RPM. (default: %(default)s).")
    parser.add_argument('-z', dest="zero", action='store_true', required=False, default=False,
                        help="fill the empty codon with [0.01 * minimal value]. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="specify the minimum count of RPFs per transcript to keep. (default: %(default)s).")
    parser.add_argument('--fdr', dest="fdr", required=False, type=str, choices=['bhfdr', 'pvalue'], default='bhfdr',
                        help="specify the multipletests. (default: %(default)s).")
    parser.add_argument('-v', dest="value", required=False, type=float, default=0.05,
                        help="specify the p value. (default: %(default)s).")
    
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="specify the number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="specify the number of codons before TTS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--stop', dest="stop", action='store_true', required=False, default=False,
                        help="remove the stop codon from table. (default: %(default)s).")
    parser.add_argument('--scale', dest="scale",  choices=['zscore', 'minmax'], required=False, type=str, default='minmax',
                        help="normalize the codon odd ratio. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def rpf_quant_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to quantify RPFs in the CDS region.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="the prefix of output file. (default: prefix + _rpf_quant.txt)")

    parser.add_argument('-f', dest="frame", choices=['0', '1', '2', 'all'], required=False, type=str, default='all',
                        help="set the reading frame for occupancy calculation. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="The number of codons after TIS will be discarded. (default: %(default)s).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="The number of codons before TES will be discarded. (default: %(default)s).")
    parser.add_argument('--utr5', dest="utr5", action='store_true', required=False, default=False,
                        help="quantification of 5'-utr. (default: %(default)s).")
    parser.add_argument('--utr3', dest="utr3", action='store_true', required=False, default=False,
                        help="quantification of 3'-utr. (default: %(default)s).")

    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args


def deseq2_args_parser():
    parser = argparse.ArgumentParser(
        description="This script is used to perform differential analysis with riboParser output.")
    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs count file in TXT format.")
    parser.add_argument("-d", dest="design", required=True, type=str, help="input whole design file")
    input_group.add_argument("-o", dest="output", required=True, type=str, help="output the DESeq2 results")

    parser.add_argument("-c", dest="control", required=True, type=str, help="input the control group name")
    parser.add_argument("-t", dest="treat", required=True, type=str, help="input the treat group name")
    parser.add_argument("-m", dest="min", required=False, type=int, default=5,
                        help="specify the minimum reads count. (default: %(default)s).")
    parser.add_argument("-l", dest="logfc", required=False, type=float, default=1,
                        help="specify the logfc threshold for volcano. (default: %(default)s).")
    parser.add_argument("-p", dest="padj", required=False, type=float, default=0.05,
                        help="specify the padj threshold for volcano. (default: %(default)s).")
    args = parser.parse_args()
    file_check(args.rpf)
    file_check(args.design)
    args_print(args)

    return args


def calc_te_args_parser():
    parser = argparse.ArgumentParser(description="Calculate translation efficiency.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="ribo", required=True, type=str,
                             help="specify the input Ribo-seq expression table generate by DESeq2.")
    input_group.add_argument('-m', dest="mrna", required=True, type=str, default=None,
                             help="specify the input RNA-seq expression table generate by DESeq2.")
    input_group.add_argument('--rd', dest="ribo_design", required=True, type=str,
                             help="specify the input Ribo-seq expression table generate by DESeq2.")
    input_group.add_argument('--md', dest="mrna_design", required=True, type=str, default=None,
                             help="specify the input RNA-seq expression table generate by DESeq2.")
    input_group.add_argument('-o', dest="output", required=True, type=str, default=None,
                             help="prefix of output filename.")

    # arguments for the RPFs calculation
    parser.add_argument("--min", dest="min", required=False, type=int, default=10,
                        help="specify the minimum reads count. (default: %(default)s).")
    parser.add_argument("-l", dest="logfc", required=False, type=float, default=1,
                        help="specify the logfc threshold for volcano. (default: %(default)s).")
    parser.add_argument("-t", dest="type", choices=['pvalue', 'padj'], required=False, type=str, default='pvalue',
                        help="specify the type of probability value ['pvalue', 'padj']. (default: %(default)s).")
    parser.add_argument("-p", dest="pvalue", required=False, type=float, default=0.05,
                        help="specify the padj threshold for volcano. (default: %(default)s).")
    args = parser.parse_args()

    file_check(args.ribo)
    file_check(args.mrna)
    file_check(args.ribo_design)
    file_check(args.mrna_design)
    args_print(args)

    return args


def ppi_args_parser():
    parser = argparse.ArgumentParser(description="Draw the protein protein interaction network.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-i', dest="input", required=True, type=str, help="specify the input gene list.")
    input_group.add_argument("-x", dest="taxonomy", required=True, type=str,
                             default=None,
                             help="specify the NCBI taxonomy (eg. human: 9606). (default: %(default)s).")
    input_group.add_argument("-p", dest="scopes", required=False,
                             choices=['name', 'entrezgene', 'ensembl', 'symbol', 'uniprot', 'refseq'],
                             default='ensembl',
                             help="specify the input gene type. (default: %(default)s).")
    input_group.add_argument('-o', dest="output", required=True, type=str, default=None,
                             help="prefix of output filename.")

    # arguments for the RPFs calculation
    parser.add_argument('-num', dest="number", required=False, type=str, default=50,
                        help="specify the input number of input gene for network (default: %(default)s).")
    parser.add_argument("-f", dest="flavor", required=False, type=str,
                        choices=['evidence', 'confidence', 'actions'],
                        help="the style of edges in the network: evidence, confidence (default), actions. (default: %(default)s).")
    parser.add_argument("-s", dest="score", required=False, type=int,
                        default=500,
                        help="specify the taxonomy . (default: %(default)s).")
    parser.add_argument("-t", dest="type", required=False, type=str,
                        choices=['functional', 'physical'],
                        default='physical', help="network type (default: %(default)s).")
    parser.add_argument("-n", dest="nodes", required=False, type=int,
                        choices=[0, 1],
                        default=0,
                        help="hides all protein names from the picture (0 or 1) (default: %(default)s).")
    parser.add_argument("-l", dest="labels", required=False, type=int,
                        choices=[0, 1],
                        default=0,
                        help="when provided use submitted names as protein labels in the network image (0 or 1) (default: %(default)s).")
    parser.add_argument("-c", dest="connected", required=False, type=int,
                        choices=[0, 1],
                        default=0,
                        help="hides all proteins that are not connected to any other protein in your network (0 or 1) (default: %(default)s).")
    parser.add_argument("-b", dest="block", required=False, type=int,
                        choices=[0, 1],
                        default=0,
                        help="disables structure pictures inside the bubble (0 or 1) (default: %(default)s).")
    parser.add_argument('-g', dest="graphic", required=False, choices=['highres_image', 'svg'], type=str,
                        default='svg', help="specify the graphic format (default: %(default)s).")
    args = parser.parse_args()

    file_check(args.input)
    args_print(args)

    return args


def serp_peak_args_parser():
    parser = argparse.ArgumentParser(
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=5, width=90),
        description="This script is used to detect the selective translatome.")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    # arguments for the input files
    input_group = parser.add_argument_group('Input files arguments')
    input_group.add_argument("-r", dest="rpf", required=True, type=str, help="input the RPF coverage file")
    input_group.add_argument("-n", dest="norm", required=False, type=str, default=None,
                             help="specify the file contain RPFs number for normalise. "
                                  "(default: The total RPFs count of input samples is used.)")
    input_group.add_argument("--scale", dest="scale", required=False, type=float, default=1e6,
                             help="specify the scale level for reads normalise (default RPM %(default)s )")
    input_group.add_argument("-a", dest="anno", required=False, type=str,
                             help="specify the gene annotation file in txt format.")
    input_group.add_argument("--ck", dest="control", required=True, type=str,
                             help="specify the name of control samples, separated by commas.")
    input_group.add_argument("--ip", dest="ip", required=True, type=str,
                             help="specify the name of immunoprecipitation samples, separated by commas.")

    # arguments for the data filtering
    data_group = parser.add_argument_group('Data filtering arguments')
    data_group.add_argument("-m", dest="min", required=False, type=int, default=50,
                            help='''specify the minimum number of RPFs coverage for genes to be retained
             (default %(default)s)''')
    data_group.add_argument("--corr", dest="corr", required=False, type=float, default=0.3,
                            help="specify the minimum correlation of replicates (default %(default)s)")

    # arguments for the peaks scanning, gaps and ratio filtering
    peak_group = parser.add_argument_group('Peak scanning arguments')
    peak_group.add_argument("-f", dest="fill", required=False, type=int, default=30,
                            help='''Use the mean value of gene RPF to fill the missing values.
        (Commonly used first 30 AA, default %(default)s codon). Because the
         ribo-seq data is too noisy in single codon level and many locations 
         are not covered by RPF, it is necessary to fill the blank position 
         use the mean value of RPF as the background. Three modes are provided: 
         (1): [30] average of the first 30 AA  is used (or other length); 
         (2): [-1] the average RPFs value of current gene is used; 
         (3): [0] the average RPF of total genes is used.''')
    
    peak_group.add_argument("-s", dest="size", required=False, type=int, default=3,
                            help="Specifies the window size, the numbers must be odd. "
                                 "(default %(default)s AA)."
                                 "The sequencing data is usually noisy, and need to be smoothed."
                                 "savgol_filter is used here. ")
    peak_group.add_argument("-k", dest="k", required=False, type=int, default=1,
                            help="Specifies the polyorder, the numbers must be odd. "
                                 "(default %(default)s AA)."
                                 "The sequencing data is usually noisy, and need to be smoothed."
                                 "savgol_filter is used here. ")
    # peak_group.add_argument("--mode", dest="k", required=False, type=str, default="nearest", 
    #                         choices=["mirror", "nearest", "constant", "wrap"],
    #                         help="Specifies the smooth mode options, (default %(default)s AA).")
    
    peak_group.add_argument("-w", dest="width", required=False, type=int, default=5,
                            help="specify the width of binding peaks (default %(default)s AA)."
                                 "By default, the width of each peak is not less than 5 amino acids.")
    peak_group.add_argument("-e", dest="enrich", required=False, type=float, default=2.0,
                            help="specify the enrich threshold of each peak height (default %(default)s)")
    peak_group.add_argument("-c", dest="collision", required=False, type=float, default=1.5,
                            help="specify the enrich threshold of collision region of binding peak (default %(default)s)")
    peak_group.add_argument("-g", dest="gaps", required=False, type=int, default=1,
                            help='''specify the gap width inside each peak. (default %(default)s AA).
        The ribo-seq profiling is noisy at the single gene level (for example,
        the ribosomal density along a gene may change from a positive number
        to zero and, again, to a positive number). Therefore, gaps in the peak
        are allowed. By default, the gaps there can be no more than 2 consecutive
        amino acids in a peak.''')
    peak_group.add_argument("-p", dest="proportion", required=False, type=float, default=0.2,
                            help='''specify the proportion of gap in each peak (default %(default)s).
        Data with too much noise has low credibility, so it is necessary to
        add a limit to the gap proportion in each peak. By default, the total gaps
        length within a peak cannot exceed 20%%.''')
    peak_group.add_argument("--back", dest="background", required=False, type=int, default=0,
                            help='''use the 5 prime AA to filter the enrichment fold (default %(default)s codon,
        not applicable). Two modes are provided: 
        (1): [30](or other number > 0) Some proteins that bind to the nascent polypeptide 
        chain cannot affect the translation process of the first 30 AA. Theoretically, 
        the RPFs of the unbound area should be smaller than the RPFs of the bound area. 
        (2): [0] However, some proteins only affect the elongation of the ribosome, so 
        the peak may exist in any location.''')
    peak_group.add_argument("--bf", dest="backFold", required=False, action='store_true', default=True,
                            help='''specify the fold of background (default %(default)s).
        Enrichment after background region need greater than before (HSP70 binding model). 
        backFold used to filter the strong binding peak.
        ''')
    
    peak_group.add_argument("--up", dest="upstream", required=False, type=int, default=10,
                            help='''retrieve the upstream sequence of binding peaks, (default %(default)s codon.)''')
    peak_group.add_argument("--down", dest="downstream", required=False, type=int, default=10,
                            help='''retrieve the downstream sequence of binding peaks, (default %(default)s codon.)''')
    
    # arguments for output figures, ratio, peak annotation
    output_group = parser.add_argument_group('Output files arguments')
    output_group.add_argument("-o", dest="output", required=False, type=str, default='results',
                              help='''prefix of output file name (default %(default)s_ratio.txt).
        This prefix will be used as the output of Peak / Ratio / Matlab Script / Bed files.''')
    output_group.add_argument("--all", dest="all", required=False, action='store_true', default=False,
                              help='''output all the peak region. If this option on, all peak 
        regions include overlapped items in same gene are retained. Instead, 
        only the optimal peak for non-overlapping regions is output. (default %(default)s)''')
    output_group.add_argument("--rpm", dest="rpm", required=False, action='store_true', default=False,
                              help="output the rpm of each peak. (default %(default)s)")
    output_group.add_argument("--ratio", dest="ratio", required=False, action='store_true', default=False,
                              help="output the original ratio. (default %(default)s)")
    # output_group.add_argument(
    #     "--script", dest="script", required=False, action='store_true', default=True,
    #     help="output the MATLAB scripts for each figure. (default %(default)s)"
    # )
    output_group.add_argument("--fig", dest="fig", required=False, action='store_true', default=False,
                              help="draw the demo graph of peak scan results. (default %(default)s)."
                                   "This step may takes a lot of time.")

    args = parser.parse_args()
    file_check(args.rpf)
    file_check(args.anno)
    args_print(args)

    return args


def serp_overlap():
    parser = argparse.ArgumentParser(description="This script is used to check the overlap of mock-IP and IP.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-m', dest="mock", required=True, type=str,
                             help="the name of input mock-IP peak file in TXT format.")
    input_group.add_argument('-f', dest="flag", required=False, type=str,
                             help="the name of input flag-IP peak file in TXT format.")
    input_group.add_argument('--om', dest="out_mock", required=False, type=str,
                             help="the prefix of mock peak output file (prefix + .overlap.txt).")
    input_group.add_argument('--of', dest="out_flag", required=False, type=str,
                             help="the prefix of flag peak output file (prefix + .overlap.txt).")
    args = parser.parse_args()
    file_check(args.mock)
    file_check(args.flag)
    args_print(args)

    return args


def serp_retrieve():
    parser = argparse.ArgumentParser(description="This script is used to retrieve the peak sequence.")
    input_group = parser.add_argument_group('Required arguments')

    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument('-p', dest="peak", required=True, type=str,
                             help="the name of input peak file in TXT format.")
    input_group.add_argument('-o', dest="output", required=False, type=str, default='input',
                             help="output sequence file (default: %(default)s.fa).")
    # arguments for the peak region
    parser.add_argument('-u', dest="upstream", required=False, type=int, default=0,
                        help="set the up stream length. (default: %(default)s AA).")
    parser.add_argument('-d', dest="downstream", required=False, type=int, default=0,
                        help="set the down stream length. (default: %(default)s AA).")
    args = parser.parse_args()
    file_check(args.rpf)
    file_check(args.bed)
    args_print(args)

    return args


def serp_retrieve_seq():
    parser = argparse.ArgumentParser(description="This script is used to retrieve the peak sequence.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-p', dest="peak", required=True, type=str,
                             help="the name of input peak file in TXT format.")
    input_group.add_argument('-f', dest="fasta", required=False, type=str, help="the name of fasta sequence.")
    input_group.add_argument('-o', dest="output", required=False, type=str, help="the prefix of output file.")
    # arguments for the peak region
    parser.add_argument('-u', dest="up", required=False, type=int, default=30,
                        help="set the up stream length. (default: %(default)s nt).")
    parser.add_argument('-d', dest="down", required=False, type=int, default=30,
                        help="set the down stream length. (default: %(default)s nt).")

    args = parser.parse_args()
    file_check(args.peak)
    file_check(args.fasta)
    args_print(args)

    return args


def serp_properties():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="This script is used to evaluate the different properties of sequence.",
                                     epilog=textwrap.dedent('''
        Properties contain the following patrs:
        Part1:
            CAI  - codon adaptation index
            RSCU - Relative Synonymous Codon Usage
            Obsi - Observed number of occurrences of codon (per 1000 codon)
        Part2:
            gravy             - according to Kyte and Doolittle
            structure         - fraction of helix, turn and sheet
            flexibility       - according to Vihinen, 1994
            instability       - according to Guruprasad et al 1990
            isoelectric_point - use module IsoelectricPoint
        '''))

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-f', dest="fasta", required=True, type=str,
                             help="the name of input RNA sequence in fasta format.")
    input_group.add_argument('-o', dest="output", required=True, type=str, help="the prefix of output file,.")
    args = parser.parse_args()
    file_check(args.fasta)
    args_print(args)

    return args
