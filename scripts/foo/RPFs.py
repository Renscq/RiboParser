# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : RPFs.py


import sys
import pandas as pd
import polars as pl


def read_rpf_file(rpf_file):
    '''
    @Message  : import the rpf table.
    @Input    : rpf_file --> rpf table file generated from rpf_density.py
    @Return   : raw_rpf --> import rpf table with polars
    '''
    
    raw_rpf = pl.read_csv(rpf_file, has_header = True, separator = '\t')
    raw_rpf = raw_rpf.with_columns(raw_rpf['codon'].str.to_uppercase())
    raw_rpf = raw_rpf.filter(~raw_rpf['codon'].str.contains('N'))

    return raw_rpf


def get_frame_rpf(raw_rpf, sample_name, frame):

    '''
    @Message  : merge all frame density or retrieve one of them
    @Input    : raw_rpf --> rpf table
                merged_rpf --> merged three frame rpf table 
                sample_name --> sample name of rpf table 
                frame --> frame of  [0, 1, 2, all] to retrieve
    @Return   : merged_rpf --> retrieve specific frame of rpf table
    '''
    
    merged_rpf = raw_rpf[:, 0:6]

    if frame == 'all':
        for now_sp in sample_name:
            print('Import the density file {file_name}.'.format(file_name=now_sp), flush=True)
            now_sp_index = [now_sp + '_f0', now_sp + '_f1', now_sp + '_f2']
            
            # merged_rpf = merged_rpf.with_columns(raw_rpf[now_sp_index].sum(axis=1)).rename({now_sp_index[0]: now_sp})
            # report a bug in polars, the sum function is not working.
            merged_rpf = merged_rpf.with_columns(raw_rpf[:, now_sp_index[0]] + raw_rpf[:, now_sp_index[1]] + raw_rpf[:, now_sp_index[2]]).rename({now_sp_index[0]: now_sp})
            
    elif frame == '0':
        for now_sp in sample_name:
            print('Import the density file {file_name}.'.format(file_name=now_sp), flush=True)
            now_sp_index = now_sp + '_f0'
            merged_rpf = merged_rpf.with_columns(raw_rpf[:, now_sp_index]).rename({now_sp_index: now_sp})

    elif frame == '1':
        for now_sp in sample_name:
            print('Import the density file {file_name}.'.format(file_name=now_sp), flush=True)
            now_sp_index = now_sp + '_f1'
            merged_rpf = merged_rpf.with_columns(raw_rpf[:, now_sp_index]).rename({now_sp_index: now_sp})

    elif frame == '2':
        for now_sp in sample_name:
            print('Import the density file {file_name}.'.format(file_name=now_sp), flush=True)
            now_sp_index = now_sp + '_f2'
            merged_rpf = merged_rpf.with_columns(raw_rpf[:, now_sp_index]).rename({now_sp_index: now_sp})

    merged_rpf = merged_rpf.to_pandas()

    return merged_rpf


def set_codon_shift(sites):

    '''
    @Message  : set the codon shift by E/P/A site
    @Input    : sites --> [E, P, A] sites in the ribosome
    @Return   : shift_codon --> [-1, 0, 1] number for [E, P, A] site to shift.
    '''
    
    if sites == 'E':
        shift_codon = -1
    elif sites == 'P':
        shift_codon = 0
    elif sites == 'A':
        shift_codon = 1
    else:
        shift_codon = 0
    return shift_codon


def shift_site(merged_rpf, sample_name, shift_num):
    '''
    @Message  : to shift the E/P/A-site.
    @Input    : merged_rpf --> shift codon from default P-site to E/A site
    @Return   : merged_rpf --> shift to specific coodn site of rpf table
    '''
    
    merged_rpf.loc[:, sample_name] = merged_rpf.groupby('name')[sample_name].shift(shift_num).fillna(0).values
    # merged_rpf = merged_rpf.with_columns(merged_rpf.groupby('name').apply(lambda x: x[sample_name].shift(-1)))
    
    return merged_rpf


def sum_total_rpf(raw_rpf):
    # total_rpf_num = rpf[sample_name].sum()

    total_rpf_num = raw_rpf.select(raw_rpf.columns[6:]).sum().to_pandas().T.reset_index()
    total_rpf_num.columns = ['sample', 'total_rpf']
    total_rpf_num['sample'] = total_rpf_num['sample'].str[:-3]
    total_rpf_num = total_rpf_num.groupby('sample')['total_rpf'].sum()

    return total_rpf_num


def high_expression(merged_rpf, total_rpf_num, sample_name, rpf_num):

    '''
    @Message  : filter the high expression gene.
    @Input    : merged_rpf --> merged rpf table
                sample_name --> sample name of rpf table
                rpf_num --> minimum value for rpf table filter
    @Return   : total_rpf_num --> total rpf of each sample
                gene_rpf_sum --> total rpf of each gene
                high_gene --> retrieved gene list with more than minimum rpf value
                high_rpf --> retrieved gene rpf table with more than minimum rpf value
    '''

    gene_rpf_sum = merged_rpf.groupby('name')[sample_name].sum()
    high_gene = gene_rpf_sum.loc[gene_rpf_sum.mean(axis=1) >= rpf_num, :].index
    high_rpf = merged_rpf.loc[merged_rpf['name'].isin(high_gene), :]

    return gene_rpf_sum, high_gene, high_rpf


def sum_sample_number(raw_rpf, sample_num):
    '''
    @Message  : statistic the number of samples.
    @Input    : raw_rpf --> rpf table in polars dataframe format
                sample_num --> default [None]
    @Return   : sample_num --> number of samples
    '''
    
    if sample_num:
        pass
    else:
        sample_num = int((len(raw_rpf.columns) - 6) / 3)

    return sample_num


def retrieve_sample_name(raw_rpf, sample_name):
    '''
    @Message  : retrieve the rpf table with specific sample name.
    @Input    : sample_name --> input sample name, default [None]
    @Return   : sample_name --> retrieved sample name from rpf table
    '''
    
    if sample_name:
        pass
    else:
        # sample_name = pl.Series(raw_rpf.columns[6::]).apply(lambda s: s[:-3]).unique().to_list()
        sample_name = pd.Series(raw_rpf.columns[6::]).str[:-3].drop_duplicates().to_list()

    return sample_name


def filter_tis_tts(raw_rpf, tis, tts):
    '''
    @Message  : trim the rpf table with tis offset.
    @Input    : raw_rpf --> raw rpf table
                tis --> offset from translation initiation site
                tts --> offset from translation termination site
    @Return   : raw_rpf --> trimmed rpf table
    '''
    
    if not tis is None:
        raw_rpf = raw_rpf.filter(pl.col("from_tis") >= tis)
    else:
        pass

    if not tts is None:
        raw_rpf = raw_rpf.filter(pl.col("from_tts") <= -tts)
    else:
        pass

    return raw_rpf


def retrieve_gene_rpf(raw_rpf, gene):
    '''
    @Message  : retrieve the rpf table with specific gene ids.
    @Input    : raw_rpf --> input rpf in polars dataframe format
                gene --> input gene ids, default [None]
    @Return   : sample_name --> retrieved gene rpf from merged rpf table
    '''

    if gene:
        gene_csv = pd.read_csv(gene, header=0, index_col=None, sep='\t')
        if 'transcript_id' in gene_csv.columns:
            gene_name = gene_csv['transcript_id'].tolist()
        else:
            gene_name = gene_csv.iloc[:, 0].tolist()

        raw_rpf = raw_rpf.filter(pl.col("name").is_in(gene_name))
    else:
        pass
    
    return raw_rpf


def check_merged_rpf_results(merged_rpf):
    if merged_rpf.empty:
        print('Filtered density table is empty! Please recheck the parameters and data quality.', flush=True)
        sys.exit()


def check_high_rpf_results(high_rpf):
    if high_rpf.empty:
        print('Filtered density table is empty! Please recheck the parameters and data quality.', flush=True)
        sys.exit()


def rpf_to_rpm(raw_rpf, norm):
    '''
    @Message  : convert rpf to rpm value.
    @Input    : raw_rpf --> rpf table in polars dataframe
                norm --> default [False]
    @Return   : raw_rpf --> rpm table
    '''
    
    if norm:
        raw_rpf.iloc[:, 6:] = raw_rpf.iloc[:, 6:] * 1000000 / raw_rpf[:, 6:].sum()
    else:
        pass
    
    return raw_rpf


def import_rpf(rpf_file, sample_num, sample_name, sites, frame, gene, rpf_num, tis = None, tts = None):
    ########################################
    # import the rpf file generate by ribo_parser
    raw_rpf = read_rpf_file(rpf_file)
    total_rpf_num = sum_total_rpf(raw_rpf)

    ########################################
    # User can specify a gene list in txt format directly to exclude genes that don't want to analyze.
    # and the gene list like this (no column names are required):

    # YDL246C
    # YDL243C
    # YDR387C
    # ...
    raw_rpf = retrieve_gene_rpf(raw_rpf, gene)

    ########################################
    # filter the CDS range of density file (tis, tts)
    raw_rpf = filter_tis_tts(raw_rpf, tis, tts)
    valid_rpf_num = sum_total_rpf(raw_rpf)

    ########################################
    # retrieve the sample name and numbers
    # columns like this ['sample1_f0', 'sample1_f1', 'sample1_f2']
    sample_num = sum_sample_number(raw_rpf, sample_num)
    sample_name = retrieve_sample_name(raw_rpf, sample_name)

    ########################################
    # merge the RPFs of nucleotide to codon
    # retrieve the codon frame
    # front 6 columns ['name', 'now_nt', 'from_tis', 'from_tes', 'region', 'codon']
    merged_rpf = get_frame_rpf(raw_rpf, sample_name, frame)
    # total_rpf_num = sum_total_rpf(merged_rpf, sample_name)
    

    ########################################
    # RPFs are different at E/P/A site. Default P-site. shift -1/1 codon to E-site/A-site
    shift_num = set_codon_shift(sites)
    merged_rpf = shift_site(merged_rpf, sample_name, shift_num)

    ########################################
    # Filter the high expression genes by specified rpfs num with at least one sample.
    gene_rpf_sum, high_gene, high_rpf = high_expression(merged_rpf, total_rpf_num, sample_name, rpf_num)

    ########################################
    # Check the filtered results. exit if the dataframe is empty.
    check_merged_rpf_results(merged_rpf)
    check_high_rpf_results(high_rpf)

    return raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf


def codon_table():
    codon_anno = {
        'AAA': ['Lys', 'K'], 'AAC': ['Asn', 'N'], 'AAG': ['Lys', 'K'], 'AAT': ['Asn', 'N'],
        'ACA': ['Thr', 'T'], 'ACC': ['Thr', 'T'], 'ACG': ['Thr', 'T'], 'ACT': ['Thr', 'T'],
        'AGA': ['Arg', 'R'], 'AGC': ['Ser', 'S'], 'AGG': ['Arg', 'R'], 'AGT': ['Ser', 'S'],
        'ATA': ['Ile', 'I'], 'ATC': ['Ile', 'I'], 'ATG': ['Met', 'M'], 'ATT': ['Ile', 'I'],
        'CAA': ['Gln', 'Q'], 'CAC': ['HIS', 'H'], 'CAG': ['Gln', 'Q'], 'CAT': ['HIS', 'H'],
        'CCA': ['Pro', 'P'], 'CCC': ['Pro', 'P'], 'CCG': ['Pro', 'P'], 'CCT': ['Pro', 'P'],
        'CGA': ['Arg', 'R'], 'CGC': ['Arg', 'R'], 'CGG': ['Arg', 'R'], 'CGT': ['Arg', 'R'],
        'CTA': ['Leu', 'L'], 'CTC': ['Leu', 'L'], 'CTG': ['Leu', 'L'], 'CTT': ['Leu', 'L'],
        'GAA': ['Glu', 'E'], 'GAC': ['Asp', 'D'], 'GAG': ['Glu', 'E'], 'GAT': ['Asp', 'D'],
        'GCA': ['Ala', 'A'], 'GCC': ['Ala', 'A'], 'GCG': ['Ala', 'A'], 'GCT': ['Ala', 'A'],
        'GGA': ['Gly', 'G'], 'GGC': ['Gly', 'G'], 'GGG': ['Gly', 'G'], 'GGT': ['Gly', 'G'],
        'GTA': ['Val', 'V'], 'GTC': ['Val', 'V'], 'GTG': ['Val', 'V'], 'GTT': ['Val', 'V'],
        'TAC': ['Tyr', 'Y'], 'TAT': ['Tyr', 'Y'], 'TCA': ['Ser', 'S'], 'TCC': ['Ser', 'S'],
        'TCG': ['Ser', 'S'], 'TCT': ['Ser', 'S'], 'TGC': ['Cys', 'C'], 'TGG': ['Trp', 'W'],
        'TGT': ['Cys', 'C'], 'TTA': ['Leu', 'L'], 'TTC': ['Phe', 'F'], 'TTG': ['Leu', 'L'],
        'TTT': ['Phe', 'F'], 'TAA': ['Stop', '*'], 'TAG': ['Stop', '*'], 'TGA': ['Stop', '*']
    }

    codon_df = pd.DataFrame(codon_anno, index=['AA', 'Abbr']).T.sort_values('Abbr')

    codon_df.index.name = 'codon'
    # codon_df = codon_df.reset_index()

    return codon_anno, codon_df
