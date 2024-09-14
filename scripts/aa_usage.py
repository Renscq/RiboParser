# author: Rhondene Wint
# To compute Amino acid usage of 1 or hundreds of species

import pandas as pd
import numpy as np


def get_seqs(species):
    """species: path of fasta file of the species
    parse fasta file into a list of coding sequences
    Returns: list of each coding sequences
   """
    seqs = []
    with open(species, 'r') as f:
        for line in f:
            if line.startswith('>') is False:
                seqs.append(line.rstrip())
    return seqs


## dictionary that maps codons to amino acids
"""break up 6-codon family into 2 and 4 fold (Ser (S), L(Leu), R (Arg))"""
codon_to_aa = {
    "UUU": "Phe", "UUC": "Phe",
    "UCU": "Ser4", "UCC": "Ser4", "UCA": "Ser4", "UCG": "Ser4",
    "AGU": "Ser2", "AGC": "Ser2",
    "CUU": "Leu4", "CUC": "Leu4", "CUA": "Leu4", "CUG": "Leu4",
    "UUA": "Leu2", "UUG": "Leu2",

    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg4", "CGC": "Arg4", "CGA": "Arg4", "CGG": "Arg4",
    "AGA": "Arg2", "AGG": "Arg2",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",

    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",

    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}


def get_cod_freq(seqs):
    """ seqs: list of CDS
    	Returns a 59-dim dataframe of total absolute codon frequencies
    """

    codon_count = dict()

    for codon in list(codon_to_aa.keys()):
        codon_count[codon] = 0  ##dictionary to accumulate codon count

    for cds in seqs:
        cds = cds.upper().replace('T', 'U')
        codons = []
        ##make a list of codons
        for c in range(0, len(cds), 3):
            if len(cds) % 3 == 0:
                cod = cds[c:c + 3]
                if 'N' not in cod:  ##ignore N and seqs not multiple of 3
                    codons.append(cod)
            else:
                continue

        for c in list(codon_count.keys()):
            codon_count[c] += codons.count(c)

    df_codcount = pd.DataFrame(list(codon_count.items()))
    df_codcount.columns = ['Codon', 'Obs_Freq']
    df_codcount['Amino_Acid'] = [codon_to_aa[codon] for codon in df_codcount['Codon'].values]  ## add amino acid column

    return df_codcount


def compute_rscu_weights(df_codcount):
    """ Caclculates Relative Synonymous codon usage (RSCU) wij = RSCUij/ RSCU i,max
    Input: 59-dim codon count dataframe
    Returns: 59-dim dataframe of RSCU values for each codon """
    aa_groups = df_codcount.groupby('Amino_Acid')
    aa = df_codcount['Amino_Acid'].unique()  # make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d = aa_groups.get_group(a)
        d['RSCU'] = d['Obs_Freq'].values / d['Obs_Freq'].mean()  # obs/expected freq
        d['Relative_Adaptive_Weights'] = d['RSCU'].values / d['RSCU'].max()
        d['optimal'] = [True if rscu == d['RSCU'].max() else False for rscu in d['RSCU'].values]  # marks optimal codon
        d['Species'] = species
        df_list.append(d)
    return pd.concat(df_list)


## save codon usage tables
for species in names:  # names is a list of filenames of fasta files
    seqs = get_seqs(species)  ##formats fasta into list of sequences
    df_rscu = get_cod_freq(seqs)  ##computes absolute codon frequencies
    rscu = compute_rscu_weights(df_rscu)  ##computes RSCU and adaptive weights
    ##saves final table with RSCU, adaptive weights and codon frequencies,
    rscu.to_csv('../RSCU_genomes/{}_all_codons_rscu.csv'.format(species), index=False, sep=',')

###### then compute AA usage from codon usage table######
base_freq = {'U': 0.220, "A": 0.303, 'C': 0.217, 'G': 0.261}  ##natural frequency in nature


## to consolidate all codons of 6-box amino acids; you may or may not want to do this
def no_six(aa):
    return aa[:3]


##very fast , 460 genomes under 30 seconds

# genomes is a list of file names with the RSCU codon usage tables
for species in genomes:
    df_rscu = pd.read_table('../RSCU_genomes_all_codons/{}_all_codons_rscu.csv'.format(species), sep=',')
    df_rscu = df_rscu[df_rscu['Amino_Acid'] != 'STOP']
    df_rscu['AA_no_six'] = df_rscu['Amino_Acid'].apply(no_six)  # comment out if you choose not to
    AA = df_rscu['AA_no_six'].unique()
    AA_group = df_rscu.groupby('AA_no_six')
    aa_usage = dict()
    for amino in AA:

        df = AA_group.get_group(amino)
        ##compute expected aa frequncy
        expected_aa_usage = 0;
        for codon in df['Codon'].values:
            expected_aa_usage += base_freq[codon[0]] * base_freq[codon[1]] * base_freq[codon[2]]
        expected_aa_usage = expected_aa_usage * 1.057  ##correction factor

        ##compute observed aa frequncy (aa_freq/all_codon_freq )
        obs_freq = df['Obs_Freq'].sum()

        obs_freq = obs_freq / df_rscu['Obs_Freq'].sum()
        num_codons = df['Codon'].unique().shape[0]

        aa_usage[amino] = [expected_aa_usage * 100, obs_freq * 100, df['Obs_Freq'].sum(), num_codons]

    ##format into table
    aa_df = pd.DataFrame.from_dict(aa_usage, orient='index').reset_index()
    aa_df.columns = ['Amino_acid', 'Expected_Freq(%)', 'Obs_Freq(%)', 'Abs_Freq', 'Num_Codons']
    aa_df['Species'] = species
    aa_df.to_csv('../AA_usage/{}_aa_usage.csv'.format(species), index=False, sep=',')