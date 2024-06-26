#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Project      : riboParser
@Script       : get_negative_ppi.py
@Environment  : python 3.8.5
@Version      : 1.0
@Author       : Rensc 
@Time         : 2024/05/27 10:17:44
@E-mail       : rensc0718@163.com
@License      : (C)Copyright 2023-2025, Ren Shuchao
'''


import pandas as pd
# import polars as pl
import numpy as np




# # import the ppi data and save to a dictionary

# def import_ppi(ppi_file, number = 10000):
    
#     ppi_dict = {}

#     # read the file and save to the dictionary
#     with open (ppi_file, 'r') as ppi:
#         flag = 0

#         for lines in ppi:
#             rec = lines.strip().split('_')
            
#             try:
#                 ppi_dict[rec[0]].append(rec[1])
#             except KeyError:
#                 ppi_dict[rec[0]] = [rec[1]]

#             flag += 1
#             if flag > number:
#                 break
    
#     # # convert the list to set
#     for key, value in ppi_dict.items():
#         ppi_dict[key] = set(value)

#     return ppi_dict


# ppi_dict = import_ppi('/mnt/t64/database/gmx/ncbi/ppi/merged_paired_set.txt', number = 10000)





# import the ppi data and save to a dictionary
def import_ppi(ppi_file):
    
    ppi_dict = {}

    # read the file and save to the dictionary
    with open (ppi_file, 'r') as ppi:
        flag = 0

        for lines in ppi:
            rec = lines.strip().split('_')
            
            try:
                ppi_dict[rec[0]].append(rec[1])
            except KeyError:
                ppi_dict[rec[0]] = [rec[1]]
            
            flag += 1
            if flag % 5000000 == 0:
                print(flag)

    # convert the list to set
    flag2 = 0
    for key, value in ppi_dict.items():
        ppi_dict[key] = set(value)

        flag2 += 1
        if flag2 % 1000 == 0:
            print(flag2)

    return ppi_dict

ppi_dict = import_ppi(ppi_file = '/mnt/t64/database/gmx/ncbi/ppi/merged_paired_set.txt')





# import all the uiprot ids
# /mnt/t64/database/gmx/ncbi/ppi/uniprotkb_glycine_max_AND_taxonomy_id_3_2024_05_09.tsv
uniprot_ids = pd.read_csv("/mnt/t64/database/gmx/ncbi/ppi/uniprotkb_glycine_max_AND_taxonomy_id_3_2024_05_09.tsv", sep="\t", low_memory=False)

# uniprot_ids.iloc[0:5, 0:5]
Entry = uniprot_ids['Entry'].values
GeneID = uniprot_ids['GeneID'].values
EnsemblPlants = uniprot_ids['EnsemblPlants'].values




# create the random pairs
random_pairs = {}
flag3 = 0
for i in range(len(Entry)):
    flag3 += 1
    if flag3 % 1000 == 0:
        print(flag3)

    protein1 = Entry[i]
    rand_index = np.random.choice(len(Entry), 2000)
    protein2 = set(Entry[rand_index])
    try:
        protein2 = protein2 - ppi_dict[protein1]
    except KeyError:
        pass
    random_pairs[protein1] = protein2





# create the dataframe and save the random pairs
random_pairs_df = pd.DataFrame.from_dict(random_pairs, orient='index').stack().reset_index(level=1, drop=True).reset_index()
# random_pairs_df

# save the random pairs
random_pairs_df.to_csv("/mnt/t64/database/gmx/ncbi/ppi/random_pairs.txt", sep="\t", index=False, header=False)
