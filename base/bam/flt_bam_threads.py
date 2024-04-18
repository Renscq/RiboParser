#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Quality.py


import os.path


from collections import Counter
from collections import OrderedDict
from itertools import chain
from itertools import islice

import numpy as np
import pandas as pd
import pysam


bam_in_file = '/mnt/t64/rensc/herr/hsa/pc-b4-20230731/test/R_PC_NC_1.bam'
bam_out_file = '/mnt/t64/rensc/herr/hsa/pc-b4-20230731/test/R_PC_NC_1_flt'
bam_in_idx_file = '/mnt/t64/rensc/herr/hsa/pc-b4-20230731/test/R_PC_NC_1.idx'

bam_in_idx = pd.read_csv(bam_in_idx_file, header=None, sep='\t', index_col=None)
bam_in_idx_list = bam_in_idx.iloc[:, 0].tolist()

# split the gene list to 4 parts
splits = np.array_split(bam_in_idx_list, 4)
splits = [list(split) for split in splits]

# import the bam file
bam_in1 = pysam.AlignmentFile(bam_in_file, 'rb')
bam_in2 = pysam.AlignmentFile(bam_in_file, 'rb')
bam_in3 = pysam.AlignmentFile(bam_in_file, 'rb')
bam_in4 = pysam.AlignmentFile(bam_in_file, 'rb')

# set the output multiple bam file
bam_out1 = pysam.AlignmentFile(bam_out_file + '1.bam', 'wb', template=bam_in1)
bam_out2 = pysam.AlignmentFile(bam_out_file + '2.bam', 'wb', template=bam_in2)
bam_out3 = pysam.AlignmentFile(bam_out_file + '3.bam', 'wb', template=bam_in3)  
bam_out4 = pysam.AlignmentFile(bam_out_file + '4.bam', 'wb', template=bam_in4)

# fetch the output reads with multi-threads
from multiprocessing import Pool

def fetch_reads(bam_in, bam_out, gene_list):
    # fetch reads
    for gene in gene_list:
        reads = bam_in.fetch(gene)
        # write reads
        for read in reads:
            bam_out.write(read)

pool = Pool(processes=4)
pool.map(fetch_reads, [(bam_in1, bam_out1, splits[0]), (bam_in2, bam_out2, splits[1]), (bam_in3, bam_out3, splits[2]), (bam_in4, bam_out4, splits[3])])
pool.close()
pool.join()

# close the bam file
bam_in1.close()
bam_in2.close()
bam_in3.close()
bam_in4.close()
bam_out1.close()
bam_out2.close()
bam_out3.close()
bam_out4.close()



bam_parts = []
for i in range(10):
    bam_parts.append(pysam.AlignmentFile(bam_out_file + str(i) + '.bam', 'wb', template=bam_in))
