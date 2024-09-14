#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboparser
# @Script  : rpm_smooth.py


import argparse
import textwrap
import time
import sys
import os

import numpy as np
import pandas as pd
import polars as pl
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter


def now_time():
    """
    1. print the time
    """
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)


def args_print(args):
    """
    1. print the arguments
    """
    args_dict = vars(args)

    for k, v in args_dict.items():
        print("{:<12}:  {:<}".format(k, str(v)), flush=True)


def file_check(*files):
    """
    1. check the input files
    """
    for my_file in files:
        if os.path.exists(my_file):
            continue
        else:
            print("\nFile {file_name} is not exists!\n".format(file_name=my_file), flush=True)
            sys.exit()


def rpm_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to smooth the rpm density of bedgraph.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-i', dest="input", required=True, type=str,
                             help="input bedgraph file name.")
    input_group.add_argument('-o', dest="output", required=True, type=str,
                             help="output bedgraph file name.")
    
    # arguments for the RPM calculation
    parser.add_argument('-m', dest="method", required=False, choices=["SG", "Gaus"], 
                        type=str, default="SG",
                        help="the algorithm for smooth ['Savitzky-Golay', 'Gaussian'] (default: %(default)s).")
    parser.add_argument('-w', dest="width", required=False, type=int, default=11,
                        help="width for windows scanning (only for SG, default: %(default)s ).")
    parser.add_argument('-p', dest="power", required=False,  type=int, default=3,
                        help="power of polynomial (only for SG, default: %(default)s ).")
    parser.add_argument('-s', dest="sigma", required=False,  type=int, default=1,
                        help="sigma of standard deviation (only for Gaus, default: %(default)s ).")
    parser.add_argument('-z', dest="zero", required=False,  action="store_true", default=False,
                        help="discard the site without coverage (default: %(default)s ).")
    args = parser.parse_args()

    file_check(args.input)

    return args


class BedGraph(object):

    def __init__(self, args):

        self.input = args.input
        self.output = args.output
        
        self.method = args.method
        self.width = args.width
        self.power = args.power
        self.sigma = args.sigma

        self.zero = args.zero

        self.bg = None
        self.smooth_bg = None

    def import_file(self):
        """
        1. import the bedgraph data
        2. rename the title of bedgraph
        3. drop the NA and INF value
        4. convert the float to integer type
        """
        new_columns = ["chr", "start", "end", "rpm"]

        self.bg = pl.read_csv(self.input, separator='\t', has_header=False)
        self.bg = self.bg.to_pandas()
        self.bg = self.bg.rename(columns=dict(zip(self.bg.columns, new_columns)))
        
        self.bg.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.bg.dropna(inplace=True)
        
        self.bg[["start", "end"]] = self.bg[["start", "end"]].astype(int)

    def expand_bin(self):
        """
        1. expand the range density to each nucleotide
        eg:
        chr start end rpm
        1   20  23  1.5
        
        to

        chr start end rpm
        1   20  21  1.5
        1   21  22  1.5
        1   22  23  1.5

        """

        expand_site = []
        def expand_rows(row):
            chrom = row["chr"]
            start = row["start"]
            end = row["end"]
            rpm = row["rpm"]

            for i in range(start, end):
                expand_site.append([chrom, i, i + 1, rpm])
        
        _ = self.bg.apply(expand_rows, axis=1)

        self.bg = pd.DataFrame(expand_site, columns=["chr", "start", "end", "rpm"])

    def smooth_rpm_SG(self):
        """
        1. get the chromosome list
        
        2. fill the empty site with 0
        eg:
        chr start end rpm
        1   22  23  1.5
        1   25  26  2.6
        
        to

        chr start end rpm
        1   22  23  1.5
        1   23  24  0
        1   24  25  0
        1   25  26  2.6

        3. smooth the data rpm with Savitzky-Golay

        4. discard the smoothed data with rpm less than 0
        """
        if self.width < 3:
            self.width = 3
        else:
            pass

        chrom_bg_list = []
        chrom_list = self.bg["chr"].unique()

        for chrom in chrom_list:
            print("Processing Chromosome:" + str(chrom), flush=True)

            chrom_bg = self.bg[self.bg["chr"] == chrom]
            new_index = range(0, chrom_bg["end"].max())
            chrom_bg = chrom_bg.set_index('start').reindex(new_index, fill_value=0).reset_index()
            chrom_bg = chrom_bg[["chr", "start", "end", "rpm"]]
            chrom_bg["end"] = chrom_bg["start"] + 1
            chrom_bg["chr"] = chrom
            chrom_bg.loc[:, "rpm"] = savgol_filter(chrom_bg['rpm'], self.width, self.power, mode='nearest')
            chrom_bg["rpm"] = chrom_bg["rpm"].round(4).clip(lower=0)
            chrom_bg_list.append(chrom_bg)

        self.smooth_bg = pd.concat(chrom_bg_list)

    def smooth_rpm_Gaus(self):
        """
        1. get the chromosome list
        
        2. fill the empty site with 0
        eg:
        chr start end rpm
        1   22  23  1.5
        1   25  26  2.6
        
        to

        chr start end rpm
        1   22  23  1.5
        1   23  24  0
        1   24  25  0
        1   25  26  2.6

        3. smooth the data rpm with Gaussian fitting

        4. discard the smoothed data with rpm less than 0
        """
        chrom_bg_list = []
        chrom_list = self.bg["chr"].unique()

        for chrom in chrom_list:
            print("Processing Chromosome: " + str(chrom), flush=True)  
            
            chrom_bg = self.bg[self.bg["chr"] == chrom]
            new_index = range(0, chrom_bg["end"].max())
            chrom_bg = chrom_bg.set_index('start').reindex(new_index, fill_value=0).reset_index()
            chrom_bg = chrom_bg[["chr", "start", "end", "rpm"]]
            chrom_bg["end"] = chrom_bg["start"] + 1
            chrom_bg["chr"] = chrom
            chrom_bg.loc[:, "rpm"] = gaussian_filter(chrom_bg['rpm'], self.sigma)
            chrom_bg["rpm"] = chrom_bg["rpm"].round(4).clip(lower=0)
            chrom_bg_list.append(chrom_bg)

        self.smooth_bg = pd.concat(chrom_bg_list)

    def output_file(self):
        """
        1. remove the lines with rpm = 0
        2. output smoothed bedgraph
        """
        if self.zero:
            self.smooth_bg = self.smooth_bg[self.smooth_bg['rpm'] != 0]
        else:
            pass

        self.smooth_bg = pl.from_dataframe(self.smooth_bg)
        self.smooth_bg.write_csv(self.output, separator='\t', has_header=False)

def main():
    now_time()
    
    print('Step1: Checking the input Arguments.', flush=True)
    args = rpm_args_parser()
    args_print(args)
    bg = BedGraph(args)

    print('Step2: improt the bedgraph file.', flush=True)
    bg.import_file()

    print('Step3: expand the bins density to each nucleotide.', flush=True)
    bg.expand_bin()

    print('Step4: smooth the rpm data.', flush=True)
    if bg.method == "SG":
        bg.smooth_rpm_SG()
    else:
        bg.smooth_rpm_Gaus()

    print('Step5: output the bedgraph file.', flush=True)
    bg.output_file()

    print('All done.', flush=True)
    now_time()


if __name__ == '__main__':

    main()
