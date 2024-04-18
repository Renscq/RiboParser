from multiprocessing import Pool
import numpy as np
import pandas as pd
import pysam

def fetch_reads(args):
    bam_in_file, bam_out_file, gene_list = args
    with pysam.AlignmentFile(bam_in_file, 'rb') as bam_in:
        with pysam.AlignmentFile(bam_out_file, 'wb', template=bam_in) as bam_out:
            for gene in gene_list:
                reads = bam_in.fetch(gene)
                for read in reads:
                    # filter the mapped reads or not
                    if reads.is_unmapped:
                        continue
                    # filter the unique reads
                    if reads.has_tag('NH') and reads.get_tag('NH') <= 0:

                        if reads.reference_name in gene_list:
                            read_length = reads.infer_read_length()
                            if reads.is_reverse:
                                try:
                                    self.length_dict[read_length][1] += 1
                                except KeyError:
                                    self.length_dict[read_length] = [0, 1]
                            else:
                                try:
                                    self.length_dict[read_length][0] += 1
                                except KeyError:
                                    self.length_dict[read_length] = [1, 0]
                                # self.mrna_dict[reads.reference_name].bam.append(reads)
                                if reads.query_name not in bam_seq_dict:
                                    bam_seq_dict[reads.query_name] = set()

                                bam_seq_dict[reads.query_name].add(reads.reference_name)
                                bam_out.write(reads)

def main():
    bam_in_file = '/mnt/t64/rensc/herr/hsa/pc-b4-20230731/test/R_PC_NC_1.bam'
    bam_out_file = '/mnt/t64/rensc/herr/hsa/pc-b4-20230731/test/R_PC_NC_1_flt'
    bam_in_idx_file = '/mnt/t64/rensc/herr/hsa/pc-b4-20230731/test/R_PC_NC_1.idx'

    bam_in_idx = pd.read_csv(bam_in_idx_file, header=None, sep='\t', index_col=None)
    bam_in_idx_list = bam_in_idx.iloc[:, 0].tolist()

    # Split the gene list into 4 parts
    splits = np.array_split(bam_in_idx_list, 4)

    pool = Pool(processes=4)
    args = [(bam_in_file, bam_out_file + str(i) + '.bam', split) for i, split in enumerate(splits)]
    pool.map(fetch_reads, args)

if __name__ == "__main__":
    main()
