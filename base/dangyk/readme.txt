These scripts are created by Yunkun and used for analysis of ribosome profiling data.

Before start
require Perl(5.22) and R(3.4.0) or above
1, download and install bowtie2, tophat2, SAMTOOLS, BEDTools. make sure the installed softwares are in your PATH.
2, create the reference files (genomic references, ORF sequences, CAI and CBI values, GTF file and so on).


step 1: This script is used to trim the 3' linker or adaptor . 
step 2: This script is for the removal of contamination£ºincluding removal of rRNA and tRNA sequences, mapping the clean reads to reference, and conversion of file format. 
step 3: This script is used to assess the right 5'end flank of RPF protected fragment, which determines the codon within P site.
step 4: This script is used to trim 5'end flank of RPF protected fragment, which reveals the codon on P site and used to calculate coverage. It will also produce the 5' end prefernce of 4 nt, 2nt ahead of 5' end and 2 nt after the cleavage site, which will provide a statistics for end cleavage preference.
step 5: This script is used for conversion of file format.
step 6: This script is prepared for calculating codon occupancy. 
step 7: This script is used to normalize the ribosome profiling data by library size (total reads).
step 8: This script is used to calculate the absolute codon occupancy. RPF density of a given codon was normalized to library size (total reads) and the mRNA level (FPKM).
step 9: This script is prepared for the calculation of the relative codon occupancy. 1, avergage RPF density of a given gene = total in-frame reads / total effective codon (with PRF reads); 2, for each gene, RPF density of a given codon is normalized with the the average RPF density of the gene.
