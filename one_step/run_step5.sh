###########################################################
#         Pipeline of the Riboparser (Ribo-seq)           #
###########################################################
#
# This step is used for analyzing Ribo-seq data, utilizing 
# RiboParser to check the sequencing quality of the Ribo-seq
# data.
# 
# **Note!** 
#
# The BAM files, parameters and reference genome 
# information files in Step 0.0 may need to be modified 
# according to the files defined for your project!
#
###########################################################
# 0.0 set the paramater of the riboparser
# set the RPFs length range
min_length=27
max_length=33

# set the paramater for the offset detection
offset_mode="RSBM" # [RSBM, SSCBM]
expect_rpf=30
shift_nt=1

# set the minimum RPFs for each gene
min_reads=50

#sites="A" # [E, P, A]
tis=15
tts=5

# 
barplot="bar" # [bar, line]
genebody="10,150,10" # [utr5,cds,utr3]
background=0
frame=0 # [0, 1, 2]
around=20

##################################################################
# 1.0 run the rpf_Check
mkdir ./4.ribo-seq/5.riboparser
mkdir ./4.ribo-seq/5.riboparser/01.qc
cd ./4.ribo-seq/5.riboparser/01.qc/

for bam in ../../3.star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 12 -t ../../../1.reference/norm/gene.norm.txt \
  -o $prefix_name &> $prefix_name".log"

done

## merge the ribo-seq quality results
merge_length -l *length_distribution.txt -o RIBO
merge_saturation -l *gene_saturation.txt -o RIBO

cd ../
##################################################################
# 2.0 Enzymatic bias in NGS library preparation
mkdir 02.digestion
cd 02.digestion

for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m "$min_length" -M "$max_length" --scale \
 -s ../../../1.reference/norm/gene.norm.rna.fa \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

merge_digestion -l *pwm.txt -o RIBO

cd ../
##################################################################
# 3.0 Use RiboParser for quality checks
mkdir 03.offset
cd 03.offset

for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Offset -b $bam -m "$min_length" -M "$max_length" -p 30 -d \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

merge_offset_detail -l *end.txt -o RIBO
merge_offset -l *RSBM_offset.txt -o RIBO_RSBM
merge_offset -l *SSCBM_offset.txt -o RIBO_SSCBM

cd ../
##################################################################
# 4.0 Convert the bam file to reads density
mkdir 04.density
cd 04.density

for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Density -b $bam -m "$min_length" -M "$max_length" --period 40 -l --thread 12 \
 -p ../03.offset/$prefix_name"_"$offset_mode"_offset.txt" \
 -s ../../../1.reference/norm/gene.norm.rna.fa \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

cd ../
##################################################################
# 5.0 Integrating density file for all samples
mkdir 05.merge
cd 05.merge

merge_dst_list -l ../04.density/*_rpf.txt -o RIBO.file.list

rpf_Merge -l RIBO.file.list -o RIBO &> RIBO.log

cd ../
##################################################################
# 6.0 Check the tri-nucleotide periodicity of Ribo-seq
mkdir 06.periodicity
cd 06.periodicity

rpf_Periodicity \
 -r ../05.merge/RIBO_merged.txt \
 -m "$min_reads" --tis "$tis" --tts "$tts" -o RIBO &> RIBO.log

cd ../
##################################################################
# 7.0 Meta-gene analysis of Ribo-seq
mkdir 07.metaplot
cd 07.metaplot

rpf_Metaplot \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m "$min_reads" --mode "$barplot" -o RIBO &> RIBO.log

cd ../
##################################################################
# 8.0 Check gene density of Ribo-seq
mkdir 08.coverage
cd 08.coverage

rpf_Coverage \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m "$min_reads" --outlier \
 -b "$genebody" \
 -n --heat \
 -o RIBO &> RIBO.log

cd ../
##################################################################
# 9.0 Check the repeatability of Ribo-seq
mkdir 09.correlation
cd 09.correlation

rpf_Corr -r ../05.merge/RIBO_merged.txt -o RIBO &> RIBO.log

cd ../
##################################################################
# 10.0 Quantification of Ribo-seq
mkdir 10.quantification
cd 10.quantification

rpf_Quant -r ../05.merge/RIBO_merged.txt --tis "$tis" --tts "$tts" -o RIBO &> RIBO.log

cd ../
##################################################################
# 11.0 Calculate codon pausing score
mkdir 11.pausing_score
cd 11.pausing_score

for sites in E P A
do
rpf_Pausing \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -b 0 --stop \
 -m "$min_reads" \
 -s $sites \
 -f 0 \
 --scale minmax \
 -o "$sites"_site &> "$sites"_site.log
done

cd ../
##################################################################
# 12.0 Calculate codon-level occupancy in Ribo-seq data
mkdir 12.codon_occupancy
cd 12.codon_occupancy

for sites in E P A
do
rpf_Occupancy \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m "$min_reads" \
 -s "$sites" \
 -f 0 --stop \
 --scale minmax \
 -o "$sites"_site &> "$sites"_site.log
done

cd ../
##################################################################
# 13.0 Calculate codon-level decoding time in Ribo-seq data
mkdir 13.codon_decoding_time
cd 13.codon_decoding_time

for sites in E P A
do
rpf_CDT \
 -l ../../../1.reference/norm/gene.norm.txt \
 --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
 --rpf ../05.merge/RIBO_merged.txt \
 --stop \
 -m "$min_reads" \
 -f 0 \
 -s $sites \
 --tis "$tis" --tts "$tts" \
 -o "$sites"_site &> "$sites"_site.log
done

cd ../
##################################################################
# 14.0 Calculate codon-level selection time in Ribo-seq data
mkdir 14.codon_selection_time
cd 14.codon_selection_time

for sites in E P A
do
rpf_CST \
 -l ../../../1.reference/norm/gene.norm.txt \
 --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
 --rpf ../05.merge/RIBO_merged.txt \
 --stop \
 -m "$min_reads" \
 -f 0 \
 -s $sites \
 --tis "$tis" --tts "$tts" \
 -o "$sites"_site &> "$sites"_site.log
done

cd ../
##################################################################
# 15.0 Calculate gene coefficient of variation in Ribo-seq data
mkdir 15.coefficient_of_variation
cd 15.coefficient_of_variation

rpf_CoV \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -f 0 \
 -m "$min_reads" \
 --tis "$tis" --tts "$tts" \
 --fig \
 -g ../../../design.txt \
 -o RIBO &> RIBO.log

cd ../
##################################################################
# 16.0 Calculate meta-codon density in Ribo-seq data
mkdir 16.meta_codon
cd 16.meta_codon

rpf_Meta_Codon \
 -r ../05.merge/RIBO_merged.txt \
 -m "$min_reads" -f 0 \
 -c ../../../codon_list.txt \
 -a 15 -u -n \
 -o RIBO &> RIBO.log
