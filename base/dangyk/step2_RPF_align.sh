echo "job started"
start_time=`date +%s`
mkdir 'de_r+tRNA_stats'

for input in 'sample_name1' 'sample_name2' 'sample_name3' 'sample_name4'

do
mkdir $input
cd $input
mkdir 'de_r+tRNA'
cd ..

bowtie2 -p 8 --un=$input/de_r+tRNA/$input\_norRNA.fq -x /xxxxx/ribosomal_RNA_reference /xxxxx/fq/$input\_15.fq 2>> de_r+tRNA_stats/$input\_rRNA_stats_r1.txt > $input/de_r+tRNA/$input\_rRNAalign.aln

echo "$input rRNA removal finished"

bowtie2 -p 8 --un=$input/de_r+tRNA/$input\_no_t+rRNA.fq -x /xxxxx/tRNA_reference /xxxxx/$input/de_r+tRNA/$input\_norRNA.fq 2>> de_r+tRNA_stats/$input\_tRNA_stats_r1.txt > $input/de_r+tRNA/$input\_tRNAalign.aln

echo "$input tRNA removal finished"

bowtie2 -p 8 --un $input/b2.umapped.fq -x /xxxxx/100_reference_100 $input/de_r+tRNA/$input\_no_t+rRNA.fq -S $input/$input.sam 2>> de_r+tRNA_stats/$input.txt

cd $input

samtools view -bS $input.sam > $input.bam
samtools sort $input.bam $input.sorted
samtools index $input.sorted.bam
bedtools bamtobed -i $input.sorted.bam > $input.bed
sort -k1,1 -k2,2n $input.bed > $input.sort.bed
wc -l $input.sort.bed > $input.mapped.txt

cd ..
echo "$input alignment finished"
done
echo "job finished"
end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
