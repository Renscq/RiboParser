#!/bin/sh

# Print start status message.
echo "job started"
start_time=`date +%s`

for input in 'sample_name1' 'sample_name2' 'sample_name3' 'sample_name4'

do
bedtools genomecov -split -i $input/$input\_A.bed -d -strand + -g 100_reference_100.chrom > $input/$input\_refA.bg
echo "$input for ref is done"
#this part makes the graph for IGV
bedtools genomecov -split -i $input/$input\_A.bed -bga -strand + -g 100_reference_100.chrom > $input/$input\_A.bg
echo "$input for igv is done"
done

echo "ribosomal profiling data are all done!"

# Print end status message.

echo "job finished"
end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
