genePredToGtf file -utr -honorCdsStat -source=ribotish test.genepred test.gtf
gtfToGenePred -allErrors -genePredExt -ignoreGroupsWithoutExons test.gtf test.gp
