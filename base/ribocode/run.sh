genePredToGtf file -utr -honorCdsStat -source=ribocode test.genepred test.gtf
gtfToGenePred -allErrors -genePredExt -ignoreGroupsWithoutExons test.gtf test.gp
