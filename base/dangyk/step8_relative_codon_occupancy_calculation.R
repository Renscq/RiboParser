library(dplyr)
library(reshape2)
library(ggplot2)

setwd("X:\XXXX")

wt = "sample_name_wt"
mut = "sample_name_mut"
A = "A"

  s1 = read.table(paste(wt, "/", wt, "_site_", A, "_occupancy_norm1.txt", sep = ""), sep = "\t", header = T)
  s1$norm_reads = s1$reads*1000000/sum(s1$reads)
 
  d1 = s1 %>% group_by(codon) %>% summarize(S1_oc = mean(AVG_density))
  
  s4 = read.table(paste(mut, "/", mut, "_site_", A,  "_occupancy_norm1.txt", sep = ""), sep = "\t", header = T)
  s4$norm_reads = s4$reads*1000000/sum(s4$reads)
  
  d4 = s4 %>% group_by(codon) %>% summarize(S4_oc = mean(AVG_density))
  
  d = merge(d1, d4, by.x = "codon", by.y = "codon")
  
  d$S4_S1 = d$S4_oc/d$S1_oc
  
  cod = read.table("CAI_AA.txt", sep = "\t")
  cod = cod[,1:2]
  names(cod) = c("AA", "codon")
  
  d = merge(cod, d, by.x = "codon", by.y = "codon")
  write.csv(d, paste(wt, mut, A, "AVG.csv", sep = "_"))
  
 
