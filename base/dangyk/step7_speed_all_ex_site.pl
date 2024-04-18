#!/usr/bin/perl
use strict;
use warnings;
print "\n\nWARNING: this script require large RAM (>10GB). Make sure your computer has at least 12GB RAM or it could freeze your computer\n";
####This script is to call calculate the number of codons and the RPF of each codon in a selected group of genes. This will somewhat reflect the speed of translation at codon base level. 
### this one is different in that it will take codons that has coverage, but if no coverage, it will not be count. Important!!
##set number of codon excluded from both end of 5' end and 3' end, is codon, not nucleotides
my $ex = 20;
my $site = 0; # if A, $site is 0, # if A+1, $site is 1; # if P, $site is -1; # if E, $site is -2; so on.....
#get the CAI reference
my %cai;
open(hand1, "CAI_AA.txt") or die $!;
while (<hand1>) {
    $_ =~ s/\s+$//;
    my ($aa, $co, $rscu, $ca) = split /\t/;
    $cai{$co} = "$aa\t$rscu\t$ca";
}
close hand1;
mkdir "speed";

&speed('RPF_sample_name1', 'mRNA_sample_name1');
&speed('RPF_sample_name2', 'mRNA_sample_name2');
&speed('RPF_sample_name3', 'mRNA_sample_name3');
&speed('RPF_sample_name4', 'mRNA_sample_name4');

sub speed {

    my ($RPF, $mRNAseq) = @_;

    my $sample = $RPF;

    $sample =~ s/\_\d+$//;

    print "now processing $sample\n\n";

    my $ge;
    my %gene; #the ncu number as key and ORF sequence as input
    my %rpf;  #coverage hash for RPF with ncu and start position as key and score as value
    my %mRNA; #ditto
    my %si;   #size of each gene

    open(hand1, "100_reference_100.fa") or die $!;
    while (<hand1>) {
        $_ =~ s/\s+$//;

        if (/^>/) {
            $ge = $_;
            $ge =~ s/^>//;
            next;
        }

        $gene{$ge} = $_;
    }

    close hand1;

    ###############################
    #get the genome length
    ###############################
    foreach my $j (keys %gene) {
        $si{$j} = length($gene{$j});}
    ###############################
    #get data from RPF files and bed file
    ###############################
    open(hand2, "$RPF/$RPF\_A.bg") or die $!;
    open(hand3, "$mRNAseq/$mRNAseq.bed") or die $!;
    while (<hand2>) {
        chomp;
        my @a1 = split /\t/;
        for my $i ($a1[1] .. $a1[2] - 1) {
            my $position = $i + $site * 3;
            $rpf{$a1[0] . "\t" . $i} = $a1[3];
        }
    }
    close hand2;
    # calculate the RPKM of genes
    my $mrna_tot = 0;
    while (<hand3>) {
        chomp;
        $mrna_tot++;
        my @a2 = split /\t/;
        if (exists $si{$a2[0]}) {if ($a2[2] > 99 or $a2[1] < ($si{$a2[0]} - 200) and $a2[5] eq '+') {$mRNA{$a2[0]}++;}}
    }
    close hand3;
    # get fpkm from cufflinks

    my %fpkm = ();
    open(hand4, "FPKM/$mRNAseq\.genes.fpkm_tracking") or die $!; #this file is from sep 12 RNA seq. pair end
    while (<hand4>) {
        chomp;
        my @a3 = split /\t/;
        $fpkm{$a3[0]} = $a3[9];
    } #this is gene: NCU00004
    close hand4;

    ###############################
    #calculate coverage for each codon
    ##########################

    my %codon = ();
    my %codon_12 = ();
    my %rpfcod = ();
    my %rpfcod_norm = ();
    my %rpfcod_fpkm = ();
    my $fold = 0;
    my %can = ();
    my %can_12 = ();
    my $rpkm;

    open(RPKM, ">speed/$mRNAseq\_rpkm_site_$site.txt") or die $!;
    print RPKM "Gene\tlength\t#reads_in_ORF\tRPKM\tRPKM(2012seq)\n";

    foreach my $id (sort keys %si) {
        my $ncu = $id;
        $ncu =~ s/T[0-9]$//;; #print "$id\t$ncu\t";
        $rpkm = 0;            #this RPKM is based on the RNA seq along with the RPF
        if (exists $mRNA{$id}) {
            $rpkm = $mRNA{$id} * 1000 / $si{$id} * (1000000 / $mrna_tot);
            #RPKM= number of reads in gene/length of genes. the number of read in gene needs to multiply with norm factor (1M/total read mapped)
            print RPKM "$id\t$si{$id}\t$mRNA{$id}\t$rpkm\t";

            if (exists $fpkm{$ncu} and $fpkm{$ncu} > 0) {
                print RPKM "$fpkm{$ncu}\n";}
            else {print RPKM "NA\n";}
        }

        if (exists $gene{$id}) {for (my $i = 100 + $ex * 3; $i < ($si{$id} - 100 - $ex * 3); $i += 3) {
            my $cod = substr($gene{$id}, $i, 3);
            my $j = $i + 1;
            if (exists $rpf{$id . "\t" . $j} and $rpf{$id . "\t" . $j} != 0) {

                if ($rpkm >= 5) {
                                                            #RPKM from RPF RNA seq
                    $rpfcod{$cod} += $rpf{$id . "\t" . $i}; #ABS value
                    $codon{$cod}++;
                    $can{$id}++;
                    my $norm_cod = $rpf{$id . "\t" . $j} / $rpkm; #print "$norm_cod\t";
                    $rpfcod_norm{$cod} += $norm_cod;
                }                                                 #(rpf at certain position)/RPKM.

                if (exists $fpkm{$ncu} and $fpkm{$ncu} >= 5) {
                    #RNAseq from sep 2012
                    $codon_12{$cod}++;
                    $can_12{$id}++;
                    my $norm_cod = $rpf{$id . "\t" . $j} / $fpkm{$ncu}; #print "$norm_cod\t";
                    $rpfcod_fpkm{$cod} += $norm_cod;
                }
            }
        }}
    }

    open(hand5, ">speed/$sample\_speed_all_genes_exclude$ex\_site_$site.txt");
    my $can = scalar keys %can;
    my $can_12 = scalar keys %can_12;
    print hand5 "number of genes in $sample included in this test:\nspeed_rpkm=$can\tspeed_fpkm12=$can_12\nnumber of codons excluded from both 5' and 3' end of ORF = $ex\n";
    print hand5 "NOTE: all genes are selected as long as it has at least one codon that have reads. For a specific codon, if there is no RPF coverage, the codon will not be counted!!!! Reads are normalized by either RPKM (calculated by reads in ORF) or by FPKM generated by cufflinks\n";
    print hand5 "----------------------------------------------------------------------------------\n\n";
    print hand5 "Codon\tamino_acid\tRSCU\tCAI\tcodon_number\tspeed_abs\tspeed_rpkm\tspeed_fpkm_12seq\n";
    foreach my $c1 (sort keys %codon) {
        if (exists $cai{$c1}) {
            #print "$c1\t$cai{$c1}\n";
            if (exists $rpfcod{$c1} and $codon{$c1} != 0 and exists $codon_12{$c1} and exists $rpfcod_fpkm{$c1}) {
                my $speed = $rpfcod{$c1} / $codon{$c1};
                $speed = sprintf("%.4f", $speed);
                $speed =~ s/\.?0+$//;

                exists $rpfcod_norm{$c1} and exists $rpfcod_fpkm{$c1} and

                    my $speed_norm = $rpfcod_norm{$c1} / $codon{$c1};
                $speed_norm = sprintf("%.4f", $speed_norm);
                $speed_norm =~ s/\.?0+$//;

                my $speed_fpkm = $rpfcod_fpkm{$c1} / $codon_12{$c1};
                $speed_fpkm = sprintf("%.4f", $speed_fpkm);
                $speed_fpkm =~ s/\.?0+$//;
                print hand5 "$c1\t$cai{$c1}\t$codon{$c1}\t$speed\t$speed_norm\t$speed_fpkm\n";
            }}}
    close hand5;
}
