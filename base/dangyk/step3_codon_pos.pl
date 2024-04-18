use strict;
use warnings;
#use Statistics::Histogram;
use Number::Format qw/round/;

my $up = 35;
my $low = 25;

mkdir 'codon_stat';

my @sample = ('sample_name1', 'sample_name2', 'sample_name3', 'sample_name4');

foreach my $s (@sample) {
	sta($s);   }

sub sta                                                      {

	my ($sam)=@_;

	my %count=(); my $count=0; my %p12=(); my %p13=(); my %p14=(); my %size=(); my %strand=();

	open (hand1, "$sam/$sam.bed") or die $!;

	while (<hand1>)                                             {

	$_=~ s/\s+$//;
	my @a1=split /\t/;

	$strand{$a1[5]}++;

	if ($a1[5] eq '-') {next;}

	my $size=$a1[2]-$a1[1];

	$size{$size}++;

	if ($size<$low or $size>$up) {next;}

	$count{$size}++;
	$count++;

	my $f11=($a1[1]-11-100)%3;         #assume P site at 12nt
	my $f12=($a1[1]-12-100)%3;
	my $f13=($a1[1]-13-100)%3;

	if ($f11==0)  {$p12{$size}++;}
	if ($f12==0)  {$p13{$size}++;}
	if ($f13==0)  {$p14{$size}++;}                                     }
	close hand1;

	open (hand2, ">codon_stat/$sam.txt");

	foreach (keys %strand)  {
		print hand2 "$_\t$strand{$_}\n"; }

	foreach (sort keys %size)  {
		print hand2 "$_\t$size{$_}\n"; }

	print hand2 "The mapped reads of $low-$up nt are $count\n";

	print hand2 "length\tnumber_of_reads\t12th\t13th\t14th\n";

	for my $si ($low..$up)  {

		my $ra12=round($p12{$si}/$count{$si}*100);

		my $ra13=round($p13{$si}/$count{$si}*100);

		my $ra14=round($p14{$si}/$count{$si}*100);

		print hand2 "$si\t$size{$si}\t$ra12\t$ra13\t$ra14\n" if exists $size{$si};

                         }

	print "$sam is done\n";

	close hand2;                                                          }
    

