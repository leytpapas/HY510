#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use POSIX;

my $frame = "";
my $i = 0;
my $filename = 'yersinia_genome.fasta';
# my $reg_exp  = '/([TA][AC]AGGA[GA][GA])([ATCG]{4,10})(ATG)([ATCG]{3})*?(TGA|TAA|TAG)/g';
my $complement;

#reading whole sequence from file, bypassing 1st line (details)
open(my $fh, '<:encoding(UTF-8)', $filename);
LOOP: while (my $row = <$fh>) {
	if($i == 0){
		$i+=1;
		print $row;
		next LOOP;
	}
	chomp $row;
	$frame .= $row;	
}

$complement = reverse($frame);
$complement =~ tr/ATGC/TACG/;

my $sum;
my $numMatches;
($sum, $numMatches) = sequenceAnalyse($frame, "normal");
(my $temp1, my $temp2) = sequenceAnalyse($complement, "reverse complement");
$sum+=$temp1;
$numMatches+=$temp2;
print "\nTotal triplets:",$sum," ,in ", $numMatches," matches\n";


sub sequenceAnalyse{

	my $sequence = $_[0];
	my $type = $_[1];
	my $sum =0;
	my $i=0;

	print "\n------Analysing $type sequence------\n";
	
	while ( $sequence =~ /([TA][AC]AGGA[GA][GA])([ATCG]{4,10})(ATG)([ATCG]{3})*?(TGA|TAA|TAG)/g ) {
		my $string = substr($sequence , $-[0], $+[0] - $-[0]);
		my $start = $-[0];
		my $end = $+[0];
		
		$string =~ /(ATG)([ATCG]{3})*(TGA|TAA|TAG)$/;
		$start += $-[0];
		# print "Beginning/end :\n\tStarts@ $start(row:,", $start//71 + 1,")\n\tEnds@ $end\n";
		# print "without start/end codons: ",substr($sequence, $start+3, $end - $start - 6 ),"\n"; #without start/end codons
		# print "with start/end codons: ",substr($sequence, $start, $end - $start ),"\n\n"; #with start/end codons
		$sum += $end - $start - 6; #sum is without start/end codons
		$i += 1;
	}
	print "------Analysing $type sequence COMPLETE------\n";

	return $sum ,$i;
}