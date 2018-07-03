#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use POSIX;
use Tie::File;
use Fcntl 'O_RDONLY';

my %PFM = (
	A => [],
	T => [],
	G => [],
	C => [],
	);
my $i = 0;
my $y=0;

#get filename from args
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: name.pl filename output\n";
    exit;
}
my $filename = $ARGV[0];
my $outputFilename = $ARGV[1];
my $outOff=0;
my $STDOUT = "";

#Open file containing sequence with O_RDONLY
tie my @resolvarray, 'Tie::File', $filename, mode => 'O_RDONLY' or die "$0: open $filename: $!";

#	(default is O_CREATE | O_RDWR)
tie my @outputF, 'Tie::File', $outputFilename or die "$0: open $outputFilename: $!";


$outputF[$outOff]= "Reading from file $filename";
$outOff+=1;

#we assume every row must be of same length
#then we take as a basis the length of the 1st row
my $Length = length($resolvarray[0]);

#init PFM with zeros
foreach my $key (keys %PFM){
	for($y=0;$y<$Length;$y++){
		push($PFM{$key},0);
	}
}

#Read every row and add an 'appearence' for
#every base we find to each respective position
#on the hash table(PFM)
for($i=0; $i < scalar(@resolvarray); $i++){
	my $row = $resolvarray[$i];
	chomp $row;
	if($row =~ /[^ATCG]+/){
		$STDOUT = "Error reading file: $row (row: ".$i."+1)\n"."Found illegal char: ".(substr($row, $-[0], $+[0] - $-[0]))."\n";
		print $STDOUT;
		$outputF[$outOff] = $STDOUT;
		untie @resolvarray;
		untie @outputF;
		exit;
	}elsif(length($row)!=$Length){
		$STDOUT = "Error reading file: All rows must be of equal length (row:".$i."+1 expected $Length but got ". length($row).")\n";
		print $STDOUT;
		$outputF[$outOff] = $STDOUT;
		untie @resolvarray;
		untie @outputF;
		exit;
	}
	$y=0;
	for my $c (split //, $row){
		${$PFM{$c}}[$y] += 1;
		$y+=1;
	}
}
untie @resolvarray;

print "Position Frequency Matrix $i (PFM)\n";
$outputF[$outOff] = "Position Frequency Matrix $i (PFM)\n";
$outOff += 1;
foreach my $key (keys %PFM){
	my $message = "$key: {".join(", ", @{$PFM{$key}})."}\n";
	print "$message";
	$outputF[$outOff] = $message;
	$outOff += 1;

}

#Find the most frequent sequence
my @MFS;
for($i=0; $i < $Length; $i++){
	push(@MFS, 0);
	my $value =" ";
	foreach my $key(keys %PFM){
		if( ${$PFM{$key}}[$i] > $MFS[$i]){
			$MFS[$i]=${$PFM{$key}}[$i];
			$value = $key;
		}elsif(${$PFM{$key}}[$i] == $MFS[$i]){
			$key = "|$key";
			$value .= $key;
		}
	}
	$MFS[$i]=$value;
}

print "Most frequent sequence\n",join(" ", @MFS), "\n";
$outputF[$outOff] = "Most frequent sequence\n".join(" ", @MFS)."\n";
untie @outputF;

