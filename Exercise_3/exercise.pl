#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use POSIX;
use Tie::File;
use Fcntl 'O_RDONLY';
use Storable qw(dclone);


#get filename from args

#Open file containing sequence with O_RDONLY
print "Enter 1st sequence\n";

my $firstSeq = "ACGGTAG";#<>;
chomp $firstSeq;
print "$firstSeq\n";

print "Enter 2nd sequence\n";
my $secondSeq = "CCTAAG";#<>;
chomp $secondSeq;
print "$secondSeq\n";
my $match = 1; #int
chomp $match;
my $mismatch = -1; #int
chomp $mismatch;
my $gap = -2; #int
chomp $gap;

#init hash mpa. We implement the 2D array as a hash map with keys (i,j)\
print "Length of 1st ",length($firstSeq),"\n";
print "Length of 2nd ",length($secondSeq),"\n";
my %newHash;
for(my $i=0; $i<length($firstSeq)+1; $i++){
	# print"i:$i \n\r";
	if($i==0){
		for(my $j=0; $j<length($secondSeq)+1; $j++){
			# print "j: $j \n\r";
			$newHash{"$i,$j"} = $gap*$j;
			# print "!",$newHash{"$i,$j"},"!\n";
		}
	}else{
		$newHash{"$i,0"} = $gap*$i;
		# print "~",$newHash{"$i,0"},"~\n";
	}
}

# printHash(\%newHash,length($firstSeq)+1,length($secondSeq)+1);
print "Reading from input my 2 sequences are\n";
print "i.$firstSeq\n";
print "j.$secondSeq\n";
print "Scoring\n";
print "Match: $match\n";
print "Mismatch: $mismatch\n";
print "Gap: $gap\n";

my ($newHash2,$posHash) = &fillPos(\%newHash,$firstSeq,$secondSeq,$gap,$match,$mismatch);
#fix the references 
my %newHash2 = %{$newHash2};
my %posHash = %{$posHash}; 
print "--Matrix with scores--\n";
printHash(\%newHash2,length($firstSeq)+1, length($secondSeq)+1);
print "--BackTrack Matrix--\n";
printHash(\%posHash,length($firstSeq)+1, length($secondSeq)+1);
print "\n";
# print "$firstSeq ",length($firstSeq)," LOL $firstSeq ",length($secondSeq),"\n";
my ($path,$score) = backTrace(\%newHash2,\%posHash,length($firstSeq),length($secondSeq));
print "Path for best alignment $path with score $score\n";
my @stop = reverse(split(":", $path));

my ($starti,$startj) = split(",",$stop[0]);
my %hash;
$hash{"$firstSeq"}="";
$hash{"$secondSeq"}="";
my @temp1= split("", $firstSeq);
my @temp2= split("", $secondSeq);


foreach my $s (@stop){
	if($s eq $stop[0]){
		next;
	}


	my ($stopi,$stopj) = split(",",$s);
	# print "$stopi,$stopj\n";
	if ($startj != $stopj and $starti != $stopi){ #diagonal
		$hash{"$firstSeq"} .= "$temp1[$starti] ";
		$hash{"$secondSeq"} .= "$temp2[$startj] ";
	}elsif($startj != $stopj and $starti == $stopi){ #right
		$hash{"$firstSeq"} .= "- ";
		$hash{"$secondSeq"} .= "$temp2[$startj] ";
	}elsif($startj == $stopj and $starti != $stopi){ #down
		$hash{"$firstSeq"} .= "$temp1[$starti] ";
		$hash{"$secondSeq"} .= "- ";
	}
	$starti = $stopi;
	$startj = $stopj;
}
print "BEST ALIGNMENT \n";
print $hash{"$firstSeq"},"\n";
print $hash{"$secondSeq"},"\n";
# printHash(\%posHash,length($firstSeq)+1, length($secondSeq)+1);

sub backTrace{
	my %hash = %{$_[0]};
	my %pos = %{$_[1]};
	my $i = $_[2];
	my $j = $_[3];

	my $score = $hash{"$i,$j"};
	my $index0 = $pos{"$i,$j"};

	my $path = "$i,$j";

	if($i eq "0" or $j eq "0"){
		return $path,$score;
	}else{
		# sleep(1);
		my @paths = split(":",$index0);
		my $flag=0;
		my $posScore = undef;
		my @posPath;
		foreach my $m (@paths) {
			
			my ($posI,$posJ) = split(",", $m);
			my ($resulPath,$score) = backTrace(\%newHash2,\%posHash,$posI,$posJ);
			if(defined $posScore){
				if($score>$posScore){
					$posScore=$score;
					@posPath = ($resulPath);
				}elsif($score==$posScore){
					push @posPath,($resulPath);
				}
			}else{
				$posScore=$score;
				@posPath = ($resulPath);
			}
		}
		$score += $posScore;
		if($path eq "7,6"){
			# print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		}
		$path .=":";
		my $dummy= (scalar(@posPath)>1)?"&":"";
		foreach my $mpath (@posPath){
			$path .= ($dummy.$mpath.$dummy);
		}
	
		# print "$path AND $score\n";

		return $path,$score;
				
	}
}

sub fillPos{

	my %hash = %{$_[0]};
	my %pos;
	my @firstSeq = split("", $_[1]);
	my @secondSeq = split("", $_[2]);
	my $gap = $_[3];
	my $match = $_[4];
	my $mismatch = $_[5];
	
	for(my $i=1; $i<length($firstSeq)+1;$i++){
		for(my $j=1; $j<length($secondSeq)+1;$j++){
			# print "($i:$j) $firstSeq[$i-1] eq $secondSeq[$j-1]\n";
			my $index1 = ($i-1).",".($j-1); # plus gap (up left)
			my $penalty = ($firstSeq[$i-1] eq $secondSeq[$j-1])?$match:$mismatch;
			# print "1.checking $index1 ~ $hash{$index1} + $penalty\n";

			my $value = ($hash{$index1} + $penalty);
			$pos{"$i,$j"} = $index1;

			my $index2 = ($i).",".($j-1);	#plus match/mismatch (up)
			# my @temp = split("|",$value);
			if( $value < ($hash{$index2} + $gap)) {
				$value = ($hash{$index2}+$gap);
				$pos{"$i,$j"} = $index2;
			}elsif($value == ($hash{$index2} + $gap)){
				$pos{"$i,$j"} .= (":".$index2);
			}
			# (split("|",$value)[0]<$hash{$index2} + $penalty)?(($hash{$index2}+$penalty)."|".$index2):$value;
			# print "2.checking $index2 ~ $hash{$index2} + $gap\n";

			my $index3 = ($i-1).",".($j);	#plus match/mismatch (left)
			# @temp = split("|",$value);

			if( $value < ($hash{$index3} + $gap)){
				$value = ($hash{$index3}+$gap);
				$pos{"$i,$j"} = $index3;
			}elsif($value == ($hash{$index3} + $gap)){
				$pos{"$i,$j"} .= (":".$index3);
			}
			# $value = (split("|",$value)[0]<$hash{$index3} + $penalty)?(($hash{$index3} + $penalty)."|".$index3):$value;
			$hash{"$i,$j"} = $value;
			# print "3.checking $index3 ~ $hash{$index3} + $gap\n";

			# print "$i,$j was filled by ",$pos{"$i,$j"},"\n";
		}
	}
	return (\%hash, \%pos);
}

sub printHash{

	my %hash = %{$_[0]};
	my $sizeI = $_[1];
	my $sizeJ = $_[2];

	# print $sizeJ,"\n";
	# print $sizeI,"\n";

	for(my $i=0;$i<$sizeI;$i++){
		for(my$j=0;$j<$sizeJ;$j++){
			if(exists $hash{"$i,$j"}){
				my $plus = (ref($hash{"$i,$j"}) eq "LVALUE" and $hash{"$i,$j"}<0)?"&&  ":"&&   ";
				
				printf "|$i,$j| %-7s ",$hash{"$i,$j"}	;
			}else{
				printf "? ";
			}
		}
		print "\n";
	}
}