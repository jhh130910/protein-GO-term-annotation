use strict;
use warnings;

my $usage="$0 <treef>

the output of this program is two file with suffix
.ancest and .offsprings
\n";

die $usage if @ARGV<1;

my($treef)=@ARGV;

my(@info,%ipr,$num1,$num2,@parents,$sym,$id1,$id2,%ipr_ances);

open I,$treef;
push @parents,"a";
while(<I>){
	chomp;
	($sym,$id1)=$_=~/^(.*?)(IPR\d+):/g;
	if(length($sym)==0){
		@parents=();
		push @parents,$id1;
	}else{
		my $len=length($sym)/2;
		if($len <= $num1){
			@parents=@parents[0..($len-1)]
		}
		foreach $id2 (@parents){
			$ipr{$id2}{$id1}++;
			$ipr_ances{$id1}{$id2}++;
		}
		push @parents,$id1;
	}
	$num1=length($sym)/2;
}
close I;

open O,">","$treef.offsprings";
print O "IPR\tOffsprings\n";
foreach $id1 (sort keys %ipr){
	print O "$id1\t",join(",",sort keys %{$ipr{$id1}}),"\n";
}
close O;

open O,">","$treef.ancest";
print O "IPR\tAncestors\n";
foreach $id1 (sort keys %ipr_ances){
	print O "$id1\t",join(",",sort keys %{$ipr_ances{$id1}}),"\n";
}
close O;



