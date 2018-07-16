use strict;
use warnings;

my $usage="$0 <map.gene>\n";

die $usage if @ARGV<1;

my $file=shift;

open I,$file;
my(@info);
open O,">",$file.".forenrich";
while(<I>){
	chomp;
	@info=split /\s+/;
	$info[0] =~ s/map//;
	print O $info[0];
	foreach my $id (@info[2..$#info]){
	 my @info2=split /,/,$id;
	 print O "\t$info2[0]";
	}
	print O "\n";
}
close I;
close O;

