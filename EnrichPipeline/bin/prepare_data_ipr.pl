use strict;
use warnings;

my $usage="$0 <ipr.gene>\n";

die $usage if @ARGV<1;

my $file=shift;

open I,$file;
open O,">",$file.".forenrich";
my(@info);
while(<I>){
	chomp;
	@info=split /\t/;
	print O join("\t",@info[0,2,1,3..$#info]),"\n";
}
close I;
close O;

