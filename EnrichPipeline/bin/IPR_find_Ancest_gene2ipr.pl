use strict;
use warnings;
use FindBin qw($Bin);

my $usage="$0 <gene2iprf>  <out_gene2iprf>

this script is used to generate the gene 2 ipr information that include all the
ancestor IPRs
\n";

die $usage if @ARGV<1;

my($gene2iprf,$out_gene2iprf)=@ARGV;

my $ipr_entryf="$Bin/../data/entry.list";
my $ipr_ancesf="$Bin/../IPR_ParentChildTreeFile.txt.ancest";

my(@info,%iprEntry,%iprAncest,@info2,$ipr1,%gene2ipr);
open I,$ipr_entryf;

while(<I>){
	chomp;
	next if ! /^IPR/;
	@info=split /\s+/;
	$iprEntry{$info[0]}=join(" ",@info[1..$#info]);
}
close I;

open I,$ipr_ancesf;
while(<I>){
	chomp;
	next if ! /^IPR/;
	@info=split /\s+/;
	@info2=split /,/,$info[1];
	map{$iprAncest{$info[0]}{$_}++}@info2;
}
close I;

open I,$gene2iprf;
my(@iprs);
while(<I>){
	chomp;
	@info=split /\t/;
	(@iprs)=$_=~/(IPR\d+)/g;
	map{$gene2ipr{$info[0]}{$_}++}@iprs;
	foreach $ipr1 (@iprs){
		if(exists $iprAncest{$ipr1}){
			foreach my $ipr2 (keys %{$iprAncest{$ipr1}}){
				$gene2ipr{$info[0]}{$ipr2}++;
			}
		}
	}
}
close I;

open O,">",$out_gene2iprf;
foreach my $gene (sort keys %gene2ipr){
	print O "$gene\t",join(",",sort keys %{$gene2ipr{$gene}}),"\n";
}
close O;


