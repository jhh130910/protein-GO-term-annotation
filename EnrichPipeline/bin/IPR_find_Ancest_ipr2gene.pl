use strict;
use warnings;
use FindBin qw($Bin);

my $usage="$0 <ipr2genef> <out_ipr2genef> <out_gene2iprf>

this script is used to add the ancestory IPRs to the ipr2gene file that used to
do the enrichment analysis
\n";

die $usage if @ARGV<1;

my($ipr2genef,$out_ipr2genef,$out_gene2iprf)=@ARGV;


my $ipr_entryf="$Bin/../data/entry.list";
my $ipr_ancesf="$Bin/../data/ParentChildTreeFile.txt.ancest";

my(@info,%iprEntry,%iprAncest,@info2,%out,$ipr1,%gene2ipr);
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

open I,$ipr2genef;
while(<I>){
	chomp;
	@info=split /\t/;
	$iprEntry{$info[0]}=$info[2] if ! exists $iprEntry{$info[0]};
	map{$out{$info[0]}{$_}++}@info[3..$#info];
	map{$gene2ipr{$_}{$info[0]}++}@info[3..$#info];
	if(exists $iprAncest{$info[0]}){
		foreach $ipr1 (keys %{$iprAncest{$info[0]}}){
			map{$out{$ipr1}{$_}++}@info[3..$#info];
			map{$gene2ipr{$_}{$ipr1}++}@info[3..$#info];

		}
	}
}
close I;

open O,">",$out_ipr2genef;
foreach $ipr1 (sort keys %out){
	my $num=scalar keys %{$out{$ipr1}};
	print O "$ipr1\t$num\t$iprEntry{$ipr1}\t",join("\t",sort keys %{$out{$ipr1}}),"\n";
}
close O;

open O,">",$out_gene2iprf;
foreach my $gene1 (sort keys %gene2ipr){
	print O "$gene1\t",join(",",sort keys %{$gene2ipr{$gene1}}),"\n";
}
close O;



