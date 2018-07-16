use strict;
use warnings;

my $usage="$0 <diff> <annf>
<diff> this file is the enrichment results, columns are seperated by table,
       and the last column is the gene ids belong to the class
<annf> the annotation file with the format: the first column is the gene id,
       and the remain columns are the annotation information. the first line
       should be the head information.
NOTE: the output of this program is the diff.anno. with one line in diff and
      the following lines are the gene annotation of the genes to this class.
\n";

die $usage if @ARGV<2;

my($diff,$annf)=@ARGV;

open I,$annf;
my(@info,%ann);
my $head=<I>;chomp $head;
while(<I>){
  chomp;
  @info=split /\t/;
  $ann{$info[0]}=join("\t",@info[1..$#info]);
}
close I;

open I,$diff;
open O,">",$diff.".anno.xls";
<I>;
print O $head,"\n";
while(<I>){
  chomp;
  @info=split /\t/;
  print O ">$info[0]\n";
  map{print O "$_\t$ann{$_}\n"}split /[\s,]/,$info[$#info];
  print O "\n";
}
close O;
close I;

