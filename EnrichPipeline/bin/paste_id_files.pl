use strict;
use warnings;
use File::Basename qw(basename dirname);

my $usage="$0 <outfile> <idfile(s)>
<outfile> output file with one column store ids from one file
<idfile(s)> one or more id files, support regular expression
";

die $usage if @ARGV<2;

my @files=@ARGV;
my $outfile=shift @files;
my(%file_ids,$maxnum,@ids);
$maxnum=0;
open O,">",$outfile;
while(my $file=shift @files){
  my $file_base=basename($file);
  open I,$file;
  @ids=<I>;chomp @ids;
  $maxnum = @ids<$maxnum?$maxnum:@ids;
  map{chomp;$file_ids{$file_base}{$_}=$ids[$_]}(0..$#ids);
  close I;
}

@ids=sort keys %file_ids;
print O join("\t",@ids),"\n";
for(my $i=0;$i<$maxnum;$i++){
  my $out;
  foreach my $id (@ids){
    if(exists $file_ids{$id}{$i}){
      $out .= $file_ids{$id}{$i}."\t";
    }else{
      $out .= "\t";
    }
  }
  $out =~ s/\s$//;
  print O "$out\n";
}
close O;
