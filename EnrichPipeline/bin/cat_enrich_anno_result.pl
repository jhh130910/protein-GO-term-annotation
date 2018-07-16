use strict;
use warnings;
use File::Basename qw(basename dirname);

my $usage="$0 <resltDir> <xlsf>
<resltDir> the enrichment analysis result directory, the result with anno.xls
           will be catted into one file
<xlsf>     output result file
";

die $usage if @ARGV<2;

my($resltDir,$xlsf)=@ARGV;

my(@info,$head);
my @files=<$resltDir/*.anno.xls>;
my($file);
open O,">",$xlsf;
my $flg=0;
foreach $file (@files){
  my $file_filt=basename($file);
  $file_filt =~ s/.dif.*?filt*.anno.xls//;
  open I,$file;
  my $head = <I>;
  print O $head if $flg==0;
  print O $file_filt,"\n";
  $flg=1 if $flg==0;
  while(<I>){
    chomp;
    print O $_,"\n";
  }
  close I;
  print O "\n";
}
close O;
