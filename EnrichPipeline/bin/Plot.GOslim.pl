use strict;
use warnings;

my $usage="$0 <configf> <goslimf> <graph_pdf>
this is a encapsulation of the Plot.GOslim.R\n";

die $usage if @ARGV<3;

my($configf,$goslimf,$graph_pdf)=@ARGV;


