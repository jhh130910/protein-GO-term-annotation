use strict;
use warnings;

my $usage="$0 <RData> <gofile> <outgene>\n";


die $usage if @ARGV<3;
my($RData,$gofile,$outgene)=@ARGV;

GetGeneByGO($RData,$gofile,$outgene);


sub GetGeneByGO{
	my($RData,$gofile,$outgene)=@_;

	my $rscp=<<"RSCP";
	load("$RData")
	goids <- scan("$gofile",what="character")
	go2allgene <- 
	sapply(goids,function(x){
			paste(names(gene2go)[sapply(gene2go,function(x2){
					length(intersect(x2,all.off[[x]]))>0
					})],collapse=",")
			})
		go.onto <- sapply(goids,function(x) goOnto[[x]][2])
		out <- data.frame(GOID=goids,GO_Term=go.onto,GeneIDs=go2allgene)
		write.table(out,file="$outgene",sep="\\t",row.names=FALSE,quote=FALSE)

RSCP

	my $out_rs=$outgene.".R";
	my $out_rsout=$outgene.".Rout";
	open O,">",$out_rs;
	print O $rscp;
	close O;
	`R CMD BATCH $out_rs $out_rsout`;
	`rm $out_rs;rm $out_rsout`;
}


