#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;

my $usage="$0 [options]
--cls     <str> [GO|KEGG|IPR] class to do enrichment
--GOdata  <str> GOdata.RData file that store the objects of GO information
--wegoF   <str> wego format file store the go information for every gene
--MetaGO  <str> MetaGO.RData file store the preprocessed origininal GO information
--GOdata_outF [GOdata.RData] the output file of the first step
--goLev   <str> [3|all|1] GO enrichment at level 3 or all level default[1] both 3 and all
--mapGeneF <str> kegg.map.gene file, in the format like wego\n
--ipr2geneF <str> file IPR 2 Gene\n
--Rpath   <str> the absolute path of R program, default [/share/raid1/genome/bin/R]
--supplyF <str> supplied gene ids file
--backF   <str> background gene ids file.
--Unanno  <str> [Y|N] Default [N]. [Y]es or [N]o include the unannotated genes when counting the
          gene number. If [Y],  the background gene ids, including the unannotated can be supplied
				  with --backF parameter
--P_Adjust_Method <holm|hochberg|hommel|bonferroni|BH|BY|fdr|none> one of them. Default [fdr]
--filt    <str> method to filt the enrichment results <p|adjp>, filter the results using P value or
              the adjusted p value, default [p]
--pc      <fload> the upper value of the p/adjusted p. default [0.05]
--TestMethod <str> the hypothesis test method (FisherChiSquare, HyperGeometric),default [FisherChiSquare]
--outDir directory of output files. The default output file of GO|KEGG|IPR enrichment analysis
\tis supplyF.difgo|difkegg|difipr
\nExample:
$0 --cls GO --MetaGO MetaGO.RData --wegoF gene.wego --outDir ./
$0 --cls GO --GOdata GOdata.RData --supplyF ids1.txt
$0 --cls GO --GOdata GOdata.RData --supplyF ids1.txt --backF total_ids.txt
$0 --cls GO --GOdata GOdata.RData --supplyF ids1.txt --P_Adjust_Method bonferroni
$0 --cls IPR --ipr2geneF ipr2gene.txt --supplyF ids1.txt
$0 --cls KEGG --mapGene hsal.kegg.map.gene --supplyF ids1.txt
\n\n Information:  \n
Author: Pengcheng Yang yangpengcheng\@genomics.org.cn
The detailed introduction and versioin can be found in readme.txt file\n

Reference:
(if you used this program, any of these two papers should be cited)

1. Chen S, Yang P, Jiang F, Wei Y, Ma Z, Kang L: De Novo Analysis of Transcriptome
Dynamics in the Migratory Locust during the Development of Phase Traits. PLoS One 2010,
5(12):e15633.
2. Huang da, W., B.T. Sherman, and R.A. Lempicki. 2009. Bioinformatics
enrichment tools: paths toward the comprehensive functional analysis of large
gene lists. Nucleic Acids Res 37: 1-13.
3. Beissbarth, T. and T.P. Speed. 2004. GOstat: find statistically
overrepresented Gene Ontologies within a group of genes. Bioinformatics 20:
1464-1465.

Version 1.01 2011.12.19
add --Unanno parameter
\n";



my($cls,$GOdata,$wegoF,$MetaGO,$supplyF,$backF,$outDir,$goLev,$mapGeneF,$ipr2geneF,$GOdata_outF,
   $p_adjust_method,$TestMethod,$Rpath,$filt,$pc);
my($Unanno);
GetOptions("cls:s"=>\$cls,
		   "Rpath:s"=>\$Rpath,
		   "GOdata:s"=>\$GOdata,
		   "wegoF:s"=>\$wegoF,
		   "GOdata_outF:s"=>\$GOdata_outF,
		   "MetaGO:s"=>\$MetaGO,
		   "supplyF:s"=>\$supplyF,
		   "backF:s"=>\$backF,
			 "Unanno:s"=>\$Unanno,
		   "P_Adjust_Method:s"=>\$p_adjust_method,
		   "TestMethod:s"=>\$TestMethod,
		   "ipr2geneF:s"=>\$ipr2geneF,
		   "outDir:s"=>\$outDir,
		   "goLev:s"=>\$goLev,
		   "mapGeneF:s"=>\$mapGeneF,
		   "filt:s"=>\$filt,
		   "pc:f"=>\$pc
		   );
die $usage if !defined $cls;

$outDir ||= "./";
`mkdir $outDir` if ! -e $outDir;
$p_adjust_method ||= "fdr";
$backF ||= "NullFile";
$GOdata_outF ||= "GOdata.RData";
$TestMethod ||= "FisherChiSquare";
$filt ||= "adjp";
$pc ||= 0.05;
$Unanno ||= "N";

`export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/panfs/ANIMAL/GROUP/group001/yangpch/soft/lib/:/opt/blc/atlas-3.8.3/lib/`;
$Rpath ||= "/share/raid1/genome/bin/R";


################
#GO
if($cls eq "GO"){
    $goLev ||= 1;
	my $idfile=basename($GOdata_outF);
    if(defined $wegoF && defined $MetaGO){
	  my $GOdataRF=$outDir."/".$idfile;$GOdataRF=~s/Data$//;
	  $GOdata_outF = $outDir."/".$idfile;
	  GOdata($wegoF,$GOdataRF,$MetaGO,$outDir,$GOdata_outF);
	  my $GOdataRFout = $GOdataRF."out";
	  `$Rpath CMD BATCH $GOdataRF $GOdataRFout`;
    }
    if(defined $GOdata && defined $supplyF){
	  enrichGO($GOdata,$supplyF,$backF,$goLev,$outDir,$p_adjust_method,$TestMethod,$filt,$pc,$Unanno);
    }
}

################
#KEGG
if($cls eq "KEGG"){
    my $mapTitleF="$Bin/../data/map_title.tab";
    enrichKEGG($mapTitleF,$mapGeneF,$supplyF,$backF,$outDir,$p_adjust_method,$Bin,$TestMethod,$filt,$pc,$Unanno);
}

################
#IPR
if($cls eq "IPR"){
    enrichIPR($ipr2geneF,$supplyF,$backF,$outDir,$p_adjust_method,$Bin,$TestMethod,$filt,$pc,$Unanno);
}


sub enrichIPR{
  my($ipr2geneF,$supplyF,$backF,$outDir,$p_adjust_method,$Bin,$TestMethod,$filt,$pc,$Unanno)=@_;
  my $enrichFile=basename($supplyF);
  my $difiprF = $outDir."/$enrichFile.difipr";
  my $rscrip;
  $rscrip=<<"RSCRIPT";
  source("$Bin/EnrichIPR.R");
  source("$Bin/Fisher.Chi.test.R");
  source("$Bin/sort.data.frame.R");
  source("$Bin/hyper.test.R")
  ids <- scan("$supplyF",what="character",sep="\\n");
  if("$backF" != "NullFile"){
      backid <- scan("$backF",what="character");
      aa <- EnrichIPR("$ipr2geneF",supplyID=ids,univerID=backid,enrichFile="$difiprF",
		      p.adjust.methods="$p_adjust_method",test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno");
      }else{
	  aa <- EnrichIPR("$ipr2geneF",supplyID=ids,enrichFile="$difiprF",
			  p.adjust.methods="$p_adjust_method",test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno");
      }
  q(save="no")
RSCRIPT

  my $enrichRF=$outDir."/".$enrichFile.".enrichIPR.R";
  open R,">",$enrichRF;
  print R $rscrip,"\n"; 
  close R;
  my $enrichRFout = $enrichRF."out";
  `$Rpath CMD BATCH $enrichRF $enrichRFout`;
}

sub enrichKEGG{
  my($mapTitleF,$mapGeneF,$supplyF,$backF,$outDir,$p_adjust_method,$Bin,$TestMethod,$filt,$pc,$Unanno)=@_;
  my $enrichFile=basename($supplyF);
  my $difkeggF=$outDir."/".$enrichFile.".difkegg";
  my $rscrip=<<"RSCRIPT";
  source("$Bin/EnrichKEGG.R");
  source("$Bin/Fisher.Chi.test.R");
  source("$Bin/sort.data.frame.R");
    ids <- scan("$supplyF",what="character",sep="\\n")
    if("$backF" != "NullFile"){
       backid <- scan("$backF",what="character");
       aa <- EnrichKEGG("$mapTitleF","$mapGeneF",supplyID=ids,univerID=backid,enrichFile="$difkeggF",
	                     p.adjust.methods="$p_adjust_method",test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno");
    }else{
       aa <- EnrichKEGG("$mapTitleF","$mapGeneF",supplyID=ids,enrichFile="$difkeggF",
		                 p.adjust.methods="$p_adjust_method",test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno");
    }
    q(save="no")
RSCRIPT

  my $enrichRF=$outDir."/".$enrichFile.".enrichKEGG.R";
  open R,">",$enrichRF;
  print R $rscrip,"\n";
  close R;
  my $enrichRFout = $enrichRF."out";
  `$Rpath CMD BATCH $enrichRF $enrichRFout`;
}


#print out R script file to produce GOdata.RData
sub GOdata{
    my($wegoF,$GOdataRF,$MetaGO,$outDir,$GOdata_outF)=@_;
    my $rscrip=<<"RSCRIPT";
    metago.dir <- "$MetaGO"
	dir.wego <- "$wegoF"
	
##load in the GO Meta data.
	load(metago.dir)
##read in reference GO annotation
	
#########
##GOID to gene products at level 3
	go.anno <- unique(scan(dir.wego,what="character",sep="\\n"))
	go.anno.nam <- as.character(sapply(go.anno,function(x) strsplit(x,"\\t")[[1]][1]))
	go.anno.gos <- sapply(go.anno,function(x) strsplit(x,"\\t")[[1]][-1])
	gene2go <- go.anno.gos
	names(gene2go) <- go.anno.nam
##filter the gos without annotation in goOnto
	gene2go.filt <- sapply(gene2go,function(x) intersect(x,names(goOnto)))
	go.loc <- which(sapply(gene2go.filt,length)<1)
	if (length(go.loc)>=1) gene2go <- gene2go.filt[-go.loc] else gene2go <- gene2go.filt
	
	go2gene <- vector("list",length=length(c(mf.level3,bp.level3,cc.level3)))
	names(go2gene) <- c(mf.level3,bp.level3,cc.level3)
	aces.tt <- c(cc.aces,bp.aces,mf.aces)
	for(i in 1:length(gene2go)){
	    a <- gene2go[[i]]
		a.aces <- unique(as.character(unlist(aces.tt[a])))[-1]
		a.aces <- c(a.aces,a)
		gointer <- unique(intersect(a.aces,names(go2gene)))
		if(length(gointer)==0) next
		for(j in gointer){ go2gene[[j]] <- c(go2gene[[j]],names(gene2go)[i])}
#the next line popup error message that:number of items to replace is not a multiple of replacement length
	}
    go2gene <- go2gene[sapply(go2gene,length)>0]
	
##########
##reference
	genefreq.mf.ref <- GOfreq(genefreq.mf,gene2go,"MF",mf.aces,goOnto)
	genefreq.cc.ref <- GOfreq(genefreq.cc,gene2go,"CC",cc.aces,goOnto)
	genefreq.bp.ref <- GOfreq(genefreq.bp,gene2go,"BP",bp.aces,goOnto)
	genefreq.mf.ref  <- genefreq.mf.ref[genefreq.mf.ref>0]
	genefreq.cc.ref  <- genefreq.cc.ref[genefreq.cc.ref>0]
	genefreq.bp.ref  <- genefreq.bp.ref[genefreq.bp.ref>0]
	
	save(all.aces,all.child,all.parent,all.off,bp.level3,cc.level3,mf.level3,bp.aces,cc.aces,mf.aces,goOnto,
	     genefreq.mf,genefreq.bp,genefreq.cc,genefreq.mf.ref,genefreq.cc.ref,
	     genefreq.bp.ref,go2gene,sort.data.frame,Fisher.Chi.test,
	     DiffGOs,GOfreq,gene2go,file="$GOdata_outF")
	q(save="no")
	
RSCRIPT
open G,">",$GOdataRF;
    print G $rscrip;
    close G;
}

sub enrichGO{
    my($GOdata,$supplyF,$backF,$goLev,$outDir,$p_adjust_method,$TestMethod,$filt,$pc,$Unanno)=@_;
	my $idfile=basename($supplyF);
    my $enrichRF=$outDir."/".$idfile.".enrichGO.R";
    my $go3F=$outDir."/".$idfile.".difgo3";
    my $goallF=$outDir."/".$idfile.".difgoall";
    my $rscrip=<<"RSCRIPT";
    load("$GOdata")
    source("$Bin/EnrichGO.R")
    source("$Bin/EnrichGOallLevl.R")
    source("$Bin/MergGO.R")
    ids <- scan("$supplyF",what="character")
    if("$backF" != "NullFile") {
	backid <- scan("$backF",what="character")
	aa <- EnrichGO(ids,univerID=backid,enrichFile="$go3F",p.adjust.method="$p_adjust_method",
				test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno")
	aa <- EnrichGOallLevl(ids,univerID=backid,enrichFile="$goallF",p.adjust.method="$p_adjust_method",
				test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno")
    }else{
	aa <- EnrichGO(ids,enrichFile="$go3F",p.adjust.method="$p_adjust_method",test.method="$TestMethod",
				filt="$filt",pc=$pc,Unanno="$Unanno")
        aa <- EnrichGOallLevl(ids,enrichFile="$goallF",p.adjust.method="$p_adjust_method",
				test.method="$TestMethod",filt="$filt",pc=$pc,Unanno="$Unanno")
    }
    q(save="no")
RSCRIPT
    open E,">",$enrichRF;
    print E $rscrip;
    close E;
    my $enrichRFout = $enrichRF."out";
    `$Rpath CMD BATCH $enrichRF $enrichRFout`;
}

