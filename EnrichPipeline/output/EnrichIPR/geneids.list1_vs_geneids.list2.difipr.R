	source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/CompareTwolist.IPR.R")
  source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/Fisher.Chi.test.R");
  source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/sort.data.frame.R");
	SupplyID1 <- scan("input/geneids.list1",what="character")
	SupplyID2 <- scan("input/geneids.list2",what="character")
	CompareTwolist.IPR(supplyID1=SupplyID1,supplyID2=SupplyID2,ipr2geneF="input/ipr2gene.txt",samp1="geneids.list1",samp2="geneids.list2",
		p.adjust.method="fdr",test.method="FisherChiSquare",
		filt="adjp",pc=0.05,enrichFile="output/EnrichIPR//geneids.list1_vs_geneids.list2.difipr")
	q(save="no")
