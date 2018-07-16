	load("output/GOdata/GOdata_example.RData")
	source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/CompareTwolist.GO.R")
	source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/MergGO.pair.R")
	SupplyID1 <- scan("input/geneids.list1",what="character")
	SupplyID2 <- scan("input/geneids.list2",what="character")
	CompareTwolist.GO(SupplyID1,SupplyID2,samp1="geneids.list1",samp2="geneids.list2",p.adjust.method="fdr",test.method="FisherChiSquare",
		filt="adjp",pc=0.05,enrichFile="output/EnrichGO//geneids.list1_vs_geneids.list2.difgoall")
	q(save="no")
