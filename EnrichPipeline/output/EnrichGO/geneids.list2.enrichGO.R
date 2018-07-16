    load("output/GOdata/GOdata_example.RData")
    source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/EnrichGO.R")
    source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/EnrichGOallLevl.R")
    source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/MergGO.R")
    ids <- scan("input/geneids.list2",what="character")
    if("NullFile" != "NullFile") {
	backid <- scan("NullFile",what="character")
	aa <- EnrichGO(ids,univerID=backid,enrichFile="output/EnrichGO//geneids.list2.difgo3",p.adjust.method="fdr",test.method="FisherChiSquare")
	aa <- EnrichGOallLevl(ids,univerID=backid,enrichFile="output/EnrichGO//geneids.list2.difgoall",p.adjust.method="fdr",test.method="FisherChiSquare")
    }else{
	aa <- EnrichGO(ids,enrichFile="output/EnrichGO//geneids.list2.difgo3",p.adjust.method="fdr",test.method="FisherChiSquare")
        aa <- EnrichGOallLevl(ids,enrichFile="output/EnrichGO//geneids.list2.difgoall",p.adjust.method="fdr",test.method="FisherChiSquare")
    }
    q(save="no")
