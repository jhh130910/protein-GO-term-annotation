    load("output/GOdata_example.RData")
    source("/share/project002/yangpch/bin/EnrichPipeline/bin/EnrichGO.R")
    source("/share/project002/yangpch/bin/EnrichPipeline/bin/EnrichGOallLevl.R")
    ids <- scan("input/geneids.list1",what="character")
    if("NullFile" != "NullFile") {
	backid <- scan("NullFile",what="character")
	aa <- EnrichGO(ids,univerID=backid,enrichFile="output//geneids.list1.difgo3",p.adjust.method="fdr",test.method="FisherChiSquare")
	aa <- EnrichGOallLevl(ids,univerID=backid,enrichFile="output//geneids.list1.difgoall",p.adjust.method="fdr",test.method="FisherChiSquare")
    }else{
	aa <- EnrichGO(ids,enrichFile="output//geneids.list1.difgo3",p.adjust.method="fdr",test.method="FisherChiSquare")
        aa <- EnrichGOallLevl(ids,enrichFile="output//geneids.list1.difgoall",p.adjust.method="fdr",test.method="FisherChiSquare")
    }
    q(save="no")
