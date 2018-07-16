  source("/share/project002/yangpch/bin/EnrichPipeline/bin/EnrichIPR.R");
  source("/share/project002/yangpch/bin/EnrichPipeline/bin/Fisher.Chi.test.R");
  source("/share/project002/yangpch/bin/EnrichPipeline/bin/sort.data.frame.R");
  source("/share/project002/yangpch/bin/EnrichPipeline/bin/hyper.test.R")
  ids <- scan("input/geneids.list2",what="character",sep="\n");
  if("NullFile" != "NullFile"){
      backid <- scan("NullFile",what="character");
      aa <- EnrichIPR("input/ipr2gene.txt",supplyID=ids,univerID=backid,enrichFile="output/EnrichIPR/geneids.list2.difipr",
		      p.adjust.methods="fdr",test.method="FisherChiSquare");
      }else{
	  aa <- EnrichIPR("input/ipr2gene.txt",supplyID=ids,enrichFile="output/EnrichIPR/geneids.list2.difipr",
			  p.adjust.methods="fdr",test.method="FisherChiSquare");
      }
  q(save="no")

