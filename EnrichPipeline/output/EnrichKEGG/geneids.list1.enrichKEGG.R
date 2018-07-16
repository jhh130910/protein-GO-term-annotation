  source("/share/project002/yangpch/bin/EnrichPipeline/bin/EnrichKEGG.R");
  source("/share/project002/yangpch/bin/EnrichPipeline/bin/Fisher.Chi.test.R");
  source("/share/project002/yangpch/bin/EnrichPipeline/bin/sort.data.frame.R");
    ids <- scan("input/geneids.list1",what="character",sep="\n")
    if("NullFile" != "NullFile"){
       backid <- scan("NullFile",what="character");
       aa <- EnrichKEGG("/share/project002/yangpch/bin/EnrichPipeline/bin/../data/map_title.tab","input/KEGG.map.gene.txt",supplyID=ids,univerID=backid,enrichFile="output/EnrichKEGG/geneids.list1.difkegg",
	                     p.adjust.methods="fdr",test.method="FisherChiSquare");
    }else{
       aa <- EnrichKEGG("/share/project002/yangpch/bin/EnrichPipeline/bin/../data/map_title.tab","input/KEGG.map.gene.txt",supplyID=ids,enrichFile="output/EnrichKEGG/geneids.list1.difkegg",
		                 p.adjust.methods="fdr",test.method="FisherChiSquare");
    }
    q(save="no")

