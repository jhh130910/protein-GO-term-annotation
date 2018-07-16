  source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/EnrichIPR.R");
  source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/Fisher.Chi.test.R");
  source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/sort.data.frame.R");
  source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/hyper.test.R")
	source("/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/MergIPR.R")
  ids <- scan("input/geneids.list1",what="character",sep="\n");
  if("NullFile" != "NullFile"){
      backid <- scan("NullFile",what="character");
      aa <- EnrichIPR("input/ipr2gene.txt",supplyID=ids,univerID=backid,enrichFile="output/EnrichIPR/geneids.list1.difipr",pair.file="/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/../data/ParentChildTreeFile.txt",
		      p.adjust.methods="fdr",test.method="FisherChiSquare",filt="p",pc=0.05,Unanno="N");
      }else{
	  aa <- EnrichIPR("input/ipr2gene.txt",supplyID=ids,enrichFile="output/EnrichIPR/geneids.list1.difipr",pair.file="/panfs/ANIMAL/GROUP/group001/yangpch/bin/EnrichPipeline/bin/../data/ParentChildTreeFile.txt",
			  p.adjust.methods="fdr",test.method="FisherChiSquare",filt="p",pc=0.05,Unanno="N");
      }
  q(save="no")

