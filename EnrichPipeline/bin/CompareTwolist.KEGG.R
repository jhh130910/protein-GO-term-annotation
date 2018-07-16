CompareTwolist.KEGG <- function(supplyID1,supplyID2,map.title.file,map2gene.file,samp1,samp2,
                       p.adjust.methods=c("holm","hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                       enrichFile="enrichOutfile.difkegg",test.method=c("FisherChiSquare","HyperGeometric"),
                       filt=c("p","adjp"),pc=0.05,Unanno=c("Y","N"))
{
  ##prepare the used data, pathway id to title relation
  aa <- scan(map.title.file,sep="\n",what="character")
  kegg.pthid2ti <- as.character(sapply(aa,function(x) strsplit(x,"\t")[[1]][2]))
  names(kegg.pthid2ti) <- as.character(sapply(aa,function(x) strsplit(x,"\t")[[1]][1]))

  ##kegg pathway to genes and gene to pathways
  kegg.ann <- scan(map2gene.file,what="character",sep="\n")
  kegg.ann.nam <- 
    as.character(sapply(kegg.ann,function(x) strsplit(x,"\t")[[1]][1]))
  map2gene <- sapply(kegg.ann,function(x) strsplit(x,"\t")[[1]][-1])
  names(map2gene) <- kegg.ann.nam
  gene2map <- tapply(rep(kegg.ann.nam,sapply(map2gene,length)),
                     as.character(unlist(map2gene)),unique)
  
  if(missing(p.adjust.methods)) p.adjust.methods="fdr"
  if(missing(filt)) filt <- "adjp"
  if(missing(test.method)) test.method="FisherChiSquare"
	
	supplyID1 <- intersect(supplyID1,names(gene2map))
	supplyID2 <- intersect(supplyID2,names(gene2map))
	annoDiffNum1 <- length(supplyID1)
	annoDiffNum2 <- length(supplyID2)	

  mapNum <- length(map2gene)
  diffKegg.tbl <-
    data.frame("MapID"=names(map2gene),
               "MapTitle"=as.character(kegg.pthid2ti[names(map2gene)]),
               "Pvalue"=rep(0,mapNum),"AdjustedPv"=rep(0,mapNum),
               "x"=rep(0,mapNum),
               "y"=rep(0,mapNum),
               "n"=rep(0,mapNum),
               "N"=rep(0,mapNum),
               "WhichOver"=rep("a",mapNum),
               "GeneIDs1"=rep("a",mapNum),
               "GeneIDs2"=rep("a",mapNum),stringsAsFactors=FALSE)
  rownames(diffKegg.tbl) <- names(map2gene)
	colnames(diffKegg.tbl) <- c("MapID","MapTitle","Pvalue","AdjustedPv",
                        paste(c("x1","x2","n1","n2"),c(samp1,samp2),sep="."),"WhichOver",
                        paste("GeneIDs",c(samp1,samp2),sep="."))
  for(j in rownames(diffKegg.tbl)){
		inter.ids1 <- intersect(supplyID1,map2gene[[j]])
		inter.ids2 <- intersect(supplyID2,map2gene[[j]])
		knum1 <- length(inter.ids1)
		knum2 <- length(inter.ids2)
    if(knum1 <1 & knum2 <1) next
    diffKegg.tbl[j,5:8] <- c(knum1,knum2,annoDiffNum1,annoDiffNum2)
    dif.data <- c(knum1,knum2,annoDiffNum1-knum1,annoDiffNum2-knum2)
    mt <- matrix(dif.data,nr=2)
    test.pv <- switch(test.method,
                      FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                      HyperGeometric=hyper.test(mt)
                      )
    diffKegg.tbl[j,3] <- test.pv
    diffKegg.tbl[j,9] <- ifelse(knum1/annoDiffNum1 > knum2/annoDiffNum2,
                                samp1,samp2)
    diffKegg.tbl[j,10] <- paste(inter.ids1,collapse=" ")
		diffKegg.tbl[j,11] <- paste(inter.ids2,collapse=" ")
  }
  diffKegg.tbl <- diffKegg.tbl[diffKegg.tbl[,5]>0 & diffKegg.tbl[,6]>0,]
  diffKegg.tbl[,4] <- p.adjust(diffKegg.tbl[,3],method=p.adjust.methods)
  diffKegg.tbl <- sort.data.frame(diffKegg.tbl,key="Pvalue")
	diffKegg.tbl[,1] <- paste("map",diffKegg.tbl[,1],sep="")
  write.table(diffKegg.tbl,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
  enrichFile.filt <- paste(enrichFile,".filt",sep="")
  if(filt=="adjp")
    diffKegg.filt <- diffKegg.tbl[diffKegg.tbl[,"AdjustedPv"] < pc,]
  else
    diffKegg.filt <- diffKegg.tbl[diffKegg.tbl[,"Pvalue"] < pc,]  
  write.table(diffKegg.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
  return(diffKegg.tbl)
}
