EnrichKEGG <- function(map.title.file,map2gene.file,supplyID,univerID,
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
  if(missing(univerID)) univerID <- names(gene2map)
  if(missing(filt)) filt <- "adjp"
	if(missing(Unanno) | identical(Unanno,"N")){
  	univerID <- intersect(univerID,names(gene2map))
 	 ids <- intersect(supplyID,univerID)
 	 totAnnNum <- length(univerID)
  	annoDiffNum <- length(ids)
	}else{
		totAnnNum <- length(univerID)
		annoDiffNum <- length(supplyID)
  	univerID <- intersect(univerID,names(gene2map))
 	 ids <- intersect(supplyID,univerID)
	}
  if(length(univerID) != length(gene2map)){
    map2id <- lapply(map2gene,function(x) intersect(x,univerID))
    map2id <- map2id[sapply(map2id,length)>0]
  }else{
    map2id <- map2gene
  }
  mapNum <- length(map2id)
  diffKegg.tbl <-
    data.frame("MapID"=names(map2id),
               "MapTitle"=as.character(kegg.pthid2ti[names(map2id)]),
               "Pvalue"=rep(0,mapNum),"AdjustedPv"=rep(0,mapNum),
               "x"=rep(0,mapNum),
               "y"=rep(0,mapNum),
               "n"=rep(0,mapNum),
               "N"=rep(0,mapNum),
               "EnrichDirect"=rep("a",mapNum),
               "GeneIDs"=rep("a",mapNum),stringsAsFactors=FALSE)
  rownames(diffKegg.tbl) <- names(map2id)
  for(j in rownames(diffKegg.tbl)){
    knum <- length(intersect(ids,map2id[[j]]))
    if(knum <1) next
    dif.data <- c(knum,length(map2id[[j]]),annoDiffNum,totAnnNum)
    diffKegg.tbl[j,5:8] <- dif.data
    knum.rm <- length(map2id[[j]])-knum
    tot.rm <- totAnnNum-annoDiffNum-knum.rm
    dif.data <- c(knum,knum.rm,annoDiffNum-knum,tot.rm)
    mt <- matrix(dif.data,nr=2)
    if(missing(test.method)) test.method="FisherChiSquare"
    test.pv <- switch(test.method,
                      FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                      HyperGeometric=hyper.test(mt)
                      )
    diffKegg.tbl[j,3] <- test.pv
    diffKegg.tbl[j,9] <- ifelse((dif.data[1]/dif.data[2])>(dif.data[3]/dif.data[4]),
                                "Over","Under")
    diffKegg.tbl[j,10] <- paste(intersect(ids,map2id[[j]]),collapse=" ")
  }
  diffKegg.tbl <- diffKegg.tbl[diffKegg.tbl[,5]>0,]
  diffKegg.tbl[,4] <- p.adjust(diffKegg.tbl[,3],method=p.adjust.methods)
  diffKegg.tbl <- sort.data.frame(diffKegg.tbl,key="Pvalue")
	diffKegg.tbl[,1] <- paste("map",diffKegg.tbl[,1],sep="")
  write.table(diffKegg.tbl,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
  enrichFile.filt <- paste(enrichFile,".filt",sep="")
  if(filt=="adjp")
    diffKegg.filt <- diffKegg.tbl[diffKegg.tbl[,"AdjustedPv"] < pc & diffKegg.tbl[,"EnrichDirect"]=="Over",]
  else
    diffKegg.filt <- diffKegg.tbl[diffKegg.tbl[,"Pvalue"] < pc & diffKegg.tbl[,"EnrichDirect"]=="Over",]  
  write.table(diffKegg.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
  return(diffKegg.tbl)
}
