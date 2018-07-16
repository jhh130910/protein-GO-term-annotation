EnrichIPR <- function(ipr2geneF,supplyID,univerID,##the background gene list
                      p.adjust.methods=c("holm","hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                      enrichFile="enrichOutfile.difipr",test.method=c("FisherChiSquare","HyperGeometric"),
                      filt=c("p","adjp"),pc=0.05,Unanno=c("Y","N"),pair.file)
  {
    ipr <- scan(ipr2geneF,what="character",sep="\n")
    ipr.name <- as.character(sapply(ipr,function(x) strsplit(x,"\t")[[1]][1]))
    ipr.title <- as.character(sapply(ipr,function(x) strsplit(x,"\t")[[1]][3]))
    names(ipr.title) <- ipr.name
    ipr2gene <- sapply(ipr,function(x) strsplit(x,"\t")[[1]][-c(1:3)])
    names(ipr2gene) <- ipr.name
    gene2ipr <- tapply(rep(ipr.name,sapply(ipr2gene,length)),
                     as.character(unlist(ipr2gene)),unique)
    if(missing(p.adjust.methods)) p.adjust.methods="fdr"
    if(missing(univerID)) univerID=names(gene2ipr)
    if(missing(filt)) filt <- "adjp"
	if(missing(Unanno)){
		totAnnNum <- length(univerID)
		annoDiffNum <- length(supplyID)
	}else{
		univerID <- intersect(univerID,names(gene2ipr))
		ids <- intersect(supplyID,univerID)
		totAnnNum <- length(univerID)
		annoDiffNum <- length(ids)
	}
	univerID <- intersect(univerID,names(gene2ipr))
	ids <- intersect(supplyID,univerID)
	if(length(gene2ipr) != length(univerID)){
      ipr2id <- lapply(ipr2gene,function(x) intersect(x,univerID))
      ipr2id <- ipr2id[sapply(ipr2id,length)>0]
    }else{
      ipr2id <- ipr2gene
    }
    iprNum <- length(ipr2id)
	diffIPR.tbl <-
      data.frame(
                 "IPRID"=names(ipr2id),
                 "IPRTitle"=as.character(ipr.title[names(ipr2id)]),
                 "Pvalue"=rep(0,iprNum),"AdjustedPv"=rep(0,iprNum),
                 "x"=rep(0,iprNum),
                 "y"=rep(0,iprNum),
                 "n"=rep(0,iprNum),
                 "N"=rep(0,iprNum),
                 "EnrichDirect"=rep("a",iprNum),
                 "GeneIDs"=rep("a",iprNum),stringsAsFactors=FALSE)
    rownames(diffIPR.tbl) <- names(ipr2id)
    for(j in rownames(diffIPR.tbl)){
	  knum <- length(intersect(ids,ipr2id[[j]]))
      if(knum <1) next
      dif.data <- c(knum,length(ipr2id[[j]]),annoDiffNum,totAnnNum)
      diffIPR.tbl[j,5:8] <- dif.data
      knum.rm <- length(ipr2id[[j]])-knum
      tot.rm <- totAnnNum-annoDiffNum-knum.rm
      dif.data <- c(knum,knum.rm,annoDiffNum-knum,tot.rm)
      mt <- matrix(dif.data,nr=2)
      if(missing(test.method)) test.method="FisherChiSquare"
      test.pv <- switch(test.method,
                        FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                        HyperGeometric=hyper.test(mt)
                        )
      diffIPR.tbl[j,3] <- test.pv
      diffIPR.tbl[j,9] <- ifelse((dif.data[1]/dif.data[2])>(dif.data[3]/dif.data[4]),
                                 "Over","Under")
      diffIPR.tbl[j,10] <- paste(intersect(ids,ipr2id[[j]]),collapse=" ")
    }
    diffIPR.tbl <- diffIPR.tbl[diffIPR.tbl[,5]>0,]
	diffIPR.tbl[,4] <- p.adjust(diffIPR.tbl[,3],method=p.adjust.methods)
	diffIPR.tbl <- sort.data.frame(diffIPR.tbl,key="Pvalue")
    write.table(diffIPR.tbl,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.filt <- paste(enrichFile,".filt",sep="")
    if(filt=="adjp")
      diffIPR.filt <- diffIPR.tbl[diffIPR.tbl[,"AdjustedPv"] < pc & diffIPR.tbl[,"EnrichDirect"]=="Over",]
    else
      diffIPR.filt <- diffIPR.tbl[diffIPR.tbl[,"Pvalue"] < pc & diffIPR.tbl[,"EnrichDirect"]=="Over",]   
    write.table(diffIPR.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
		merge.file <- paste(enrichFile.filt,".merg",sep="")
		if(dim(diffIPR.filt)[1] > 1){
			diffIPR.filt.merg <- MergIPR(diffIPR.filt,pair.file)
			write.table(diffIPR.filt.merg,file=merge.file,quote=FALSE,row.names=FALSE,sep="\t")
		}
    return(diffIPR.tbl)
  }
