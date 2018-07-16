CompareTwolist.IPR <- function(ipr2geneF,supplyID1,supplyID2,samp1,samp2,
                      p.adjust.methods=c("holm","hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                      enrichFile="enrichOutfile.difipr",test.method=c("FisherChiSquare","HyperGeometric"),
                      filt=c("p","adjp"),pc=0.05,Unanno=c("Y","N"))
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
    if(missing(filt)) filt <- "adjp"
    if(missing(test.method)) test.method="FisherChiSquare"

		supplyID1 <- intersect(supplyID1,names(gene2ipr))
		supplyID2 <- intersect(supplyID2,names(gene2ipr))
		annoDiffNum1 <- length(supplyID1)
		annoDiffNum2 <- length(supplyID2)
    iprNum <- length(ipr2gene)
	diffIPR.tbl <-
      data.frame(
                 "IPRID"=names(ipr2gene),
                 "IPRTitle"=as.character(ipr.title[names(ipr2gene)]),
                 "Pvalue"=rep(0,iprNum),"AdjustedPv"=rep(0,iprNum),
                 "x"=rep(0,iprNum),
                 "y"=rep(0,iprNum),
                 "n"=rep(0,iprNum),
                 "N"=rep(0,iprNum),
                 "WhichOver"=rep("a",iprNum),
                 "GeneIDs"=rep("a",iprNum),
                 "GeneIDs"=rep("a",iprNum),stringsAsFactors=FALSE)
    rownames(diffIPR.tbl) <- names(ipr2gene)
		colnames(diffIPR.tbl) <- c("IPRID","IPRTitle","Pvalue","AdjustedPv",
                  paste(c("x1","x2","n1","n2"),c(samp1,samp2),sep="."),"WhichOver",
                  paste("GeneIDs",c(samp1,samp2),sep="."))
    for(j in rownames(diffIPR.tbl)){
			inter.ids1 <- intersect(supplyID1,ipr2gene[[j]])
			inter.ids2 <- intersect(supplyID2,ipr2gene[[j]])
	  	knum1 <- length(inter.ids1)
			knum2 <- length(inter.ids2)
      if(knum1 <1 & knum2 <1) next
      diffIPR.tbl[j,5:8] <- c(knum1,knum2,annoDiffNum1,annoDiffNum2)
      dif.data <- c(knum1,knum2,annoDiffNum1-knum1,annoDiffNum2-knum2)
      mt <- matrix(dif.data,nr=2)
      test.pv <- switch(test.method,
                        FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                        HyperGeometric=hyper.test(mt)
                        )
      diffIPR.tbl[j,3] <- test.pv
      diffIPR.tbl[j,9] <- ifelse(knum1/annoDiffNum1 > knum2/annoDiffNum2,
                                samp1,samp2)
      diffIPR.tbl[j,10] <- paste(inter.ids1,collapse=" ")
			diffIPR.tbl[j,11] <- paste(inter.ids2,collapse=" ")
    }
    diffIPR.tbl <- diffIPR.tbl[diffIPR.tbl[,5]>0 | diffIPR.tbl[,6]>0,]
	diffIPR.tbl[,4] <- p.adjust(diffIPR.tbl[,3],method=p.adjust.methods)
	diffIPR.tbl <- sort.data.frame(diffIPR.tbl,key="Pvalue")
    write.table(diffIPR.tbl,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.filt <- paste(enrichFile,".filt",sep="")
    if(filt=="adjp")
      diffIPR.filt <- diffIPR.tbl[diffIPR.tbl[,"AdjustedPv"] < pc,]
    else
      diffIPR.filt <- diffIPR.tbl[diffIPR.tbl[,"Pvalue"] < pc,]   
    write.table(diffIPR.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
    return(diffIPR.tbl)
  }
