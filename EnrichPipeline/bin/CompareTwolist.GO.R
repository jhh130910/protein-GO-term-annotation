CompareTwolist.GO <-
    function(supplyID1,supplyID2,samp1,samp2,
           p.adjust.method=c("holm","hochberg","hommel","bonferroni","BH"
           ,"BY","fdr","none"),enrichFile="enrichOut.difgoall",
             test.method=c("FisherChiSquare","HyperGeometric"),filt=c("p","adjp"),pc=0.05
             ){
        opt.str <- options("stringsAsFactors")
        options(stringsAsFactors=FALSE)

        supplyID1 <- intersect(supplyID1,names(gene2go))
        supplyID2 <- intersect(supplyID2,names(gene2go))
        if(any(c(length(supplyID1),length(supplyID2))==0)) stop("Supplied IDs are not annotated!")

        go2gene1.bottom <- tapply(rep(names(gene2go[supplyID1]),sapply(gene2go[supplyID1],length)),
                                  as.character(unlist(gene2go[supplyID1])),unique)
        go2gene2.bottom <- tapply(rep(names(gene2go[supplyID2]),sapply(gene2go[supplyID2],length)),
                              as.character(unlist(gene2go[supplyID2])),unique)

        supply1.goid.self <- unique(as.character(unlist(gene2go[supplyID1])))
        supply2.goid.self <- unique(as.character(unlist(gene2go[supplyID2])))
        supply1.goid.aces <- unique(as.character(unlist(all.aces[supply1.goid.self])))
        supply1.goid <- unique(c(supply1.goid.self,supply1.goid.aces))
        supply1.goid <- supply1.goid[-grep("all",supply1.goid)]

        supply2.goid.aces <- unique(as.character(unlist(all.aces[supply2.goid.self])))
        supply2.goid <- unique(c(supply2.goid.self,supply2.goid.aces))
        supply2.goid <- supply2.goid[-grep("all",supply2.goid)]

        supply.goid <- unique(c(supply1.goid,supply2.goid))
    rowN <- length(supply.goid)

    res.df <-
      data.frame("GO_ID"=supply.goid,
                 "GO_Term"=as.character(sapply(goOnto[supply.goid],function(x) x[[2]])),
                 "GO_Class"=as.character(sapply(goOnto[supply.goid],function(x) x[[1]])),
                 "Pvalue"=rep(0,rowN),
                 "AdjustedPv"=rep(0,rowN),
                 "x1"=rep(0,rowN),
                 "x2"=rep(0,rowN),
                 "n1"=rep(0,rowN),
                 "n2"=rep(0,rowN),
                 "WhichOver"=rep("",rowN),
                 "GOlevl"=as.character(sapply(goOnto[supply.goid],function(x) x[[3]])),
                 "GeneID1"=rep("",rowN),
								 "GeneID2"=rep("",rowN))
    rownames(res.df) <- supply.goid
		colnames(res.df) <- c("GO_ID","GO_Term","GO_Class","Pvalue","AdjustedPv",
										paste(c("x1","x2","n1","n2"),c(samp1,samp2),sep="."),"WhichOver","GOlevl",
										paste("GeneID",c(samp1,samp2),sep="."))
    if(missing(test.method)) test.method="FisherChiSquare"

		num.tot.supply1 <- length(supplyID1)
		num.tot.supply2 <- length(supplyID2)

    for(i in supply.goid){
			inter1.goid <- intersect(c(all.off[[i]],i),names(go2gene1.bottom))
			id.clst.uni1 <- unique(as.character(unlist(go2gene1.bottom[inter1.goid])))
			num.clst.uni1 <- length(id.clst.uni1)

			inter2.goid <- intersect(c(all.off[[i]],i),names(go2gene2.bottom))
			id.clst.uni2 <- unique(as.character(unlist(go2gene2.bottom[inter2.goid])))
			num.clst.uni2 <- length(id.clst.uni2)

      if(all(c(num.clst.uni1,num.clst.uni2)==0)){
        res.df[i,4:9] <- rep(0,6)
        res.df[i,10:13] <- rep("a",3)
        next
      }

      res.df[i,6:9] <- c(num.clst.uni1,num.clst.uni2,num.tot.supply1,num.tot.supply2)
      dif.data <- c(num.clst.uni1,num.clst.uni2,num.tot.supply1-num.clst.uni1,num.tot.supply2-num.clst.uni2)
      mt <- matrix(dif.data,nr=2)

      test.pv <- switch(test.method,
                        FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                        HyperGeometric=hyper.test(mt)
                        )
      res.df[i,"Pvalue"] <- test.pv
      res.df[i,12] <- paste(id.clst.uni1,collapse=",")
			res.df[i,13] <- paste(id.clst.uni2,collapse=",")
      if(num.clst.uni1/num.tot.supply1 > num.clst.uni2/num.tot.supply2)
        res.df[i,"WhichOver"] <- samp1
      else
        res.df[i,"WhichOver"] <- samp2
    }
    ##filter the rows without ids
    res.df <- res.df[res.df[,6]>0 | res.df[,7]>0,]
    if(missing(p.adjust.method)) p.adjust.method="fdr"
    if(missing(filt)) filt <- "adjp"
    res.df[,"AdjustedPv"] <- p.adjust(res.df[,"Pvalue"],method=p.adjust.method)
    res.df <- sort.data.frame(res.df,key="Pvalue")
    options(opt.str)
    write.table(res.df,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.filt <- paste(enrichFile,".filt",sep="")
    if(filt=="adjp")
      res.df.filt <- res.df[res.df[,"AdjustedPv"] < pc,]
    else
      res.df.filt <- res.df[res.df[,"Pvalue"] < pc,]
    write.table(res.df.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.merg <- paste(enrichFile.filt,".merg",sep="")
    if(dim(res.df.filt)[1] >= 1){
      diffgo.merg <- MergGO.pair(res.df.filt)
      write.table(diffgo.merg,file=enrichFile.merg,sep="\t",quote=FALSE,row.names=FALSE)
    }
    return(res.df)
  }
