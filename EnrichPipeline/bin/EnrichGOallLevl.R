EnrichGOallLevl <-
  function(supplyID,univerID,
           p.adjust.method=c("holm","hochberg","hommel","bonferroni","BH"
           ,"BY","fdr","none"),enrichFile="enrichOut.difgoall",
           test.method=c("FisherChiSquare","HyperGeometric"),filt=c("p","adjp"),pc=0.05,Unanno=c("Y","N")
           ){
    go2gene.bottom <- tapply(rep(names(gene2go),sapply(gene2go,length)),
                             as.character(unlist(gene2go)),unique)
    opt.str <- options("stringsAsFactors")
    options(stringsAsFactors=FALSE)
    if(missing(univerID)) univerID <- names(gene2go)
		if(missing(Unanno) | identical(Unanno,"N")){
   	 univerID <- intersect(univerID,names(gene2go))
   	 id.ann.supply <- intersect(supplyID,univerID)
    	num.tot.supply <- length(id.ann.supply)
    	num.tot.uni <- length(univerID)
		}else{
			num.tot.supply <- length(supplyID)
			num.tot.uni <- length(univerID)
   	 	univerID <- intersect(univerID,names(gene2go))
   	 	id.ann.supply <- intersect(supplyID,univerID)
		}
    if(length(id.ann.supply)==0) stop("All supplied IDs are not annotated!")
    supply.goid.self <- unique(as.character(unlist(gene2go[id.ann.supply])))
    supply.goid.aces <- unique(as.character(unlist(all.aces[supply.goid.self])))
    supply.goid <- unique(c(supply.goid.self,supply.goid.aces))
    supply.goid <- supply.goid[-grep("all",supply.goid)]
    rowN <- length(supply.goid)
    res.df <-
      data.frame("GO_ID"=supply.goid,
                 "GO_Term"=as.character(sapply(goOnto[supply.goid],function(x) x[[2]])),
                 "GO_Class"=as.character(sapply(goOnto[supply.goid],function(x) x[[1]])),
                 "Pvalue"=rep(0,rowN),
                 "AdjustedPv"=rep(0,rowN),
                 "x1"=rep(0,rowN),
                 "x2"=rep(0,rowN),
                 "n"=rep(0,rowN),
                 "N"=rep(0,rowN),
                 "EnrichDirect"=rep("",rowN),
                 "GOlevl"=as.character(sapply(goOnto[supply.goid],function(x) x[[3]])),
                 "GeneID"=rep("",rowN))
    rownames(res.df) <- supply.goid
    if(missing(test.method)) test.method="FisherChiSquare"
    for(i in supply.goid){
      inter.goid <- intersect(c(all.off[[i]],i),names(go2gene.bottom))
      id.clst.uni <- as.character(unlist(go2gene.bottom[inter.goid]))
      id.clst.uni <- intersect(id.clst.uni,univerID)
      num.clst.uni <- length(id.clst.uni)
      if(num.clst.uni == 0){
        res.df[i,4:9] <- rep(0,6)
        res.df[i,10:12] <- rep("a",3)
        next
      }
      id.clst.supply <- intersect(id.ann.supply,id.clst.uni)
      num.clst.supply <- length(id.clst.supply)
      res.df[i,6:9] <- c(num.clst.supply,num.clst.uni,num.tot.supply,num.tot.uni)
      dif.data <- c(num.clst.supply,num.clst.uni - num.clst.supply,
                    num.tot.supply-num.clst.supply,
                    num.tot.uni - num.tot.supply-(num.clst.uni-num.clst.supply))
      mt <- matrix(dif.data,nr=2)

      test.pv <- switch(test.method,
                        FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                        HyperGeometric=hyper.test(mt)
                        )
      res.df[i,"Pvalue"] <- test.pv
      res.df[i,"GeneID"] <- paste(id.clst.supply,collapse=",")
      if(num.clst.supply/num.tot.supply > num.clst.uni/num.tot.uni)
        res.df[i,"EnrichDirect"] <- "Over"
      else
        res.df[i,"EnrichDirect"] <- "Under"
    }
    ##filter the rows without ids
    res.df <- res.df[res.df[,"N"]>0,]
    if(missing(p.adjust.method)) p.adjust.method="fdr"
    if(missing(filt)) filt <- "adjp"
    res.df[,"AdjustedPv"] <- p.adjust(res.df[,"Pvalue"],method=p.adjust.method)
    res.df <- sort.data.frame(res.df,key="Pvalue")
    options(opt.str)
    write.table(res.df,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.filt <- paste(enrichFile,".filt",sep="")
    if(filt=="adjp")
      res.df.filt <- res.df[res.df[,"AdjustedPv"] < pc & res.df[,"EnrichDirect"]=="Over",]
    else
      res.df.filt <- res.df[res.df[,"Pvalue"] < pc & res.df[,"EnrichDirect"]=="Over",]
    write.table(res.df.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.merg <- paste(enrichFile.filt,".merg",sep="")
    if(dim(res.df.filt)[1] >= 1){
      diffgo.merg <- MergGO(res.df.filt)
      write.table(diffgo.merg,file=enrichFile.merg,sep="\t",quote=FALSE,row.names=FALSE)
    }
    return(res.df)
  }
