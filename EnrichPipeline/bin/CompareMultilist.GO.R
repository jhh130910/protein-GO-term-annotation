CompareTwolist.GO <-
    function(conf.file,
             p.adjust.method=c("holm","hochberg","hommel","bonferroni","BH"
             ,"BY","fdr","none"),enrichFile="enrichOut.difgoall",
             test.method=c("FisherChiSquare","HyperGeometric"),filt=c("p","adjp"),pc=0.05
             ){
        opt.str <- options("stringsAsFactors")
        options(stringsAsFactors=FALSE)
        ids.list <- list()
        dt <- read.table(conf.file,header=FALSE,sep="\t")
        for(i in 1:dim(dt)[1]){
            ids.list[[dt[i,1]]] <- scan(dt[i,2],what="character")
        }
        supplyID.list <- sapply(ids.list,function(x) intersect(x,names(gene2go)))
        if(any(sapply(supplyID.list,length)==0)) stop("Supplied IDs are not annotated!")

        go2gene.bottom.list <- list()
        for(i in names(supplyID.list)){
            go2gene.bottom.list[i] <-
                tapply(rep(names(gene2go[supplyID.list[[i]]]),sapply(gene2go[supplyID.list[[i]]],length)),
                       as.character(unlist(gene2go[supplyID.list[[i]]])),unique)
        }

        goid.self.list <- sapply(supplyID.list,function(x) unique(as.character(unlist(gene2go[supplyID.list[[x]]]))))
        goid.aces <- sapply(supplyID.list,function(x) unique(as.character(unlist(all.aces[supplyID.list[[x]]]))))
        goid.list <- sapply(names(supplyID.list),function(x){
            x1 <- unique(c(goid.self.list[[x]],goid.aces[[x]]))
            x1[-grep("all",x1)]
        })

        supply.goid <- unique(as.character(unlist(goid.list)))
        rowN <- length(supply.goid)

        num.dt <- matrix(0,nr=rowN,nc=2*length(supplyID.list))
        colnames(num.dt) <- c(paste("x",names(supplyID.list),sep="_"),paste("n",names(supplyID.list),sep="_"))
        GeneID.dt <- matrix("",nr=rowN,nc=length(supplyID.list))
        colnames(GeneID.dt) <- paste("GeneID",names(supplyID.list),sep="_")

        ##initialize the result table
        res.df <-
            data.frame("GO_ID"=supply.goid,
                       "GO_Term"=as.character(sapply(goOnto[supply.goid],function(x) x[[2]])),
                       "GO_Class"=as.character(sapply(goOnto[supply.goid],function(x) x[[1]])),
                       "Pvalue"=rep(0,rowN),
                       "AdjustedPv"=rep(0,rowN),
                       "WhichOver"=rep("",rowN),
                       "GOlevl"=as.character(sapply(goOnto[supply.goid],function(x) x[[3]])),
                       num.dt,GeneID.dt)
        rownames(res.df) <- supply.goid
        colnames(res.df) <- c("GO_ID","GO_Term","GO_Class","Pvalue","AdjustedPv","WhichOver","GOlevl",
                              colnames(num.dt),colnames(GeneID.dt))
        if(missing(test.method)) test.method="FisherChiSquare"
        num.tot.list <- sapply(supplyID.list,function(x) length)

        for(i in supply.goid){
            id.clst.list <- sapply(names(go2gene.bottom.list),function(x){
                    x1 <- intersect(c(all.off[[i]],i),names(go2gene.bottom.list[[x]]))
                    x2 <- unique(as.character(unlist(go2gene.bottom.list[[x]][x1])))
                    x2
                })
            id.clst.num <- sapply(id.clst.list,length)
            if(all(id.clst.num==0)){
                res.df[i,"Pvalue"] <- 1
                res.df[i,colnames(num.dt)] <- c(rep(0,length(num.tot.list)),num.tot.list)
                res.df[i,colnames(GeneID.dt)] <- rep("",length(num.tot.list))
                res.df[i,"WhichOver"] <- "NA"
            }
            if(dim(dt)[1] ==2 ){
                dif.data <- c(id.clst.num[1:2],num.tot.list[1]-id.clst.num[1],num.tot.list[2]-id.clst.num[2])
                mt <- matrix(dif.data,nr=2)
                test.pv <- switch(test.method,
                                  FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                                  HyperGeometric=hyper.test(mt))
            }else{
                test.pv <- chisq.test(id.clst.num,num.tot.list/sum(num.tot.list))$p.v
            }
            res.df[i,"Pvalue"] <- test.pv
            res.df[i,colnames(num.dt)] <- c(id.clst.num,num.tot.list)
            res.df[i,colnames(GeneID.dt)] <- as.character(sapply(id.clst.list,function(x) paste(x,collapse=",")))
            res.df[i,"WhichOver"] <- names(num.tot.list)[which.max(id.clst.num/num.tot.list)]
        }
        ##filter the rows without ids
        res.df <- res.df[apply(res.df[,paste("x",names(supplyID.list),sep="_")],1,function(x) any(x>0)),]
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
