EnrichGO <- function(supplyID,univerID=names(gene2go),annoGene=names(gene2go),
                     GOclass=c("BP","CC","MF"),
                     p.adjust.method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                     enrichFile="enrichOutfile.difgo3",
                     test.method=c("FisherChiSquare","HyperGeometric"),filt=c("p","adjp"),pc=0.05,Unanno=c("Y","N"))
{
	if(missing(p.adjust.method)) p.adjust.method="fdr"
    if(missing(filt)) filt <- "adjp"
	ids <- intersect(supplyID,annoGene)
	if(missing(Unanno) | identical(Unanno,"N")){
		univerID <- intersect(univerID,annoGene)
		annoNum <- length(univerID)
		ids.num <- length(ids)
	}else{
		annoNum <- length(univerID)
		ids.num <- length(supplyID)
	}
    if(length(ids)==0) stop("All supplied IDs are not annotated!")
    opt.str <- options("stringsAsFactors")
    options(stringsAsFactors=FALSE)
	if(length(annoGene)>length(univerID)){
		genefreq.mf.uni <- GOfreq(genefreq.mf,gene2go[univerID],"MF",mf.aces,goOnto)
		genefreq.cc.uni <- GOfreq(genefreq.cc,gene2go[univerID],"CC",cc.aces,goOnto)
		genefreq.bp.uni <- GOfreq(genefreq.bp,gene2go[univerID],"BP",bp.aces,goOnto)
		genefreq.mf.uni  <- genefreq.mf.uni[genefreq.mf.uni>0]
		genefreq.cc.uni  <- genefreq.cc.uni[genefreq.cc.uni>0]
		genefreq.bp.uni  <- genefreq.bp.uni[genefreq.bp.uni>0]
	}else {
		genefreq.mf.uni <- genefreq.mf.ref
		genefreq.bp.uni <- genefreq.bp.ref
		genefreq.cc.uni <- genefreq.cc.ref
	}
	genefreq.mf.supply <- GOfreq(genefreq.mf,gene2go[ids],"MF",mf.aces,goOnto)
	genefreq.cc.supply <- GOfreq(genefreq.cc,gene2go[ids],"CC",cc.aces,goOnto)
	genefreq.bp.supply <- GOfreq(genefreq.bp,gene2go[ids],"BP",bp.aces,goOnto)
	genefreq.mf.supply  <- genefreq.mf.supply[genefreq.mf.supply>0]
	genefreq.cc.supply  <- genefreq.cc.supply[genefreq.cc.supply>0]
	genefreq.bp.supply  <- genefreq.bp.supply[genefreq.bp.supply>0]
	diffgo.mf <- DiffGOs(annoNum,ids.num,goOnto,genefreq.mf.supply,genefreq.mf.uni,"MF")
	diffgo.cc <- DiffGOs(annoNum,ids.num,goOnto,genefreq.cc.supply,genefreq.cc.uni,"CC")
	diffgo.bp <- DiffGOs(annoNum,ids.num,goOnto,genefreq.bp.supply,genefreq.bp.uni,"BP")
	diffgo <- rbind(diffgo.mf,diffgo.cc,diffgo.bp)
	diffgo <- diffgo[,c(3,8,9,1,2,4:7)]
	colnames(diffgo) <- c("GO_ID","GO_Term","GO_Class","Pvalue","AdjustedPv","x1","x2","n","N")
	EnrichDirect <- vector("character",dim(diffgo)[1])
	EnrichDirect[which(as.numeric(diffgo[,6])/as.numeric(diffgo[,7])>as.numeric(diffgo[,8])/as.numeric(diffgo[,9]))] <- "Over"
	EnrichDirect[which(as.numeric(diffgo[,6])/as.numeric(diffgo[,7])<as.numeric(diffgo[,8])/as.numeric(diffgo[,9]))] <- "Under"
	Genes <- sapply(diffgo[,1],function(x) paste(intersect(go2gene[[x]],ids),collapse=","))
	diffgo <- cbind(diffgo,EnrichDirect,Genes)
	diffgo <- sort.data.frame(diffgo,key="Pvalue")
	diffgo[,5] <- p.adjust(diffgo[,4],method=p.adjust.method)
        
	diffgo[,6:9] <- apply(diffgo[,6:9],2,as.numeric)
    write.table(diffgo,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
    enrichFile.filt <- paste(enrichFile,".filt",sep="")
    if(filt=="adjp")
      diffgo.filt <- diffgo[diffgo[,"AdjustedPv"] < pc & diffgo[,"EnrichDirect"]=="Over",]
    else
      diffgo.filt <- diffgo[diffgo[,"Pvalue"] < pc & diffgo[,"EnrichDirect"]=="Over",] 
    write.table(diffgo.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
    options(opt.str)
	return(diffgo)
}
