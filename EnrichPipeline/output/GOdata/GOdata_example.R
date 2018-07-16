
    metago.dir <- "data/MetaGO_200908.RData"
	dir.wego <- "input/gene.wego"
	
##load in the GO Meta data.
	load(metago.dir)
##read in reference GO annotation
	
#########
##GOID to gene products at level 3
	go.anno <- scan(dir.wego,what="character",sep="\n")
	go.anno.nam <- as.character(sapply(go.anno,function(x) strsplit(x,"\t")[[1]][1]))
	go.anno.gos <- sapply(go.anno,function(x) strsplit(x,"\t")[[1]][-1])
	gene2go <- go.anno.gos
	names(gene2go) <- go.anno.nam
##filter the gos without annotation in goOnto
	gene2go.filt <- sapply(gene2go,function(x) intersect(x,names(goOnto)))
	go.loc <- which(sapply(gene2go.filt,length)<1)
	if (length(go.loc)>=1) gene2go <- gene2go.filt[-go.loc] else gene2go <- gene2go.filt
	
	go2gene <- vector("list",length=length(c(mf.level3,bp.level3,cc.level3)))
	names(go2gene) <- c(mf.level3,bp.level3,cc.level3)
	aces.tt <- c(cc.aces,bp.aces,mf.aces)
	for(i in 1:length(gene2go)){
	    a <- gene2go[[i]]
		a.aces <- unique(as.character(unlist(aces.tt[a])))[-1]
		a.aces <- c(a.aces,a)
		gointer <- unique(intersect(a.aces,names(go2gene)))
		if(length(gointer)==0) next
		for(j in gointer){ go2gene[[j]] <- c(go2gene[[j]],names(gene2go)[i])}
#the next line popup error message that:number of items to replace is not a multiple of replacement length
	}
    go2gene <- go2gene[sapply(go2gene,length)>0]
	
##########
##reference
	genefreq.mf.ref <- GOfreq(genefreq.mf,gene2go,"MF",mf.aces,goOnto)
	genefreq.cc.ref <- GOfreq(genefreq.cc,gene2go,"CC",cc.aces,goOnto)
	genefreq.bp.ref <- GOfreq(genefreq.bp,gene2go,"BP",bp.aces,goOnto)
	genefreq.mf.ref  <- genefreq.mf.ref[genefreq.mf.ref>0]
	genefreq.cc.ref  <- genefreq.cc.ref[genefreq.cc.ref>0]
	genefreq.bp.ref  <- genefreq.bp.ref[genefreq.bp.ref>0]
	
	save(all.aces,all.child,all.parent,all.off,bp.level3,cc.level3,mf.level3,bp.aces,cc.aces,mf.aces,goOnto,
	     genefreq.mf,genefreq.bp,genefreq.cc,genefreq.mf.ref,genefreq.cc.ref,
	     genefreq.bp.ref,go2gene,sort.data.frame,Fisher.Chi.test,MergGO,
	     DiffGOs,GOfreq,gene2go,file="output/GOdata/GOdata_example.RData")
	q(save="no")
	
