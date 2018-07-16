GOfreq <- function(genefreq,ann,OntoClass=c("MF","BP","CC"),aces,goOnto){
	for(i in 1:length(ann)){
		a <- goOnto[ann[[i]]]
		a <- a[sapply(a,length)>1]
		a <- sapply(a,function(x) x[1])
		a <- a[as.character(a)==OntoClass]
		a.aces <- unique(as.character(unlist(aces[names(a)])))[-1] 
		a.aces <- c(a.aces,names(a))
		gointer <- intersect(a.aces,names(genefreq))
		genefreq[gointer] <- genefreq[gointer] + 1
	}
	return(genefreq)
}