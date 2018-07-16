DiffGOs <- function(totNum,geneNum,goOnto,freq.supply,freq.ref,OnClass,
                    test.method=c("FisherChiSquare","HyperGeometric") #hypothesis test methods
                    ){
	dif.tbl <- matrix(0,nr=length(freq.supply),nc=9)
	dif.tbl <- as.data.frame.matrix(dif.tbl)
	ini <- 0
	for(i in names(freq.supply)){
		ini <- ini+1
        x1 <- as.numeric(freq.supply[i])
        x2 <- as.numeric(freq.ref[i])
        n <- geneNum
        N <- totNum
		dif.data <- c(x1,x2-x1,n-x1,N-n-(x2-x1))
		dif.tbl[ini,c(3:8)] <- c(i,freq.supply[i],freq.ref[i],geneNum,totNum,goOnto[[i]][2])
        mt <- matrix(dif.data,nr=2)
        if(missing(test.method)) test.method="FisherChiSquare"
        test.pv <- switch(test.method,
                          FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                          HyperGeometric=hyper.test(mt)
                          )
		dif.tbl[ini,1] <- test.pv
	}
        if(length(freq.supply)==0) {
          return(dif.tbl)
        }else{
          dif.tbl[,9] <- OnClass
          return(dif.tbl)
        }
}
