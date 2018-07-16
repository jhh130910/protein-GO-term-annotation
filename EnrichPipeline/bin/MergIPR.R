##used to merge the results that two IPR share the same genes
MergIPR <- function(res.dt,pair.file){
	ipr <- scan(pair.file,what="character",sep="\n")
	ipr.aces <- list()
	ipr.hirar <- vector("character")
	for(i in ipr){
		i.num <- nchar(strsplit(i,"I")[[1]][1])
		i.ipr <- gsub("\\-+","",strsplit(i,"::")[[1]][1])
		if(i.num == 0){
			ipr.hirar[1] <- i.ipr
		}else if (i.num==2){
			ipr.hirar[2] <- i.ipr
			ipr.aces[[i.ipr]] <- ipr.hirar[1]
		}else if (i.num==4){
			ipr.hirar[3] <- i.ipr
			ipr.aces[[i.ipr]] <- ipr.hirar[1:2]
		}else if (i.num==6){
			ipr.hirar[4] <- i.ipr
			ipr.aces[[i.ipr]] <- ipr.hirar[1:3]
		}else if (i.num==8){
			ipr.hirar[5] <- i.ipr
			ipr.aces[[i.ipr]] <- ipr.hirar[1:4]
		}
	}

  go.res <- list()
  go.term <- res.dt[,1]
  go.len <- length(go.term)
  go2ids <- sapply(res.dt[,"GeneIDs"],function(x) strsplit(x," ")[[1]])
  names(go2ids) <- res.dt[,"IPRID"]
	names(go.term) <- res.dt[,"IPRID"]
  if(dim(res.dt)[1]==1) return(res.dt)
  gotorem <- vector()
  gotorem.pair <- matrix(0,nr=2,nc=2)
  for(i in 1:(go.len-1)){
    go1 <- names(go.term)[i]
    for(j in (i+1):go.len){
      go2 <- names(go.term)[j]
      interlen <- length(intersect(go2ids[[go1]],go2ids[[go2]]))
      ttlen1 <- length(go2ids[[go1]]);ttlen2 <- length(go2ids[[go2]])
      if(identical(interlen,ttlen1) & identical(interlen,ttlen2)){
        if(is.element(go1,ipr.aces[[go2]])){ gotorem <- append(gotorem,go1);gotorem.pair <- rbind(gotorem.pair,c(go1,go2))}
        if(is.element(go2,ipr.aces[[go1]])){ gotorem <- append(gotorem,go2);gotorem.pair <- rbind(gotorem.pair,c(go2,go1))}
      }
    }
  }
  gotorem <- unique(gotorem)
  if(length(gotorem)>0) go.res <- res.dt[-match(gotorem,res.dt[,1]),]
  else go.res <- res.dt
  return(go.res)
}
