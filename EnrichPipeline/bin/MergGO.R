##used to merge the results that two GO share the same genes
MergGO <- function(res.dt){
  go.res <- list()
  go.term <- res.dt[,1]
  go.len <- length(go.term)
  go.level <- res.dt[,"GOlevl"]
  names(go.level) <- res.dt[,"GO_ID"]
  go2ids <- sapply(res.dt[,"GeneID"],function(x) strsplit(x,",")[[1]])
  names(go2ids) <- res.dt[,"GO_ID"]
  if(dim(res.dt)[1]==1) return(res.dt)
  gotorem <- vector()
  gotorem.pair <- matrix(0,nr=2,nc=2)
  for(i in 1:(go.len-1)){
    go1 <- names(go.level)[i]
    if(go.level[[go1]]==1) gotorem <- append(gotorem,go1)
    for(j in (i+1):go.len){
      go2 <- names(go.level)[j]
      if(go.level[[go2]]==1) gotorem <- append(gotorem,go2)
      interlen <- length(intersect(go2ids[[go1]],go2ids[[go2]]))
      ttlen1 <- length(go2ids[[go1]]);ttlen2 <- length(go2ids[[go2]])
      if(identical(interlen,ttlen1) & identical(interlen,ttlen2)){
        if(is.element(go1,all.aces[[go2]])){ gotorem <- append(gotorem,go1);gotorem.pair <- rbind(gotorem.pair,c(go1,go2))}
        if(is.element(go2,all.aces[[go1]])){ gotorem <- append(gotorem,go2);gotorem.pair <- rbind(gotorem.pair,c(go2,go1))}
      }
    }
  }
  gotorem <- unique(gotorem)
  if(length(gotorem)>0) go.res <- res.dt[-match(gotorem,res.dt[,1]),]
  else go.res <- res.dt
  return(go.res)
}
