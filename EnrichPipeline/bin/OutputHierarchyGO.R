OutputHierarchyGO <- function(go.ids,out.file){
  ##this function is used to output the hierarchial GO information
  ## input: GO ids at various level
  ## output: one GO id one line, all the parents hierarchically displayed from the first column to the last.
  go.parents <- sapply(go.ids,function(x){
    x.lv <- as.numeric(goOnto[[x]][3])
    x.go <- paste(x,goOnto[[x]][2],sep=" ")
    if(x.lv>1){
      x.par <- as.character(all.parent[[x]])
      x.par.lv.min.loc <- which.min(as.numeric(sapply(x.par,function(y) goOnto[[y]][3])))
      x.par.min.go <- x.par[x.par.lv.min.loc]
      x.par.lv.min <- as.numeric(goOnto[[x.par.min.go]][3])
      x.go <- paste(paste(x.par.min.go,goOnto[[x.par.min.go]][2],sep=" "),x.go,sep="\t")
    }else x.par.lv.min <- 1
    while(x.par.lv.min>1){
      x.par <- as.character(all.parent[[x.par[x.par.lv.min.loc]]])
      x.par.lv.min.loc <- which.min(as.numeric(sapply(x.par,function(y) goOnto[[y]][3])))
      x.par.min.go <- x.par[x.par.lv.min.loc]
      x.par.lv.min <- as.numeric(goOnto[[x.par.min.go]][3])
      x.go <- paste(paste(x.par.min.go,goOnto[[x.par.min.go]][2],sep=" "),x.go,sep="\t")
    }
    x.go
  })
  go.parents <- c(paste("Lev",1:11,collapse="\t",sep=""),go.parents)
  write.table(go.parents,file=out.file,sep="\t",col.names=FALSE,quote=FALSE)
}
