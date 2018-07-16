Plot.GOslim <- function(config.file,slim.file,graph.file.pdf){
  config <- read.table(config.file,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  conf.dm <- dim(config)
  go.slim <- read.table(slim.file,sep="\t",quote="",stringsAsFactors=FALSE)[,2]
  names(go.slim) <- read.table(slim.file,sep="\t",quote="",stringsAsFactors=FALSE)[,1]
  layout.num <- conf.dm[1]
  go.dt <- data.frame()
  dif.go <- vector("character")
  for(i in 1:conf.dm[1]){
    for(j in 1:conf.dm[2]){
      dt <- read.table(config[i,j],sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
      colnames(dt) <- dt[1,]
      rownames(dt) <- dt[,1]
      dt[dt[,"EnrichDirect"]=="Under","AdjustedPv"] <- 1
      dt[dt[,"AdjustedPv"]=="0","AdjustedPv"] <- "1e-16"
      dt <- dt[-1,-1]
      mt <- data.frame(go.term=as.character(go.slim),
                       log.p=rep(0,length(go.slim)),
                       var=rep(rownames(config)[i],length(go.slim)),
                       time=rep(colnames(config)[j],length(go.slim)))
      rownames(mt) <- names(go.slim)
      go.inter <- intersect(rownames(mt),rownames(dt))
      dt[,"AdjustedPv"] <- as.numeric(dt[,"AdjustedPv"])
      mt[go.inter,"log.p"] <- -log10(dt[go.inter,"AdjustedPv"])
      pv.go <- go.inter[mt[go.inter,"log.p"] > -log10(0.05)]
      dif.go <- append(dif.go,pv.go)
      if(dim(go.dt)[1]==0)
        go.dt <- mt
      else
        go.dt <- rbind(go.dt,mt)
    }
  }
  dif.go <- unique(dif.go)
  go.dt <- go.dt[go.dt[,"go.term"] %in% as.character(go.slim[dif.go]),]
	#go.ord <- scan("go_ord.txt",sep="\n",what="character")
	#go.dt$go.term <- factor(go.dt$go.term,leveles=go.ord)
  library(lattice)
 ## library(latticeExtra)
  pdf(graph.file.pdf,width=13,height=10)
  par.settings=list(layout.heights=list(top.padding=-1))
  bar.pt <-
    barchart(log.p ~ go.term|var,data=go.dt,groups=time,
             layout=c(1,layout.num),ylim=c(0,15),
             scales=list(x=list(rot=60,cex=0.9,font=2)),
             ylab=list("-log10(adjustedPv)",font=2),
             ##axis=axis.grid,
             auto.key=list(space="top",columns=2)
             )
  print(bar.pt)
  dev.off()
}
