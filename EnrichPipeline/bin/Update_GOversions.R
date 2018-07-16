#####################
##1. 
##this script is used to construct the MetaGO data
##which parse the GO hierarchal structure
##prerequisite package: GO.db
#####################
setwd("d:/bin/Pipeline/EnrichPipeline/bin")
library(GO.db)
##preparing the meta.data
bpchild <- as.list(GOBPCHILDREN)
bp.level3 <- 
unique(as.character(unlist(bpchild[as.character(bpchild[["GO:0008150"]])])))
ccchild <- as.list(GOCCCHILDREN)
cc.level3 <- 
unique(as.character(unlist(ccchild[as.character(ccchild[["GO:0005575"]])])))
mfchild <- as.list(GOMFCHILDREN)
mf.level3 <- mfchild[as.character(mfchild[["GO:0003674"]])]
naloc <- which(sapply(mf.level3,function(x) any(is.na(x))))
mf.level3[naloc] <- names(naloc)
mf.level3 <- unique(as.character(unlist(mf.level3)))
bp.aces <- as.list(GOBPANCESTOR)
cc.aces <- as.list(GOCCANCESTOR)
mf.aces <- as.list(GOMFANCESTOR)

bp.parent <- as.list(GOBPPARENTS)
cc.parent <- as.list(GOCCPARENTS)
mf.parent <- as.list(GOMFPARENTS)
all.parent <- c(bp.parent,cc.parent,mf.parent)

goterm <- as.list(GOTERM)
##GO level
goLev <- rep("a",length=length(goterm))
names(goLev) <- names(goterm)
Childs <- c("GO:0003674","GO:0008150","GO:0005575")
lev <- 1
all.child <- c(bpchild,ccchild,mfchild)
while(length(Childs)>=1){
  loc <- match(Childs,names(goLev))
  goLev[loc][goLev[loc] == "a"] <- lev
  Childs <- unique(as.character(unlist(all.child[Childs])))
  Childs <- Childs[! is.na(Childs)]
  lev <- lev+1  
}
goOnto <- lapply(names(goterm),function(x) c(Ontology(goterm[[x]]),
	Term(goterm[[x]]),as.character(goLev[x])))
names(goOnto) <- names(goterm)
##initialize gene frequency at level3
genefreq.mf <- vector("numeric",length=length(mf.level3))
names(genefreq.mf) <- mf.level3
genefreq.bp <- vector("numeric",length=length(bp.level3))
names(genefreq.bp) <- bp.level3
genefreq.cc <- vector("numeric",length=length(cc.level3))
names(genefreq.cc) <- cc.level3

source("GOfreq.R")
source("sort.data.frame.R")
source("DiffGOs.R")
source("EnrichGO.R")
source("EnrichGOallLevl.R")
source("Fisher.Chi.test.R")
bp.off <- as.list(GOBPOFFSPRING)
cc.off <- as.list(GOCCOFFSPRING)
mf.off <- as.list(GOMFOFFSPRING)
all.off <- c(bp.off,cc.off,mf.off)
all.aces <- c(bp.aces,cc.aces,mf.aces)


save(all.child,all.aces,all.off,all.parent,bp.level3,cc.level3,mf.level3,bp.aces,cc.aces,mf.aces,goOnto,
     genefreq.mf,genefreq.bp,genefreq.cc,bp.off,cc.off,mf.off,bpchild,ccchild,mfchild,Fisher.Chi.test,
     sort.data.frame,DiffGOs,GOfreq,EnrichGO,EnrichGOallLevl,file="MetaGO_20140308.RData")



###############
##2. 
##output the table of all the GO information
godt <- sapply(goOnto,function(x){
  as.character(x)
})
godt <- t(godt)
godt <- cbind(rownames(godt),godt[,1:3])
colnames(godt) <- c("GOID","Ontology","Definition","Level")

write.table(godt,file="MetaGO_table.txt",sep="\t",quote=FALSE,
            row.names=FALSE)




##########################
##3. 
##output .class file for go_svg.pl
bp.lev2 <- unique(as.character(bpchild[["GO:0008150"]]))
cc.lev2 <- unique(as.character(ccchild[["GO:0005575"]]))
mf.lev2 <- unique(as.character(mfchild[["GO:0003674"]]))
bp.lev2.off <- as.character(unlist(sapply(bp.lev2,function(x) all.off[[x]])))
cc.lev2.off <- as.character(unlist(sapply(cc.lev2,function(x) all.off[[x]])))
mf.lev2.off <- as.character(unlist(sapply(mf.lev2,function(x) all.off[[x]])))
bp.lev2.term <- as.character(sapply(bp.lev2,function(x) goOnto[[x]][2]))
cc.lev2.term <- as.character(sapply(cc.lev2,function(x) goOnto[[x]][2]))
mf.lev2.term <- as.character(sapply(mf.lev2,function(x) goOnto[[x]][2]))

go.lev2.tb <-
  data.frame(Lev1=c(rep("biological_process",length(bp.lev2.off)),
               rep("cellular_component",length(cc.lev2.off)),
               rep("molecular_function",length(mf.lev2.off))),
             Lev2=c(rep(bp.lev2.term,sapply(bp.lev2,function(x) length(all.off[[x]]))),
               rep(cc.lev2.term,sapply(cc.lev2,function(x) length(all.off[[x]]))),
               rep(mf.lev2.term,sapply(mf.lev2,function(x) length(all.off[[x]])))),
             GOID=c(bp.lev2.off,cc.lev2.off,mf.lev2.off),
             GOTerm=c(as.character(sapply(bp.lev2.off,function(x) goOnto[[x]][2])),
               as.character(sapply(cc.lev2.off,function(x) goOnto[[x]][2])),
               as.character(sapply(mf.lev2.off,function(x) goOnto[[x]][2]))))

write.table(go.lev2.tb,file="GO.class",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
