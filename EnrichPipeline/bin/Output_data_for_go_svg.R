##output .class file for go_svg.pl
load("../data/MetaGO_20120901.RData")
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
