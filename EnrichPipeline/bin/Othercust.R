
bp.lev2 <- unique(as.character(bpchild[["GO:0008150"]]))
bp.lev2.go <- sapply(bp.lev2,function(x){
	paste(x,goOnto[[x]][2],sep="\t")
})
write(bp.lev2.go,file="GO_lev2_BP.txt")

cc.lev2 <- unique(as.character(ccchild[["GO:0005575"]]))
cc.lev2.go <- sapply(cc.lev2,function(x){
	paste(x,goOnto[[x]][2],sep="\t")
})
write(cc.lev2.go,file="GO_lev2_CC.txt")
mf.lev2 <- unique(as.character(mfchild[["GO:0003674"]]))
mf.lev2.go <- sapply(mf.lev2,function(x){
	paste(x,goOnto[[x]][2],sep="\t")
})
write(mf.lev2.go,file="GO_lev2_MF.txt")

