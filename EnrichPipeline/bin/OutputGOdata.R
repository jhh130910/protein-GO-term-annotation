##output the table of all the GO information
godt <- sapply(goOnto,function(x){
  as.character(x)
})
godt <- t(godt)
godt <- cbind(rownames(godt),godt[,1:3])
colnames(godt) <- c("GOID","Ontology","Definition","Level")

write.table(godt,file="MetaGO_table.txt",sep="\t",quote=FALSE,
            row.names=FALSE)
            
