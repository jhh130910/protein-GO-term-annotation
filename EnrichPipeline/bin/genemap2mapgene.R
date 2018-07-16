##this function will convert the gene2map relation to map2gene relation and write the result to file.mapgene
genemap2mapgene <- function(file.genemap,file.mapgene){
  file.genemap <- "cqu_gene_map.tab"
  file.map2gene <- "cqu_map2gene.tab"
  cc <- scan(file.genemap,what="character",sep="\n")
  gene2map <- sapply(cc,function(x) {
    x1 <- strsplit(x,"\t")[[1]][-1]
    x2 <- strsplit(x1," ")[[1]]
    return(x2)
  })
  names(gene2map) <- as.character(sapply(cc,function(x){
    strsplit(x,"\t")[[1]][1]
  }))
  map2gene <- tapply(rep(names(gene2map),sapply(gene2map,length)),
                     as.character(unlist(gene2map)),as.character)
  map2gene.tab <- sapply(names(map2gene),function(x){
    paste(c(x,map2gene[[x]]),collapse="\t")
  })
  write(map2gene.tab,file=file.map2gene)
}
