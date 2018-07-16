ExtractGOslim <- function(slim.file,out.file){
  lin <- scan(slim.file,what="character",sep="\n")
  go.id <- lin[grep("^id: GO",lin)]
  go.id <- gsub("id: ","",go.id)
  term.id <- lin[grep("name: ",lin)]
  term.id <- gsub("name: ","",term.id)[1:length(go.id)]
  write.table(cbind(go.id,term.id),file=out.file,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
