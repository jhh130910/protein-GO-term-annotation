hyper.test <- function(mt){
  r.sum <- rowSums(mt)
  c.sum <- colSums(mt)
  m <- r.sum[1];n <- r.sum[2]
  k <- c.sum[1]
  x <- mt[1,1]
  phyper(x,m,n,k)
}

