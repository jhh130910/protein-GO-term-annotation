Fisher.Chi.test <- function(mt){
  r.sum <- rowSums(mt)
  c.sum <- colSums(mt)
  t.sum <- sum(mt)
  exp <- c(r.sum[1]*c.sum[1]/t.sum,r.sum[2]*c.sum[2]/t.sum,
           r.sum[1]*c.sum[2]/t.sum,r.sum[2]*c.sum[1]/t.sum)
  if(any(exp<=5)) test.res <- fisher.test(mt)
  if(all(exp>5)) test.res <- chisq.test(mt)
  test.res
}
