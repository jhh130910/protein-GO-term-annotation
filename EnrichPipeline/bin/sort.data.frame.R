sort.data.frame <- function(x, key, ...) {
    if (missing(key)) {
        rn <- rownames(x)
        if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
        x[order(rn, ...), , drop=FALSE]
    } else {
        x[do.call("order", c(x[key], ...)), , drop=FALSE]
    }
}