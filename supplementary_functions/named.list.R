named.list <- function(...) { 
  l <- list(...)
  names(l) <- sapply(substitute(list(...)), deparse)[-1]
  l 
}