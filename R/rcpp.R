#' @useDynLib cyclotomic, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

R_gcdList <- function(ivec) {
  do.call(Rgcd, as.list(ivec))
}
