#' @useDynLib cyclotomic, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

# These functions are used in Rcpp ####

R_gcdList <- function(ivec) {
  do.call(Rgcd, as.list(ivec))
}

R_phiNrpSqfree <- function(n) {
  if(n == 1L) {
    return(list(
      "phi"    = 1L,
      "nrp"    = 0L,
      "sqfree" = TRUE
    ))
  }
  fctrs <- factors(as.vli(n), iter = 100L, output = "list")
  if(is.vli(fctrs)) { # workaround: fctrs is not a list if only one factor
    fctrs <- list(fctrs)
  }
  fcts <- vapply(
    fctrs, as.integer, FUN.VALUE = integer(1L)
  )
  tbl    <- table(fcts)
  powers <- as.integer(c(tbl))
  list(
    "phi"    = as.integer(phi(as.vli(n))),
    "nrp"    = sum(powers),
    "sqfree" = all(powers == 1L)
  )
}

R_extraneousPowers <- function(n) {
  pairs <- pqPairs(n)
  x <- do.call(rbind, apply(pairs, 2L, function(pq) {
    p <- pq[1L]
    q <- pq[2L]
    r <- Rcpp_removeExps(n, p, q)
    cbind(p, r)
  }, simplify = FALSE))
  x[!duplicated(x), , drop = FALSE]
}

R_squareFreeOddFactors <- function(n) {
  fctr <- factorise(n)
  primes <- fctr[["primes"]]
  powers <- fctr[["k"]]
  ok <- powers == 1L & primes != 2L
  primes[ok]
}
