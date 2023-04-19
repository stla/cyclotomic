#' @importFrom primes gcd scm coprime
NULL

## | the zero cyclotomic number ####
zeroCyc <- function() {
  new("cyclotomic", order = "1", terms = intmap$new())
}

## | check if it's the zero cyclotomic number ####
isZeroCyc <- function(cyc) {
  cyc@terms$size() == 0L
}

## | sum of two cyclotomic numbers ####
sumCyc <- function(cyc1, cyc2) {
  if(isZeroCyc(cyc1)) return(cyc2)
  if(isZeroCyc(cyc2)) return(cyc1)
  o1    <- cyc1@order
  trms1 <- cyc1@terms
  o2    <- cyc2@order
  trms2 <- cyc2@terms
  ord <- scm(o1, o2)
  m1  <- as.integer(ord %/% o1)
  m2  <- as.integer(ord %/% o2)
  mp1 <- mapKeys(
    function(k) { m1 * k }, trms1
  )
  mp2 <- mapKeys(
    function(k) { m2 * k }, trms2
  )
  mkCyclotomic(ord, unionWith(`+`, mp1, mp2))
}

## | product of two cyclotomic numbers ####
prodCyc <- function(cyc1, cyc2) {
  o1    <- cyc1@order
  trms1 <- cyc1@terms
  o2    <- cyc2@order
  trms2 <- cyc2@terms
  ord <- scm(o1, o2)
  m1  <- as.integer(ord %/% o1)
  m2  <- as.integer(ord %/% o2)
  keys1 <- trms1$keys()
  keys2 <- trms2$keys()
  mp <- intmap$new()
  for(k1 in keys1) {
    e1 <- k1
    c1 <- trms1$get(k1)
    for(k2 in keys2) {
      e2 <- k2
      c2 <- trms2$get(k2)
      k <- (m1*e1 + m2*e2) %% ord
      insertWith(`+`, mp, k, c1*c2)
    }
  }
  mkCyclotomic(ord, mp)
}

## | power of a cyclotomic number ####
powerCyc <- function(cyc, n) {
  stopifnot(isInteger(n))
  if(n == 0L) {
    return(fromInteger(1L))
  }
  if(n == 1L) {
    return(cyc)
  }
  if(n > 0L) {
    out <- cyc
    for(. in 2L:n) {
      out <- prodCyc(out, cyc)
    }
    return(out)
  }
  # if n < 0:
  powerCyc(invCyc(cyc), -n)
}

## | product rational and cyclotomic ####
prodRatCyc <- function(rat, cyc) {
  if(rat == 0L) {
    zeroCyc()
  } else {
    new(
      "cyclotomic",
      order = cyc@order,
      terms = mapValues(function(x) {rat * x}, cyc@terms)
    )
  }
}

## | opposite of a cyclotomic number ####
minusCyc <- function(cyc) {
  prodRatCyc(as.bigq(-1L), cyc)
}

## helper 1 for multiplicative inverse
multiplyExponents <- function(j, cyc) { # j is integer
  n <- cyc@order
  if(gcd(j, n) != 1L) {
    stop("multiplyExponents needs gcd == 1")
  }
  mkCyclotomic(
    n,
    mapKeys(
      function(k) {
        (j * k) %% n
      },
      cyc@terms
    )
  )
}

## helper 2 for multiplicative inverse
productOfGaloisConjugates <- function(cyc) {
  ord <- cyc@order
  if(ord <= 2L) return(NULL)
  x <- seq(2L, length.out = ord - 2L)
  coprimes <- as.list(x[coprime(x, ord)])
  tomultiply <- lapply(coprimes, function(j) {
    multiplyExponents(j, cyc)
  })
  Reduce(prodCyc, tomultiply)
}

## | multiplicative inverse ####
invCyc <- function(cyc) {
  pgc <- productOfGaloisConjugates(cyc)
  if(is.null(pgc)) {
    pgc <- fromInteger(1L)
  }
  maybe_rat <- maybeRational(prodCyc(cyc, pgc))
  if(is_just(maybe_rat)) {
    r <- from_just(maybe_rat)
    prodRatCyc(1L/r, pgc)
  } else {
    stop("this is a bug!")
  }
}
