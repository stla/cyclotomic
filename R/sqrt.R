eb <- function(n) { # n is always odd
  if(n == 1L) {
    return(zeroCyc())
  }
  en <- zeta(n)
  rng <- 1L:((n-1L) %/% 2L)
  toadd <- lapply(rng, function(k) {
    powerCyc(en, (k*k) %% n)
  })
  Reduce(sumCyc, toadd)
}

zeta4 <- function() {
  new("cyclotomic", order = 4L, terms = intmap$new(1L, list(as.bigq(1L))))
}

sqrt2 <- function() {
  new(
    "cyclotomic", order = 8L,
    terms = intmap$new(c(1L, 3L), list(as.bigq(1L), as.bigq(-1L)))
  )
}

sqrtPositiveInteger <- function(n) {
  n <- as.integer(n)
  fctrs <- factorise(n)
  primes <- as.integer(fctrs[["primes"]])
  powers <- fctrs[["k"]]
  if(length(primes) != 0L) {
    fact <- as.integer(prod(intpow(primes, as.integer(powers %/% 2L))))
    nn   <- as.integer(prod(intpow(primes, powers %% 2L)))
  } else {
    fact <- nn <- 1L
  }
  switch(
    nn %% 4L,
    "1" = prodIntCyc(fact, sumCyc(prodIntCyc(2L, eb(nn)), oneCyc())),
    "2" = prodIntCyc(fact, sqrt2() * sqrtPositiveInteger(nn %/% 2L)),
    "3" = prodIntCyc(
      -fact,
      prodCyc(zeta4(), sumCyc(prodIntCyc(2L, eb(nn)), oneCyc()))
    ),
    "0" = prodIntCyc(2L*fact, sqrtPositiveInteger(nn %/% 4L))
  )
}

sqrtInteger <- function(n) { # n is integer
  if(n == 0L) {
    zeroCyc()
  } else if(n < 0L) {
    zeta(4L) * sqrtPositiveInteger(-n)
  } else {
    sqrtPositiveInteger(n)
  }
}

sqrtRational <- function(rat) {
  num <- numerator(rat)
  den <- denominator(rat)
  prodRatCyc(1L / den, sqrtInteger(as.integer(num * den)))
}

#' @title Square root as a cyclotomic number
#' @description Square root of an integer or a rational number as a cyclotomic
#'   number. This is slow.
#'
#' @param x an integer, a \strong{gmp} rational number (\code{bigq} object), or
#'   a fraction given as a string (e.g. \code{"5/3"})
#'
#' @return The square root of \code{x} as a cyclotomic number.
#' @export
#'
#' @examples
#' cycSqrt(2)
#' phi <- (1 + cycSqrt(5)) / 2 # the golden ratio
#' phi^2 - phi # should be 1
cycSqrt <- function(x) {
  if(isFraction(x) || is.bigq(x)) {
    sqrtRational(x)
  } else if(isInteger(x)) {
    sqrtInteger(as.integer(x))
  } else {
    stop("Invalid argument `x`.")
  }
}
