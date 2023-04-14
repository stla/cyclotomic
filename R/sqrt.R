eb <- function(n) {
  if(n == 1L) {
    return(0L)
  }
  en <- zeta(n)
  rng <- 1L:((n-1L) %/% 2L)
  toadd <- lapply(rng, function(k) {
    powerCyc(en, (k*k) %% n)
  })
  Reduce(sumCyc, toadd)
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
    "1" = fromInteger(fact) * (2L * eb(nn) + 1L),
    "2" = fromInteger(fact) *
      (zeta(8L) - zeta(8L)^3L) * sqrtPositiveInteger(nn %/% 2L),
    "3" = fromInteger(-fact) * zeta(4L) * (2L * eb(nn) + 1L),
    "0" = fromInteger(2L*fact) * sqrtPositiveInteger(nn %/% 4L)
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
