eb <- function(n) {
  if(n == 1L) {
    return(0L)
  }
  en <- e(n)
  rng <- 1L:((n-1L) %/% 2L)
  toadd <- lapply(rng, function(k) {
    powerCyc(en, (k*k) %% n)
  })
  Reduce(sumCyc, toadd)
}

sqrtPositiveInteger <- function(n) {
  fctrs <- factorise(n)
  primes <- as.integer(fctrs[["primes"]])
  powers <- fctrs[["k"]]
  fact <- prod(primes^(powers %/% 2L))
  nn   <- prod(primes^(powers %% 2L))
  switch(
    nn %% 4L,
    "1" = fromInteger(fact) * (2L * eb(nn) + 1L),
    "2" = fromInteger(fact) * (e(8L) - e(8L)^3L) * sqrtPositiveInteger(nn %/% 2L),
    "3" = fromInteger(-fact) * e(4L) * (2L * eb(nn) + 1L),
    "0" = fromInteger(2L*fact) * sqrtPositiveInteger(nn %/% 4L)
  )
}

sqrtInteger <- function(n) {
  stopifnot(isInteger(n))
  n <- as.integer(n)
  if(n == 0L) {
    zeroCyc()
  } else if(n < 0L) {
    e(4L) * sqrtPositiveInteger(-n)
  } else {
    sqrtPositiveInteger(n)
  }
}

sqrtRational <- function(rat) {
  num <- numerator(rat)
  den <- denominator(rat)
  prodRatCyc(1L / den, sqrtInteger(as.integer(num * den)))
}

cycSqrt <- function(x) {

}
