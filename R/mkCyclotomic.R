#' @importFrom VeryLargeIntegers as.vli factors is.vli phi
#' @importFrom primes Rgcd
#' @importFrom maybe is_just is_nothing nothing just
#' @importFrom gmp as.bigq numerator denominator is.bigq asNumeric
NULL

## | integer^integer
intpow <- Vectorize(function(p, k) {
  result <- 1L
  while(k) {
    if(bitwAnd(k, 1L)) {
      result <- result * p
    }
    k <- bitwShiftR(k, 1L)
    p <- p * p
  }
  result
})

#' @title The primitive n-th root of unity.
#' @description For example, `zeta(4) = i` is the primitive 4th root of unity,
#'   and `zeta(5) = exp(2*pi*i/5)` is the primitive 5th root of unity.
#'   In general, `zeta(n) = exp(2*pi*i/n)`.
#' @param n a positive integer
#' @return A cyclotomic number.
#' @export
#' @examples
#' zeta(4)
zeta <- function(n) {
  stopifnot(isStrictlyPositiveInteger(n))
  n <- as.integer(n)
  if(n == 1L) {
    trms <- intmap$new()
    trms$insert(0L, as.bigq(1L))
    new("cyclotomic", order = 1L, terms = trms)
  } else {
    trms <- intmap$new()
    trms$insert(1L, as.bigq(1L))
    cyclotomic(n, convertToBase(n, trms))
  }
}

replacements <- function(n, p, r) {
  n <- as.integer(n)
  s <- as.integer(n %/% p)
  rpl1 <- integer(0L)
  x <- r - s
  while(x >= 0L) {
    rpl1 <- c(rpl1, x)
    x <- x - s
  }
  rpl2 <- integer(0L)
  x <- r + s
  while(x < n) {
    rpl2 <- c(rpl2, x)
    x <- x + s
  }
  c(rpl1, rpl2)
}

replace <- function(n, p, r, trms) {
  if(!trms$has_key(r)) {
    return(trms)
  }
  minusrat <- -trms$get(r)
  trms$erase(r)
  rpl <- rev(replacements(n, p, r))
  for(k in rpl) {
    insertWith(`+`, trms, k, minusrat)
  }
  trms
}

factorise <- function(n) {
  if(n == 1L) {
    return(list(primes = integer(0L), k = integer(0L)))
  }
  fctrs <- factors(as.vli(n), iter = 10L, output = "list")
  if(is.vli(fctrs)) { # workaround: fctrs is not a list if only one factor
    fctrs <- list(fctrs)
  }
  fcts <- vapply(
    fctrs, as.integer, FUN.VALUE = integer(1L)
  )
  tbl <- table(fcts)
  list(
    "primes" = as.integer(names(tbl)),
    "k"      = as.integer(c(tbl))
  )
}

phiNrpSqfree <- function(n) {
  powers <- factorise(n)[["k"]]
  list(
    phi    = as.integer(phi(as.vli(n))),
    nrp    = sum(powers),
    sqfree = all(powers == 1L)
  )
}

pqPairs <- function(n) {
  fctr <- factorise(n)
  primes <- fctr[["primes"]]
  powers <- fctr[["k"]]
  vapply(seq_along(primes), function(i) {
    p <- primes[i]
    c(p, intpow(p, powers[i]))
  }, FUN.VALUE = integer(2L))
}

includeMods <- function(n, q, start) { # n, q and start are integers
  out <- start
  x <- start - q
  while(x >= 0L) {
    out <- c(out, x)
    x <- x - q
  }
  x <- start + q
  while(x < n) {
    out <- c(out, x)
    x <- x + q
  }
  out
}

removeExps <- function(n, p, q) {
  n <- as.integer(n)
  q <- as.integer(q)
  f <- function(start) includeMods(n, q, start)
  ndivq <- n %/% q
  if(p == 2L) {
    x <- as.integer(q %/% 2L)
    out <- f(ndivq * x)
    qm1 <- q - 1L
    while(x < qm1) {
      x <- x + 1L
      out <- c(out, f(ndivq * x))
    }
    return(out)
  }
  m <- as.integer((q %/% p - 1L) %/% 2L)
  x <- - m
  out <- f(ndivq * x)
  while(x < m) {
    x <- x + 1L
    out <- c(out, f(ndivq * x))
  }
  out
}

extraneousPowers <- function(n) {
  pairs <- pqPairs(n)
  x <- do.call(rbind, apply(pairs, 2L, function(pq) {
    p <- pq[1L]
    q <- pq[2L]
    r <- removeExps(n, p, q)
    cbind(p, r)
  }, simplify = FALSE))
  return(x)
  #x[!duplicated(x), , drop = FALSE] # ça revient au même de faire unique(r)
}

convertToBase <- function(n, trms) {
  if(n == 1L) {
    return(trms)
  }
  epows <- extraneousPowers(n)
  for(i in nrow(epows):1L) {
    pr <- epows[i, ]
    trms <- replace(n, pr[1L], pr[2L], trms)
  }
  trms
}

equalReplacements <- function(p, r, cyc) {
  xs <- vapply(replacements(cyc@order, p, r), function(k) {
    as.character(cyc@terms$get(k, default = as.bigq(0L)))
  }, FUN.VALUE = character(1L))
  x1 <- xs[1L]
  for(x in xs[-1L]) {
    if(x != x1) {
      return(NULL)
    }
  }
  x1
}

reduceByPrime <- function(p, cyc) { # p: integer; cyc: cyclotomic; output: cyclotomic
  n <- cyc@order
  nminusp <- n - p
  cfs <- as.bigq(integer(0L))
  r <- 0L
  x <- equalReplacements(p, r, cyc)
  while(r <= nminusp && !is.null(x)) {
    rat <- as.bigq(x)
    cfs <- c(cfs, -rat)
    r <- r + p
    x <- equalReplacements(p, r, cyc)
  }
  if(is.null(x)) {
    return(cyc)
  }
  ndivp <- as.integer(n %/% p)
  trms <- intmap$new()
  ii <- 0L
  i <- 1L
  l <- length(cfs)
  while(ii < ndivp && i <= l) {
    coef <- cfs[i]
    if(coef != 0L) trms$insert(ii, coef)
    ii <- ii + 1L
    i <- i + 1L
  }
  new(
    "cyclotomic",
    order = ndivp,
    terms = trms
  )
}

gcdList <- function(lst) {
  do.call(Rgcd, lst)
}

gcdCyc <- function(cyc) {
  gcdList(as.list(c(cyc@order, cyc@terms$keys())))
}

gcdReduce <- function(cyc) {
  d <- gcdCyc(cyc)
  if(d == 1L) {
    cyc
  } else {
    f <- function(n) {
      as.integer(n %/% d)
    }
    neworder <- f(cyc@order)
    newterms <- mapKeys(f, cyc@terms)
    new("cyclotomic", order = neworder, terms = newterms)
  }
}

removeZeros <- function(mp) {
  filterMap(function(x) {x != 0L}, mp)
}

lenCyc <- function(cyc) {
  trms <- removeZeros(cyc@terms)
  trms$size()
}

equalCoefficients <- function(cyc) {
  trms <- cyc@terms
  if(trms$size() == 0L) {
    return(NULL)
  }
  keys <- trms$keys()
  firstelem <- trms$get(keys[1L])
  for(key in keys[-1L]) {
    elem <- trms$get(key)
    if(elem != firstelem) {
      return(NULL)
    }
  }
  firstelem
}

tryRational <- function(cyc) {
  pns <- phiNrpSqfree(cyc@order)
  if(pns[["sqfree"]] && lenCyc(cyc) == pns[["phi"]]) {
    rat <- equalCoefficients(cyc)
    if(is.null(rat)) {
      cyc
    } else {
      fromRational(if(pns[["nrp"]] %% 2L == 0L) rat else -rat)
    }
  } else {
    cyc
  }
}

squareFreeOddFactors <- function(n) {
  fctr <- factorise(n)
  primes <- fctr[["primes"]]
  powers <- fctr[["k"]]
  primes[powers == 1L & primes != 2L]
}

tryReduce <- function(cyc) {
  sfoFactors <- squareFreeOddFactors(cyc@order)
  if(length(sfoFactors) == 0L) {
    return(cyc)
  }
  Reduce(reduceByPrime, sfoFactors, init = cyc, right = TRUE)
}

cyclotomic <- function(ord, trms) {
  cyc <- new(
    "cyclotomic",
    order = ord,
    terms = removeZeros(trms)
  )
  tryReduce(tryRational(gcdReduce(cyc)))
}

mkCyclotomic <- function(ord, trms) {
  cyclotomic(
    ord, convertToBase(ord, trms)
  )
}
