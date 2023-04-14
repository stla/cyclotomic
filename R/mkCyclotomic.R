#' @importFrom VeryLargeIntegers as.vli factors is.vli phi
#' @importFrom primes Rgcd
#' @importFrom maybe is_just is_nothing from_just nothing just
#' @importFrom gmp as.bigq numerator denominator is.bigq asNumeric
NULL

#' @title The primitive n-th root of unity.
#' @description For example, `e(4) = i` is the primitive 4th root of unity,
#'   and `e(5) = exp(2*pi*i/5)` is the primitive 5th root of unity.
#'   In general, `e(n) = exp(2*pi*i/n)`.
#' @param n positive integer
#' @return A cyclotomic number.
#' @export
#' @examples
#' e(4)
e <- function(n) {
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
  out <- trms #$copy()
  if(!out$has_key(r)) {
    return(out)
  }
  minusrat <- -out$get(r)
  out$erase(r)
  rpl <- rev(replacements(n, p, r))
  for(k in rpl) {
    insertWith(`+`, out, k, minusrat)
  }
  out
}

factorise <- function(n) {
  if(n == 1L) {
    return(list(primes = integer(0L), k = integer(0L)))
  }
  fctrs <- factors(as.vli(n), iter = 100L, output = "list")
  if(is.vli(fctrs)) { # workaround pb si un seul facteur
    fctrs <- list(fctrs)
  }
  fcts <- vapply(
    fctrs, as.integer, FUN.VALUE = integer(1L)
  )
  tbl <- table(fcts)
  list(
    primes = as.integer(names(tbl)),
    k      = as.integer(c(tbl))
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
    c(p, as.integer(p^powers[i]))
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
  m <- as.integer((q %/% p - 1L) %/% 2L) # Y'AVAIT ERREUR ICI
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
  x[!duplicated(x), , drop = FALSE]
}

convertToBase <- function(n, trms) {
  out <- trms #$copy()
  if(n == 1L) {
    return(out)
  }
  epows <- extraneousPowers(n)
  for(i in nrow(epows):1L) {
    pr <- epows[i, ]
    out <- replace(n, pr[1L], pr[2L], out)
  }
  out
}

equalReplacements <- function(p, r, cyc) {
  xx <- vapply(replacements(cyc@order, p, r), function(k) {
    as.character(cyc@terms$get(k, default = as.bigq(0L)))
  }, FUN.VALUE = character(1L))
  x1 <- xx[1L]
  for(x in xx[-1L]) {
    if(x != x1) {
      return(NULL)
    }
  }
  return(x1)
}

reduceByPrime <- function(p, cyc) { # p: integer; cyc: cyclotomic; output: cyclotomic
  n <- cyc@order
  cfs <- as.bigq(integer(0L))
  x <- equalReplacements(p, 0L, cyc)
  r <- p
  nminusp <- n - p
  while(r <= nminusp && !is.null(x)) {
    rat <- as.bigq(x)
    cfs <- c(cfs, -rat)
    x <- equalReplacements(p, r, cyc)
    r <- r + p # Y'AVAIT ERREUR
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
    if(coef != 0L) trms$insert(ii, coef) # Y'AVAIT ERREUR
    ii <- ii + 1L
    i <- i + 1L
  }
  new(
    "cyclotomic",
    order = ndivp,
    terms = removeZeros(trms)
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
  return(firstelem)
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

tryReduce <- function(cyc) {
  fctr <- factorise(cyc@order)
  primes <- fctr[["primes"]]
  powers <- fctr[["k"]]
  ok <- powers == 1L & primes != 2L
  squareFreeOddFactors <- primes[ok]
  if(length(squareFreeOddFactors) == 0L) {
    return(cyc)
  }
  Reduce(reduceByPrime, rev(squareFreeOddFactors), init = cyc, right = TRUE)
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
    as.integer(ord),
    convertToBase(as.integer(ord), trms)
  )
}
