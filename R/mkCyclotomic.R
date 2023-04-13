#' @importFrom VeryLargeIntegers as.vli factors gcd lcmul is.vli phi
#' @importFrom maybe is_just is_nothing from_just
#' @importFrom gmp as.bigq
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
    trms <- fastmap()
    trms$set("0", as.bigq(1L))
    new("cyclotomic", order = "1", terms = trms)
  } else {
    trms <- fastmap()
    trms$set("1", as.bigq(1L))
    cyclotomic(n, convertToBase(n, trms))
  }
}

replacements <- function(n, p, r) {
  n <- as.vli(n)
  s <- n / as.vli(p)
  r <- as.vli(r)
  rpl1 <- character(0L)
  x <- r - s
  while(x >= 0L) {
    rpl1 <- c(rpl1, as.character(x))
    x <- x - s
  }
  rpl2 <- character(0L)
  x <- r + s
  while(x < n) {
    rpl2 <- c(rpl2, as.character(x))
    x <- x + s
  }
  c(rpl1, rpl2)
}

replace <- function(n, p, r, trms) {
  out <- trms$clone()
  if(!out$has(r)) {
    return(out)
  }
  minusrat <- -out$get(r)
  out$remove(r)
  rpl <- rev(replacements(n, p, r))
  for(k in rpl) {
    insertWith(`+`, out, k, minusrat)
  }
  out
}

factorise <- function(n) {
  n <- as.vli(n)
  if(n == 1L) {
    return(list(primes = character(0L), k = integer(0L)))
  }
  fctrs <- factors(n, iter = 100L, output = "list")
  if(is.vli(fctrs)) { # workaround pb si un seul facteur
    fctrs <- list(fctrs)
  }
  fcts <- vapply(
    fctrs, as.character, FUN.VALUE = character(1L)
  )
  tbl <- table(fcts)
  list(
    primes = names(tbl),
    k      = c(tbl)
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
    c(p, as.character(as.vli(p)^powers[i]))
  }, FUN.VALUE = character(2L))
}

includeMods <- function(n, q, start) { # n and q are vli
  out <- start
  start <- as.vli(start)
  x <- start - q
  while(x >= 0L) {
    out <- c(out, as.character(x))
    x <- x - q
  }
  x <- start + q
  while(x < n) {
    out <- c(out, as.character(x))
    x <- x + q
  }
  out
}

removeExps <- function(n, p, q) {
  n <- as.vli(n)
  q <- as.vli(q)
  f <- function(start) includeMods(n, q, start)
  ndivq <- n / q
  if(p == "2") {
    x <- q / 2L
    out <- f(as.character(ndivq * x))
    qm1 <- q - 1L
    while(x < qm1) {
      x <- x + 1L
      out <- c(out, f(as.character(ndivq * x)))
    }
    return(out)
  }
  p <- as.vli(p)
  m <- (q / (p - 1L)) / 2L
  x <- 0L - m
  out <- f(as.character(ndivq * x))
  while(x < m) {
    x <- x + 1L
    out <- c(out, f(as.character(ndivq * x)))
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
  out <- trms$clone()
  if(n == 1L) {
    return(out)
  }
  epows <- extraneousPowers(n)
  for(i in 1L:nrow(epows)) {
    pr <- epows[i, ]
    out <- replace(n, pr[1L], pr[2L], out)
  }
  out
}

equalReplacements <- function(p, r, cyc) {
  xx <- vapply(replacements(cyc@order, p, r), function(k) {
    as.character(cyc@terms$get(k, missing = as.bigq("0")))
  }, FUN.VALUE = character(1L))
  x1 <- xx[1L]
  for(x in xx[-1L]) {
    if(x != x1) {
      return(nothing())
    }
  }
  just(x1)
}

reduceByPrime <- function(p, cyc) { # p: integer; cyc: cyclotomic; output: cyclotomic
  n <- as.vli(cyc@order)
  cfs <- as.bigq(integer(0L))
  x <- equalReplacements(p, "0", cyc)
  if(is_just(x)) {
    rat <- as.bigq(from_just(x))
    cfs <- c(cfs, -rat)
  }
  r <- p <- as.vli(p)
  nminusp <- n - p
  while(is_just(x) && r <= nminusp) {
    rat <- as.bigq(from_just(x))
    cfs <- c(cfs, -rat)
    x <- equalReplacements(p, as.character(r), cyc)
    r <- r + 1L
  }
  if(is_nothing(x)) {
    return(cyc)
  }
  ndivp <- n / p
  trms <- fastmap()
  ii <- as.vli("0")
  i <- 1L
  l <- length(cfs)
  while(ii < ndivp && i <= l) {
    coef <- cfs[i]
    if(coef != 0L) trms$set(as.character(i), coef)
    ii <- ii + 1L
    i <- i + 1L
  }
  new(
    "cyclotomic",
    order = as.character(ndivp),
    terms = removeZeros(trms)
  )
}

gcdList <- function(lst) {
  n  <- lst[1L]
  ns <- lst[-1L]
  Reduce(
    function(x, y) {
      as.character(gcd(as.vli(x), as.vli(y)))
    },
    x = ns,
    init = n,
    right = TRUE
  )
}

gcdCyc <- function(cyc) {
  keys <- cyc@terms$keys()
  gcdList(c(cyc@order, keys))
}

gcdReduce <- function(cyc) {
  d <- as.vli(gcdCyc(cyc))
  if(d == 1L) {
    cyc
  } else {
    f <- function(n) {
      as.character(as.vli(n) / d)
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
    return(nothing())
  }
  keys <- trms$keys()
  firstelem <- trms$get(keys[1L])
  for(key in keys[-1L]) {
    elem <- trms$get(key)
    if(elem != firstelem) {
      return(nothing())
    }
  }
  return(just(firstelem))
}

tryRational <- function(cyc) {
  pns <- phiNrpSqfree(cyc@order)
  if(pns[["sqfree"]] && lenCyc(cyc) == pns[["phi"]]) {
    eqtrms <- equalCoefficients(cyc)
    if(is_nothing(eqtrms)) {
      cyc
    } else {
      rat <- from_just(eqtrms)
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
  ok <- powers == 1L & primes != "2"
  squareFreeOddFactors <- primes[ok]
  if(length(squareFreeOddFactors) == 0L) {
    return(cyc)
  }
  Reduce(reduceByPrime, squareFreeOddFactors, init = cyc, right = TRUE)
}

cyclotomic <- function(ord, trms) {
  cyc <- new(
    "cyclotomic",
    order = as.character(ord),
    terms = removeZeros(trms)
  )
  tryReduce(tryRational(gcdReduce(cyc)))
}

mkCyclotomic <- function(ord, trms) {
  cyclotomic(
    as.character(ord),
    convertToBase(as.character(ord), trms)
  )
}
