## | sum of two cyclotomic numbers ####
sumCyc <- function(cyc1, cyc2) {
  o1    <- as.vli(cyc1@order)
  trms1 <- cyc1@terms
  o2    <- as.vli(cyc2@order)
  trms2 <- cyc2@terms
  ord <- lcmul(o1, o2)
  m1  <- as.integer(ord / o1)
  m2  <- as.integer(ord / o2)
  mp1 <- mapKeys(
    function(k) { as.character(m1 * as.integer(k)) }, trms1
  )
  mp2 <- mapKeys(
    function(k) { as.character(m2 * as.integer(k)) }, trms2
  )
  mkCyclotomic(ord, unionWith(`+`, mp1, mp2))
}

## | product of two cyclotomic numbers ####
prodCyc <- function(cyc1, cyc2) {
  o1    <- as.vli(cyc1@order)
  trms1 <- cyc1@terms
  o2    <- as.vli(cyc2@order)
  trms2 <- cyc2@terms
  ord <- lcmul(o1, o2)
  m1  <- ord / o1
  m2  <- ord / o2
  keys1 <- trms1$keys()
  keys2 <- trms2$keys()
  mp <- fastmap()
  for(k1 in keys1) {
    e1 <- as.vli(k1)
    c1 <- trms1$get(k1)
    for(k2 in keys2) {
      e2 <- as.vli(k2)
      c2 <- trms2$get(k2)
      k <- as.character((m1*e1 + m2*e2) %% ord)
      insertWith(`+`, mp, k, c1*c2)
    }
  }
  mkCyclotomic(ord, mp)
}

## | the zero cyclotomic number ####
zeroCyc <- function() {
  new("Cyclotomic", order = "1", terms = fastmap())
}

## | product rational and cyclotomic ####
prodRatCyc <- function(rat, cyc) {
  if(rat == 0L) {
    zeroCyc()
  } else {
    new(
      "Cyclotomic",
      order = cyc@order,
      terms = mapValues(function(x) {rat * x}, cyc@terms)
    )
  }
}

## | rational as cyclotomic ####
fromRational <- function(rat) {
  if(rat == 0L) {
    zeroCyc()
  } else {
    trms <- fastmap()
    trms$set("0", rat)
    new("Cyclotomic", order = "1", terms = trms)
  }
}

## | integer as cyclotomic ####
fromInteger <- function(n) {
  fromRational(as.bigq(n))
}

## | cyclotomic as exact rational number if possible ####
toRat <- function(cyc) {
  if(cyc@order == "1") {
    trms <- cyc@terms
    if(trms$size() == 0L) {
      just(as.bigq(0L))
    } else {
      lookup("0", trms)
    }
  } else {
    nothing()
  }
}

## helper 1 for multiplicative inverse
multiplyExponents <- function(j, cyc) { # j is vli
  n <- as.vli(cyc@order)
  if(gcd(j, n) != 1L) {
    stop("multiplyExponents needs gcd == 1")
  }
  mkCyclotomic(
    n,
    mapKeys(
      function(k) {
        as.character((j * as.vli(k)) %% n)
      },
      cyc@terms
    )
  )
}

## helper 2 for multiplicative inverse
productOfGaloisConjugates <- function(cyc) {
  ord <- as.vli(cyc@order)
  coprimes <- list() # find numbers relatively prime with ord
  j <- as.vli(2L)
  while(j < ord) {
    if(gcd(j, ord) == 1L) {
      coprimes <- c(coprimes, list(j))
    }
    j <- j + 1L
  }
  tomultiply <- lapply(coprimes, function(j) {
    multiplyExponents(j, cyc)
  })
  Reduce(prodCyc, tomultiply)
}

## | multiplicative inverse ####
invCyc <- function(cyc) {
  pgc <- productOfGaloisConjugates(cyc)
  maybe_rat <- toRat(prodCyc(cyc, pgc))
  if(is_just(maybe_rat)) {
    r <- from_just(maybe_rat)
    prodRatCyc(1L/r, pgc)
  } else {
    stop("this is a bug")
  }
}
