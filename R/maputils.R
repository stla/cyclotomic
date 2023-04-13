#' @importFrom fastmap fastmap
#' @importFrom maybe just nothing
NULL

insertWith <- function(f, mp, key, value) {
  if(mp$has(key)) {
    mp$set(key, f(value, mp$get(key)))
  } else {
    mp$set(key, value)
  }
  invisible(NULL)
}

mapKeys <- function(f, mp) {
  out <- fastmap()
  for(k in mp$keys()) {
    out$set(f(k), mp$get(k))
  }
  out
}

mapValues <- function(f, mp) {
  out <- fastmap()
  for(k in mp$keys()) {
    out$set(k, f(mp$get(k)))
  }
  out
}

filterMap <- function(f, mp) {
  out <- fastmap()
  for(k in mp$keys()) {
    v <- mp$get(k)
    if(f(v)) {
      out$set(k, v)
    }
  }
  out
}

unionWith <- function(f, mp1, mp2) {
  keys1 <- mp1$keys()
  keys2 <- mp2$keys()
  commonKeys <- intersect(keys1, keys2)
  mp <- fastmap()
  for(k in commonKeys) {
    mp$set(k, f(mp1$get(k), mp2$get(k)))
  }
  for(k in setdiff(keys1, commonKeys)) {
    mp$set(k, mp1$get(k))
  }
  for(k in setdiff(keys2, commonKeys)) {
    mp$set(k, mp2$get(k))
  }
  mp
}

lookup <- function(key, mp) {
  if(mp$has(key)) {
    just(mp$get(key))
  } else {
    nothing()
  }
}
