#' @importFrom intmap intmap
NULL

insertWith <- function(f, mp, key, value) {
  if(mp$has_key(key)) {
    mp$insert(key, f(value, mp$get(key)), replace = TRUE)
  } else {
    mp$insert(key, value)
  }
  invisible(NULL)
}

mapKeys <- function(f, mp) {
  keys <- vapply(mp$keys(), f, FUN.VALUE = integer(1L))
  intmap$new(keys, mp$values())
}

mapValues <- function(f, mp) {
  values <- lapply(mp$values(), f)
  intmap$new(mp$keys(), values)
}

filterMap <- function(f, mp) {
  out <- intmap$new()
  for(k in mp$keys()) {
    v <- mp$get(k)
    if(f(v)) {
      out$insert(k, v)
    }
  }
  out
}

unionWith <- function(f, mp1, mp2) {
  keys1 <- mp1$keys()
  keys2 <- mp2$keys()
  commonKeys <- intersect(keys1, keys2)
  mp <- intmap$new()
  for(k in commonKeys) {
    mp$insert(k, f(mp1$get(k), mp2$get(k)))
  }
  for(k in setdiff(keys1, commonKeys)) {
    mp$insert(k, mp1$get(k))
  }
  for(k in setdiff(keys2, commonKeys)) {
    mp$insert(k, mp2$get(k))
  }
  mp
}
