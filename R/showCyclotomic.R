showCyclotomic <- function(cyc) {
  n  <- cyc@order
  mp <- cyc@terms
  if(mp$size() == 0L) {
    return("0")
  }
  xxs <- mp$toList()
  ex  <- names(xxs)[1L]
  rat <- xxs[[1L]]
  xs  <- xxs[-1L]
  paste0(
    leadingTerm(rat, n, ex), followingTerms(n, xs)
  )
}

showBaseExp <- function(n, ex) {
  en <- sprintf("zeta(%s)", n)
  if(ex == 1L) {
    en
  } else {
    sprintf("%s^%s", en, ex)
  }
}

leadingTerm <- function(rat, n, ex) {
  if(ex == 0L) {
    as.character(rat)
  } else {
    t <- showBaseExp(n, ex)
    if(rat == 1L) {
      t
    } else if(rat == -1L) {
      paste0("-", t)
    } else if(rat != 0L) {
      paste0(as.character(rat), "*", t)
    } else {
      ""
    }
  }
}

followingTerms <- function(n, lst) {
  if(length(lst) == 0L) {
    ""
  } else {
    ex <- names(lst)[1L]
    rat <- lst[[1L]]
    paste0(
      followingTerm(rat, n, ex),
      followingTerms(n, lst[-1L])
    )
  }
}

followingTerm <- function(rat, n, ex) {
  t <- showBaseExp(n, ex)
  if(rat == 1L) {
    paste0(" + ", t)
  } else if(rat == -1L) {
    paste0(" - ", t)
  } else if(rat > 0L) {
    paste0(" + ", as.character(rat), "*", t)
  } else if(rat < 0L) {
    paste0(" + ", as.character(-rat), "*", t)
  } else {
    ""
  }
}
