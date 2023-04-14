#' @title Cyclotomic as exact rational number if possible
#' @description Cyclotomic number as exact rational number if possible.
#'
#' @param cyc a cyclotomic number
#'
#' @return A \code{maybe} value, \code{just} a rational number if \code{cyc}
#'   is a rational number, \code{nothing} otherwise.
#' @export
#' @seealso \code{\link{isRational}}
#'
#' @examples
#' maybeRational(zeta(4))
#' maybeRational(cosDeg(60)) # use `from_just` to get the value
maybeRational <- function(cyc) {
  if(cyc@order == 1L) {
    trms <- cyc@terms
    if(trms$size() == 0L) {
      just(as.bigq(0L))
    } else {
      trms$at(0L)
    }
  } else {
    nothing()
  }
}

#' @title Is the cyclotomic a rational number?
#' @description Checks whether a cyclotomic number is a rational number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A Boolean value.
#' @export
#' @seealso \code{\link{maybeRational}}
isRational <- function(cyc) {
  cyc@order == 1L
}

#' @title Is the cyclotomic a Gaussian rational?
#' @description Checks whether a cyclotomic number is a Gaussian rational number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A Boolean value.
#' @export
isGaussianRational <- function(cyc) {
  isRational(realPart(cyc)) && isRational(imaginaryPart(cyc))
}

## | rational as cyclotomic ####
fromRational <- function(rat) {
  if(rat == 0L) {
    zeroCyc()
  } else {
    trms <- intmap$new()
    trms$insert(0L, rat)
    new("cyclotomic", order = 1L, terms = trms)
  }
}

## | integer as cyclotomic ####
fromInteger <- function(n) {
  fromRational(as.bigq(n))
}

