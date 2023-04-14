#' @title Cosine and sine of a rational number
#' @description Cosine and sine of a rational number as a cyclotomic number.
#'
#' @param rat an integer number, a \strong{gmp} rational number, or a
#'   fraction given as a character string (e.g. \code{"2/7"})
#'
#' @return A cyclotomic number.
#' @export
#' @name trigonometry
#' @rdname trigonometry
#'
#' @details The function \code{cosDeg}, resp. \code{sinDeg}, returns the cosine,
#'   resp. the sine, of its argument assumed to be given in degrees.
#'   The function \code{cosRev}, resp. \code{sinRev}, returns the cosine,
#'   resp. the sine, of its argument assumed to be given in revolutions.
#'
#' @examples
#' cosDeg(60)
#' cosDeg("2/3")^2 + sinDeg("2/3")^2 == 1
cosDeg <- function(rat) {
  stopifnot(is.bigq(rat) || isFraction(rat) || isInteger(rat))
  cosRev(as.bigq(rat) / 360L)
}

#' @rdname trigonometry
#' @export
sinDeg <- function(rat) {
  stopifnot(is.bigq(rat) || isFraction(rat) || isInteger(rat))
  sinRev(as.bigq(rat) / 360L)
}

#' @rdname trigonometry
#' @export
cosRev <- function(rat) {
  stopifnot(is.bigq(rat) || isFraction(rat) || isInteger(rat))
  n <- as.bigq(rat)
  num <- abs(as.integer(numerator(n)))
  den <- as.integer(denominator(n))
  a <- zeta(den)^num
  realPart(a)
}

#' @rdname trigonometry
#' @export
sinRev <- function(rat) {
  stopifnot(is.bigq(rat) || isFraction(rat) || isInteger(rat))
  n <- as.bigq(rat)
  num <- abs(as.integer(numerator(n)))
  den <- as.integer(denominator(n))
  a <- zeta(den)^num
  out <- imaginaryPart(a)
  if(n < 0L) {
    out <- -out
  }
  out
}
