#' @title Cosine and sine of a rational number
#' @description Cosine and sine of a rational angle as a cyclotomic number.
#'
#' @param theta an integer number, a \strong{gmp} rational number, or a
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
cosDeg <- function(theta) {
  stopifnot(is.bigq(theta) || isFraction(theta) || isInteger(theta))
  cosRev(as.bigq(theta) / 360L)
}

#' @rdname trigonometry
#' @export
sinDeg <- function(theta) {
  stopifnot(is.bigq(theta) || isFraction(theta) || isInteger(theta))
  sinRev(as.bigq(theta) / 360L)
}

#' @rdname trigonometry
#' @export
cosRev <- function(theta) {
  stopifnot(is.bigq(theta) || isFraction(theta) || isInteger(theta))
  n <- as.bigq(theta)
  num <- abs(as.integer(numerator(n)))
  den <- as.integer(denominator(n))
  a <- zeta(den)^num
  realPart(a)
}

#' @rdname trigonometry
#' @export
sinRev <- function(theta) {
  stopifnot(is.bigq(theta) || isFraction(theta) || isInteger(theta))
  n <- as.bigq(theta)
  num <- abs(as.integer(numerator(n)))
  den <- as.integer(denominator(n))
  a <- zeta(den)^num
  if(n >= 0L) imaginaryPart(a) else -imaginaryPart(a)
}
