#' @title Polar complex number with rational magnitude and angle
#' @description Complex number in polar form with rational magnitude and
#'   rational angle as a cyclotomic number.
#'
#' @param r magnitude, an integer number, a \strong{gmp} rational number, or a
#'   fraction given as a character string (e.g. \code{"2/7"})
#' @param theta angle, an integer number, a \strong{gmp} rational number, or a
#'   fraction given as a character string (e.g. \code{"2/7"}); for
#'   \code{polarDeg} the angle is given in degrees and for \code{polarRev}
#'   it is given in revolutions
#'
#' @return A cyclotomic number.
#' @export
#'
#' @name polar
#' @rdname polar
#'
#' @examples
#' polarDeg(1, 90)    # should be zeta(4)
#' polarRev(1, "1/4") # should be zeta(4) as well
polarDeg <- function(r, theta) {
  stopifnot(is.bigq(r) || isFraction(r) || isInteger(r))
  stopifnot(is.bigq(theta) || isFraction(theta) || isInteger(theta))
  polarRev(r, as.bigq(theta) / 360L)
}

#' @rdname polar
#' @export
polarRev <- function(r, theta) {
  stopifnot(is.bigq(r) || isFraction(r) || isInteger(r))
  stopifnot(is.bigq(theta) || isFraction(theta) || isInteger(theta))
  r     <- as.bigq(r)
  theta <- as.bigq(theta)
  p <- as.integer(numerator(theta))
  q <- as.integer(denominator(theta))
  if(p >= 0L) {
    fromRational(r) * zeta(q)^p
  } else {
    conjugate(fromRational(r) * zeta(q)^(-p))
  }
}
