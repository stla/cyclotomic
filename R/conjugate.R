#' @title Conjugate cyclotomic number
#' @description Complex conjugate of a cyclotomic number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A cyclotomic number, the complex conjugate of \code{cyc}.
#' @export
#'
#' @examples
#' conjugate(zeta(4)) # should be -zeta(4)
conjugate <- function(cyc) {
  n    <- cyc@order
  trms <- cyc@terms
  conjtrms <- mapKeys(
    function(k) {
      (n - k) %% n
    },
    trms
  )
  mkCyclotomic(n, conjtrms)
}

#' @title Real part of cyclotomic number.
#' @description The real part of a cyclotomic number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A cyclotomic number.
#' @export
#'
#' @examples
#' realPart(zeta(9))
realPart <- function(cyc) {
  as.bigq("1/2") * (cyc + conjugate(cyc))
}

#' @title Imaginary part of cyclotomic number.
#' @description The imaginary part of a cyclotomic number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A cyclotomic number.
#' @export
#'
#' @examples
#' imaginaryPart(zeta(9))
imaginaryPart <- function(cyc) {
  - as.bigq("1/2") * zeta(4) * (cyc - conjugate(cyc))
}

#' @title Is the cyclotomic a real number?
#' @description Checks whether a cyclotomic number is a real number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A Boolean value.
#' @export
isReal <- function(cyc) {
  cyc == conjugate(cyc)
}
