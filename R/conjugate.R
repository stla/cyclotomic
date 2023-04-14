#' @title Conjugate cyclotomic number
#' @description Conjugate of a cyclotomic number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A cyclotomic number, the conjugate of \code{cyc}.
#' @export
#'
#' @examples
#' conjugate(e(4)) # should be -e(4)
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
#' realPart(e(9))
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
#' imaginaryPart(e(9))
imaginaryPart <- function(cyc) {
  - as.bigq("1/2") * e(4) * (cyc - conjugate(cyc))
}
