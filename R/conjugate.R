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
  ord <- cyc@order
  n <- as.vli(ord)
  trms <- cyc@terms
  conjtrms <- mapKeys(
    function(k) {
      as.character((n - as.vli(k)) %% n)
    },
    trms
  )
  mkCyclotomic(ord, conjtrms)
}
