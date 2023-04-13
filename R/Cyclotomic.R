#' @importFrom methods new setClass setMethod
#' @include Cyclotomic.R
NULL

setClass(
  "cyclotomic",
  slots = c(order = "character", terms = "list")
)

setMethod(
  "show", "cyclotomic",
  function(object) {
    cat(showCyclotomic(object), "\n")
  }
)

setGeneric(
  "as.cyclotomic", function(x) {
    NULL
  }
)

#' @name as.cyclotomic
#' @aliases as.cyclotomic,character-method as.cyclotomic,cyclotomic-method as.cyclotomic,numeric-method as.cyclotomic,bigz-method as.cyclotomic,bigq-method
#' @exportMethod as.cyclotomic
#' @docType methods
#' @title Coercion to a 'cyclotomic' object
#'
#' @param x a \code{cyclotomic} object or an object yielding a quoted integer or a
#'   quoted fraction after an application of \code{as.character}
#'
#' @return A \code{cyclotomic} object.
#' @export
#'
#' @examples
#' as.cyclotomic(2)
#' as.cyclotomic("1/3")
setMethod(
  "as.cyclotomic", "character",
  function(x) {
    stopifnot(isFraction(x))
    fromRational(as.bigq(x))
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "cyclotomic",
  function(x) {
    x
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "numeric",
  function(x) {
    stopifnot(isInteger(x))
    fromInteger(x)
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "bigz",
  function(x) {
    fromRational(as.bigq(x))
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "bigq",
  function(x) {
    fromRational(x)
  }
)


#' @name cyclotomic-unary
#' @title Unary operators for cyclotomic objects
#' @description Unary operators for cyclotomic objects.
#' @aliases +,cyclotomic,missing-method -,cyclotomic,missing-method
#' @param e1 object of class \code{cyclotomic}
#' @param e2 nothing
#' @return A \code{cyclotomic} object.
setMethod(
  "+",
  signature(e1 = "cyclotomic", e2 = "missing"),
  function(e1, e2) e1
)
#' @rdname cyclotomic-unary
setMethod(
  "-",
  signature(e1 = "cyclotomic", e2 = "missing"),
  function(e1, e2) {
    minusCyc(e1)
  }
)

